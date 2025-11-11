module Vicsek

using LinearAlgebra
using Random
using Plots

import Base: copy

export VicsekSimulation, step!, run!, positions, velocities, visualize

mutable struct NeighborList{T}
    heads::Vector{Int}
    lscl::Vector{Int}
    ncell::NTuple{2, Int}
    cellsize::T
    inv_cellsize::T
    box_length::T
end

NeighborList(box_length::T, cutoff::T, N::Int) where {T<:AbstractFloat} = begin
    box_length <= zero(T) && throw(ArgumentError("box_length must be positive"))
    cutoff <= zero(T) && throw(ArgumentError("cutoff must be positive"))
    cells_per_side = max(1, Int(floor(box_length / cutoff)))
    cellsize = cells_per_side > 0 ? box_length / cells_per_side : box_length
    heads = fill(0, cells_per_side * cells_per_side)
    lscl = zeros(Int, N)
    NeighborList{T}(heads, lscl, (cells_per_side, cells_per_side), cellsize, one(T) / cellsize, box_length)
end

@inline function cell_index(nl::NeighborList{T}, x::T, y::T) where {T}
    ncellx, ncelly = nl.ncell
    cx = Int(floor(x * nl.inv_cellsize))
    cy = Int(floor(y * nl.inv_cellsize))
    cx = ifelse(cx >= ncellx, ncellx - 1, ifelse(cx < 0, 0, cx))
    cy = ifelse(cy >= ncelly, ncelly - 1, ifelse(cy < 0, 0, cy))
    return cx + cy * ncellx + 1
end

function rebuild!(nl::NeighborList{T}, positions::Matrix{T}) where {T}
    fill!(nl.heads, 0)
    N = size(positions, 2)
    if length(nl.lscl) != N
        resize!(nl.lscl, N)
    end
    @inbounds for i in 1:N
        cell = cell_index(nl, positions[1, i], positions[2, i])
        nl.lscl[i] = nl.heads[cell]
        nl.heads[cell] = i
    end
    return nothing
end

mutable struct VicsekSimulation{T, RNGType<:AbstractRNG}
    positions::Matrix{T}
    velocities::Matrix{T}
    buffer_velocities::Matrix{T}
    speed::T
    dt::T
    alignment_radius2::T
    lj_sigma::T
    lj_sigma2::T
    lj_sigma6::T
    lj_epsilon::T
    lj_cutoff2::T
    alignment_strength::T
    force_factor::T
    noise_amplitude::T
    box_length::T
    neighbor_cutoff2::T
    neighborlist::NeighborList{T}
    rng::RNGType
end

function VicsekSimulation(N::Integer; kwargs...)
    return VicsekSimulation{Float64}(N; kwargs...)
end

function VicsekSimulation{T}(N::Integer; box_length::Real=10.0, speed::Real=0.03, dt::Real=1.0,
        alignment_radius::Real=0.5, noise_amplitude::Real=0.2,
        alignment_strength::Real=1.0, force_factor::Real=0.03,
        lj_sigma::Real=0.1, lj_epsilon::Real=1.0,
        lj_cutoff::Real=0.25 * lj_sigma, rng::AbstractRNG=Random.default_rng()) where {T<:AbstractFloat}
    N <= 0 && throw(ArgumentError("N must be positive"))
    bl = T(box_length)
    spd = T(speed)
    timestep = T(dt)
    align_r = T(alignment_radius)
    noise = T(noise_amplitude)
    align_strength = T(alignment_strength)
    force_fac = T(force_factor)
    sigma = T(lj_sigma)
    epsilon = T(lj_epsilon)
    cutoff = T(lj_cutoff)
    bl <= zero(T) && throw(ArgumentError("box_length must be positive"))
    spd <= zero(T) && throw(ArgumentError("speed must be positive"))
    timestep <= zero(T) && throw(ArgumentError("dt must be positive"))
    align_r <= zero(T) && throw(ArgumentError("alignment_radius must be positive"))
    noise < zero(T) && throw(ArgumentError("noise_amplitude must be non-negative"))
    sigma <= zero(T) && throw(ArgumentError("lj_sigma must be positive"))
    epsilon < zero(T) && throw(ArgumentError("lj_epsilon must be positive"))
    cutoff <= zero(T) && throw(ArgumentError("lj_cutoff must be positive"))
    neighbor_cutoff = max(align_r, cutoff)
    positions = rand(rng, T, 2, N) .* bl
    velocities = zeros(T, 2, N)
    buffer_velocities = similar(velocities)
    @inbounds for i in 1:N
        θ = T(2π) * rand(rng, T)
        velocities[1, i] = spd * cos(θ)
        velocities[2, i] = spd * sin(θ)
    end
    # remove COM motion
    vcmx = sum(view(velocities, 1, :)) / T(N)
    vcmy = sum(view(velocities, 2, :)) / T(N)
    @inbounds for i in 1:N
        velocities[1, i] -= vcmx
        velocities[2, i] -= vcmy
    end
    # normaize velocities to exact speed
    @inbounds for i in 1:N
        vx = velocities[1, i]
        vy = velocities[2, i]
        mag = sqrt(vx * vx + vy * vy)
        invmag = spd / mag
        velocities[1, i] = vx * invmag
        velocities[2, i] = vy * invmag
    end

    nl = NeighborList(T(bl), T(neighbor_cutoff), N)
    rebuild!(nl, positions)
    sim = VicsekSimulation{T, typeof(rng)}(positions, velocities, buffer_velocities,
        spd, timestep, align_r^2, sigma, sigma^2, sigma^6, epsilon,
        cutoff^2, align_strength, force_fac, noise, bl, neighbor_cutoff^2, nl, rng)
    return sim
end

@inline function wrap_coord(x::T, L::T) where {T}
    y = mod(x, L)
    return y == L ? zero(T) : y
end

@inline function min_image(dx::T, L::T) where {T}
    return dx - L * round(T, dx / L)
end

function step!(sim::VicsekSimulation{T, RNGType}, steps::Integer=1) where {T, RNGType}
    steps < 0 && throw(ArgumentError("steps must be non-negative"))
    steps == 0 && return sim
    @inbounds for _ in 1:steps
        rebuild!(sim.neighborlist, sim.positions)
        update_velocities!(sim)
        update_positions!(sim)
    end
    return sim
end

function update_velocities!(sim::VicsekSimulation{T, RNGType}) where {T, RNGType}
    pos = sim.positions
    vel = sim.velocities
    next_vel = sim.buffer_velocities
    rng = sim.rng
    nl = sim.neighborlist
    N = size(pos, 2)
    ncellx, ncelly = nl.ncell
    heads = nl.heads
    lscl = nl.lscl
    box = sim.box_length
    align_r2 = sim.alignment_radius2
    lj_cut2 = sim.lj_cutoff2
    neighbor_cut2 = sim.neighbor_cutoff2
    sigma2 = sim.lj_sigma2
    sigma6 = sim.lj_sigma6
    epsilon = sim.lj_epsilon
    align_strength = sim.alignment_strength
    force_factor = sim.force_factor
    speed = sim.speed
    noise = sim.noise_amplitude
    inv_cellsize = nl.inv_cellsize
    for i in 1:N
        xi = pos[1, i]
        yi = pos[2, i]
        cx = Int(floor(xi * inv_cellsize))
        cy = Int(floor(yi * inv_cellsize))
        cx = ifelse(cx >= ncellx, ncellx - 1, ifelse(cx < 0, 0, cx))
        cy = ifelse(cy >= ncelly, ncelly - 1, ifelse(cy < 0, 0, cy))
        ax = vel[1, i]
        ay = vel[2, i]
        fx = zero(T)
        fy = zero(T)
        for dy in -1:1
            ny = cy + dy
            ny = ny < 0 ? ncelly - 1 : ny
            ny = ny >= ncelly ? 0 : ny
            for dx in -1:1
                nx = cx + dx
                nx = nx < 0 ? ncellx - 1 : nx
                nx = nx >= ncellx ? 0 : nx
                cell = nx + ny * ncellx + 1
                j = heads[cell]
                while j != 0
                    if j != i
                        dxp = min_image(pos[1, j] - xi, box)
                        dyp = min_image(pos[2, j] - yi, box)
                        r2 = dxp * dxp + dyp * dyp
                        if r2 <= neighbor_cut2 && r2 > zero(T)
                            if r2 <= align_r2
                                ax += vel[1, j]
                                ay += vel[2, j]
                            end
                            if r2 <= lj_cut2
                                inv_r2 = inv(r2)
                                sr2 = sigma2 * inv_r2
                                sr6 = sr2 * sr2 * sr2
                                sr12 = sr6 * sr6
                                fscalar = T(24) * epsilon * (T(2) * sr12 - sr6) * inv_r2
                                # Add Lennard-Jones force contribution (correct sign)
                                fx += fscalar * dxp
                                fy += fscalar * dyp
                            end
                        end
                    end
                    j = lscl[j]
                end
            end
        end
        ax *= align_strength
        ay *= align_strength
        desired_x = ax + force_factor * fx 
        desired_y = ay + force_factor * fy 
        mag2 = desired_x * desired_x + desired_y * desired_y
        if mag2 == zero(T)
            desired_x = vel[1, i]
            desired_y = vel[2, i]
            mag2 = desired_x * desired_x + desired_y * desired_y
        end
        invmag = inv(sqrt(mag2))
        desired_x *= invmag
        desired_y *= invmag
        θ = atan(desired_y, desired_x)  + noise * (rand(rng, T) - T(0.5))
        next_vel[1, i] = speed * cos(θ) 
        next_vel[2, i] = speed * sin(θ)
    end
    copyto!(vel, next_vel)
    return nothing
end

function update_positions!(sim::VicsekSimulation{T, RNGType}) where {T, RNGType}
    pos = sim.positions
    vel = sim.velocities
    box = sim.box_length
    dt = sim.dt
    N = size(pos, 2)
    @inbounds for i in 1:N
        pos[1, i] = wrap_coord(pos[1, i] + vel[1, i] * dt, box)
        pos[2, i] = wrap_coord(pos[2, i] + vel[2, i] * dt, box)
    end
    return nothing
end

function run!(sim::VicsekSimulation, steps::Integer)
    step!(sim, steps)
    return sim
end

@inline positions(sim::VicsekSimulation) = sim.positions
@inline velocities(sim::VicsekSimulation) = sim.velocities

function copy(sim::VicsekSimulation{T, RNGType}) where {T, RNGType}
    pos = copy(sim.positions)
    vel = copy(sim.velocities)
    buffer = copy(sim.buffer_velocities)
    rng_copy = deepcopy(sim.rng)
    nl = NeighborList(T(sim.box_length), sqrt(sim.neighbor_cutoff2), size(pos, 2))
    rebuild!(nl, pos)
    return VicsekSimulation{T, typeof(rng_copy)}(pos, vel, buffer, sim.speed, sim.dt,
        sim.alignment_radius2, sim.lj_sigma, sim.lj_sigma2, sim.lj_sigma6,
        sim.lj_epsilon, sim.lj_cutoff2, sim.alignment_strength, sim.force_factor,
        sim.noise_amplitude, sim.box_length, sim.neighbor_cutoff2, nl, rng_copy)
end

function visualize(sim::VicsekSimulation; frames::Integer=240, stride::Integer=1,
        keep_state::Bool=true, show_velocity::Bool=true, markersize::Real=0.1,
        particle_color=:navy, background=:white, framerate::Integer=30,
        savepath::Union{Nothing, AbstractString}=nothing, title::Union{Nothing, AbstractString}=nothing)
    frames <= 0 && throw(ArgumentError("frames must be positive"))
    stride <= 0 && throw(ArgumentError("stride must be positive"))
    local_sim = keep_state ? copy(sim) : sim
    box = local_sim.box_length
    ttl = title === nothing ? "Vicsek model" : String(title)
    anim = @animate for frame in 1:frames
        if frame > 1
            step!(local_sim, stride)
        end
        pos = local_sim.positions
        vel = local_sim.velocities
        plt = scatter(view(pos, 1, :), view(pos, 2, :);
            xlim=(0, box), ylim=(0, box), aspect_ratio=1,
            color=particle_color, markersize=markersize, markerstrokewidth=0,
            legend=false, background_color=background, title="$(ttl) – frame $(frame)",
            xlabel="x", ylabel="y")
        if show_velocity
            # Compute scaled velocity arrows without mutating simulation state
            velx_v = view(vel, 1, :)
            vely_v = view(vel, 2, :)
            speed = sqrt.(velx_v .^ 2 .+ vely_v .^ 2)
            arrow_x = velx_v .* (0.1 ./ speed)
            arrow_y = vely_v .* (0.1 ./ speed)
            quiver!(view(pos, 1, :), view(pos, 2, :),
                quiver=(arrow_x, arrow_y); color=:gray70, lw=0.75)
        end
        plt
    end
    if savepath === nothing
        return anim
    else
        gif(anim, savepath; fps=framerate)
        return savepath
    end
end

end # module

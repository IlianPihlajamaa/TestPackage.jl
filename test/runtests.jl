using TestPackage
using Test, Random, Base.Threads


const _skip_heavy_flag = any(arg -> arg == "--skip-heavy", ARGS)
const _run_heavy_env = lowercase(get(ENV, "RUN_HEAVY_TESTS", "true"))
const _run_heavy = !_skip_heavy_flag && !(_run_heavy_env in ("false", "0", "no"))

@testset "TestPackage.jl" begin
    @testset "Vicsek simulation basics" begin
        rng = Random.MersenneTwister(1234)
        sim = VicsekSimulation(64; box_length=32.0, alignment_radius=1.5, noise_amplitude=0.05, rng=rng)
        @test size(positions(sim)) == (2, 64)
        @test size(velocities(sim)) == (2, 64)
        step!(sim, 10)
        pos = positions(sim)
        box = sim.box_length
        @test all(@view(pos[1, :]) .>= 0) && all(@view(pos[1, :]) .< box)
        @test all(@view(pos[2, :]) .>= 0) && all(@view(pos[2, :]) .< box)
        vels = velocities(sim)
        for j in 1:size(vels, 2)
            speed = hypot(vels[1, j], vels[2, j])
            @test isapprox(speed, sim.speed; atol=1e-10)
        end
    end
    @testset "Visualize does not mutate velocities" begin
        rng = Random.MersenneTwister(42)
        sim = VicsekSimulation(16; box_length=8.0, alignment_radius=1.0, noise_amplitude=0.0, rng=rng)
        # one frame, no stepping inside visualize when frame == 1
        visualize(sim; frames=1, stride=1, keep_state=false, show_velocity=true, markersize=0.5)
        vels = velocities(sim)
        for j in 1:size(vels, 2)
            @test isapprox(hypot(vels[1, j], vels[2, j]), sim.speed; atol=1e-10)
        end
    end
   
end

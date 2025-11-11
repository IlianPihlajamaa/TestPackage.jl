module TestPackage

include("vicsek.jl")

using .Vicsek

export VicsekSimulation, step!, run!, positions, velocities, visualize, Vicsek

end # module TestPackage

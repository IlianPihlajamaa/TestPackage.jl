using TestPackage
using Random

rng = Random.MersenneTwister(543)
sim = VicsekSimulation(3000; lj_epsilon=1.0)

# Relax the system a bit before visualizing
# @time step!(sim, 1000)

# Generate an animation (writes vicsek.gif by default)
@time visualize(sim; frames=100, stride=1, savepath="vicsek.gif", markersize=1,
    particle_color=:dodgerblue)

# TestPackage

A Julia package that experiments with a 2D Vicsek model augmented with Lennard-Jones interactions and a fast linked-cell neighbor list. The package also contains a deliberately heavy stress-test that produces and sums one billion floating-point numbers.

## Features

- Fast 2D Vicsek model simulation with periodic boundaries
- Lennard-Jones forces coupled into the Vicsek alignment update
- Linked-cell neighbor list for \(\mathcal{O}(N)\) neighbor discovery
- Simple visualization helper built on `Plots.jl`
- Example script under `examples/vicsek_demo.jl`

## Getting started

Activate and instantiate the project (Julia 1.8 or newer recommended):

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Run the lightweight Vicsek regression test (skipping the billion-float stress case):

```julia
using Pkg
Pkg.test("TestPackage"; test_args=["--skip-heavy"])
```

> **Note**
> The default `runtests.jl` also contains a long-running stress test that fabricates one billion floats across 100 threads. This can take many minutes and allocate hundreds of gigabytes of memory. Passing `--skip-heavy` (as above) or setting `RUN_HEAVY_TESTS=false` will bypass it. Leave the flag unset or set `RUN_HEAVY_TESTS=true` to run the stress test.

### Simulation usage

```julia
using TestPackage

sim = VicsekSimulation(512; box_length=64.0, alignment_radius=1.5, noise_amplitude=0.2)
step!(sim, 1000)

# Animate (saves `vicsek.gif`)
visualize(sim; frames=240, stride=2, savepath="vicsek.gif")
```

### Example script

```julia
julia --project examples/vicsek_demo.jl
```

This will save an animated GIF in the working directory.

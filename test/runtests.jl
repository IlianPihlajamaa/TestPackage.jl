using TestPackage
using Test, Random
@show Random.GLOBAL_SEED

function foo(rng, n)
    println("Starting to process $n floats...")
    buffer = 0
    @time for i in 1:n
        # Create a random float
        x = rand(rng) * 100.0
        # Round it and convert to Int128
        rounded_val = ceil(Int, x)
        # Add to buffer
        buffer += rounded_val
    end
    return buffer
end


@testset "TestPackage.jl" begin
    # Create a billion floats in a loop, round them, and add to buffer
    @testset "Billion floats test" begin
        n = 1_000_000_000
        seed = rand(UInt)
        @show seed
        rng = MersenneTwister(seed)
        buffer = foo(rng, n)
        println("Final buffer value: $buffer")
        
        # Basic sanity check - buffer should be positive and non-zero
        @test buffer > 0
        @test buffer isa Int
    end
end

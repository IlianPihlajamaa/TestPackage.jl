using TestPackage
using Test, Random

function foo(rng, n, M)
    println("Starting to process $(n*M) floats...")
    buffer = 0
    cache = zeros(n)
    @time for j in 1:M
        rand!(rng, cache)
        cache .*= 10.0
        for i in eachindex(cache)
            # Create a random float
            x = cache[i]
            # Round it and convert to Int128
            rounded_val = ceil(Int, x)
            # Add to buffer
            buffer += rounded_val
        end
    end
    return buffer
end


@testset "TestPackage.jl" begin
    # Create a billion floats in a loop, round them, and add to buffer
    @testset "Billion floats test" begin
        n = 1000
        M = 1000000
        seed = rand(UInt)
        @show seed
        Random.seed!(seed)
        rng = Random.GLOBAL_RNG
        buffer = foo(rng, n, M)
        println("Final buffer value: $buffer")
        
        # Basic sanity check - buffer should be positive and non-zero
        @test buffer > 0
        @test buffer isa Int
    end
end

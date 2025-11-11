using TestPackage
using Test, Random, Base.Threads

function foo(n, M)
    println("Starting to process $(n*M) floats...")
    buffers = zeros(Int, M)
    seed = rand(UInt)
    @show seed
    @threads for j in 1:M
        rng = Random.MersenneTwister(j+seed)
        buffer = 0
        for i in 1:n
            # Create a random float
            x = rand(rng) * 10.0
            # Round it and convert to Int
            rounded_val = ceil(Int, x)
            # Add to buffer
            buffer += rounded_val
        end
        buffers[j] = buffer
    end
    return sum(buffers)
end


@testset "TestPackage.jl" begin
    # Create a billion floats in a loop, round them, and add to buffer
    @testset "Billion floats test" begin
        n = 100000000
        M = 100
        buffer = foo(n, M)
        println("Final buffer value: $buffer")
        
        # Basic sanity check - buffer should be positive and non-zero
        @test buffer > 0
        @test buffer isa Int
    end
end

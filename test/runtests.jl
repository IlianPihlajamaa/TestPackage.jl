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

        x = 4.0
        x1 = prevfloat(x)
        x11 = prevfloat(x1)
        x2 = nextfloat(x)
        x22 = nextfloat(x2)
        c = ceil(Int, x)
        c1 = ceil(Int, x1)
        c11 = ceil(Int, x11)
        c2 = ceil(Int, x2)
        c22 = ceil(Int, x22)
        @test c == 4
        @test c1 == 4
        @test c2 == 5
        @test c11 == 4
        @test c22 == 5
    end
end

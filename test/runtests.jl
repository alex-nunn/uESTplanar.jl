using uESTplanar
using Test

import uESTplanar: I, I_fast, ∂ρI_fast, ∂zI_fast, ∂ρI, ∂zI

@testset "uESTplanar.jl" begin
    @testset "kernels" begin
        @testset "I(ρ, z)" begin
            @test I(1.0, 0.5) ≈ 15.590913944908 / 4
            @test I(1.0, 0.01) ≈ 8.57496164339392 / 4
            @test I(1.0, 0.001) ≈ 8.3869774521349 / 4
            @test I(0.01, 0.5) ≈ 107.99360818005 / 4
            @test I(1000, 0.1) ≈ 0.012566367619371 / 4
        end
        
        @testset "K(ρ, z)" begin 
            @test K(0.00, 0.51) ≈ 0.55707347701202
            @test K(1.00, 0.01) ≈ 0.0010518318237137 atol=1e-6 rtol=1e-6
        end

        @testset "I_fast(ρ, z)" begin
            ρs = range(0, 10, 100)
            z = 1.0
            @test I_fast.(ρs, z) ≈ I.(ρs, z) atol=1e-3 rtol=1e-3
        end

        @testset "∂ρI_fast(ρ, z)" begin
            ρs = range(0, 10, 100)
            z = 1.0
            @test ∂ρI_fast.(ρs, z) ≈ ∂ρI.(ρs, z) atol=1e-3 rtol=1e-3
        end

        @testset "∂zI_fast(ρ, z)" begin
            ρs = range(0, 10, 100)
            z = 1.0
            @test ∂zI_fast.(ρs, z) ≈ ∂zI.(ρs, z) atol=1e-3 rtol=1e-3
        end
    end

    @testset "convolutions" begin
        @testset "round_pow2" begin
            @test round_pow2(2) == 2
            @test round_pow2(3) == 2
            @test round_pow2(4) == 4
            @test round_pow2(5) == 4
            @test round_pow2(7) == 8
        end
    end

end

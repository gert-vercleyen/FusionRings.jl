using Test
using FusionRings

@testset "Demo zoo" begin
    F = fibonacci_ring()
    @test tensor_product(F, :τ, :τ) == Dict(Symbol("1")=>1, :τ=>1)
    I = ising_ring()
    @test tensor_product(I, :σ, :σ) == Dict(:1=>1, :ψ=>1)
    Z5 = zn_fusion_ring(5)
    @test is_commutative(Z5)
end

@testset "SU(2)_k examples" begin
    S = su2_fusion_ring(4)
    @test is_commutative(S)
    @test length(labels(S)) == 5
end

@testset "TY & near-group" begin
    TY = ty_fusion_ring(["0","1","2"])
    @test tensor_product(TY, Symbol("m"), Symbol("m"))[Symbol("0")] == 1
    NG = near_group_ring(["0","1"]; n=1)
    @test tensor_product(NG, Symbol("X"), Symbol("X"))[Symbol("0")] == 1
end

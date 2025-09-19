
using Test
using FusionRings

@testset "Demo zoo" begin
    F = fibonacci_ring()
    @test tensor_product(F, :τ, :τ) == Dict(Symbol("1")=>1, :τ=>1)
    I = ising_ring()
    @test any(x->x==(Symbol("1"),1), decompose(I, :σ, :σ))
end

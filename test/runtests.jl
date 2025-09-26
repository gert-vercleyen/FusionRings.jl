
using Test
using FusionRings

@testset "Demo zoo" begin
    F = fibonacci_ring()
    # tensor_product should accept symbols (converted) and return Dict keyed by String
    @test tensor_product(F, :τ, :τ) == Dict("1"=>1, "τ"=>1)
    @test tensor_product(F, "τ", "τ") == Dict("1"=>1, "τ"=>1)
    I = ising_ring()
    dec = decompose(I, :σ, :σ)
    @test any(x->x==("1",1), dec)
    @test any(x->x==("ψ",1) || x==("σ",1), dec)  # sanity: at least one other term appears
end

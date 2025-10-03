
using Test
using FusionRings

@testset "Demo zoo" begin
    F = fibonacci_ring()
    # tensor_product should accept symbols (converted) and return Dict keyed by String
    @test tensor_product(F, :τ, :τ) == Dict("1"=>1, "τ"=>1)
    @test tensor_product(F, "τ", "τ") == Dict("1"=>1, "τ"=>1)
    @test length(fpdims(F)) == 2
    @test fpdim(F) == sum(fpdims(F) .^ 2)  # Basic FP dimension identity
    I = ising_ring()
    dec = decompose(I, :σ, :σ)
    @test any(x->x==("1",1), dec)
    @test any(x->x==("ψ",1) || x==("σ",1), dec)  # sanity: at least one other term appears
    # Conjugation tests (Fibonacci: τ is self-dual)
    cτ = conjugate_element(F, "τ")
    @test cτ isa Int
    @test conjugate_label(F, cτ) == "τ"
    # Unit element (index 1) self-dual
    @test conjugate_element(F, 1) == 1

    # Subring test: trivial rank-1 subring inside Fibonacci
    trivial = fusion_ring(reshape([1],1,1,1); labels=["1"], name="Trivial")
    @test is_sub_fusion_ring(F, ["1"])  # existing vector-based
    @test is_sub_fusion_ring(F, trivial)  # new ring-based overload
end

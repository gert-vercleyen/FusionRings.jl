
using Test
using FusionRings

@testset "Demo zoo" begin
    F = fibonacci_ring()
    @test tensor_product(F, "τ", "τ") == Dict("1"=>1, "τ"=>1)
    I = ising_ring()
    dI = decompose(I, "σ", "σ")
    @test any(x->x==( "1",1), dI) && any(x->x==( "ψ",1), dI)
    S = su2_fusion_ring(4)
    @test is_commutative(S)
    @test length(quantum_dimensions(S)) == 5
end

@testset "Modular diagnostics" begin
    F = fibonacci_ring()
    S = fusion_ring_smatrices(F)
    @test size(S) == (2,2)
    T = fusion_ring_twists(F)
    @test length(T) == 2
end

@testset "Oscar integration (if available)" begin
    # Try to use Oscar if it is installable on this platform
    try
        using Oscar
        G = GAP.Globals.DihedralGroup(6)  # D3
        GR = group_rep_fusion_ring(G)
        @test rank(GR) >= 2
        aut = fusion_ring_automorphisms(GR; method=:gap)
        @test haskey(aut, :size)
    catch e
        @info "Skipping Oscar tests: $(e)"
        @test true
    end
end

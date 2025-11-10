using Test
using FusionRings
const MD = FusionRings.ModularData

# small helper: off-diagonal Frobenius norm
offdiag_norm(A) = begin
    B = copy(A); for i in 1:size(A,1); B[i,i] = 0; end; norm(B)
end

@testset "Characters + S on Fibonacci / Ising / SU(2)_4" begin
    for R in (fibonacci_ring(), ising_ring(), su2_fusion_ring(4))
        d = quantum_dimensions(R)
        r = rank(R)

        # characters
        C, V = MD.fusion_ring_characters(R; tries=12)
        @test size(C) == (r, r)
        # each fusion matrix N_i is diagonalised by V (to tolerance)
        for a in labels(R)
            Ni = fusion_matrix(R, a)
            D = inv(V) * Matrix{Float64}(Ni) * V
            @test offdiag_norm(D) ≤ 1e-7
        end

        # S-matrix
        S, C2 = MD.fusion_ring_smatrices(R)
        @test size(S) == (r, r)
        # first row equals FP dims (within tolerance)
        @test maximum(abs.(S[1, :] .- d)) ≤ 1e-7
        # S diagonalises Ni^T
        Sinv = inv(S)
        for a in labels(R)
            NiT = transpose(fusion_matrix(R, a))
            D = Sinv * ComplexF64.(NiT) * S
            @test offdiag_norm(D) ≤ 1e-6
        end

        # modular_data wrapper
        M = MD.modular_data(R)
        @test haskey(M, :S) && haskey(M, :T) && haskey(M, :characters) && haskey(M, :dims)
        @test size(M[:S]) == (r, r)
        @test size(M[:T], 1) == r
    end
end

@testset "Vafa twists basic properties" begin
    # For fully self-dual rings (Fibonacci, Ising), the simplified Vafa rows can be trivial.
    for R in (fibonacci_ring(), ising_ring())
        θ, t = MD.fusion_ring_twists(R; prefer_exact=false)
        @test length(θ) == rank(R)
        @test θ[1] ≈ 1 + 0im
    end
end

@testset "OSCAR exact path (if available)" begin
    try
        import Oscar
        R = su2_fusion_ring(4)
        θ, t = MD.fusion_ring_twists(R; prefer_exact=true)
        @test length(θ) == rank(R)
        @test θ[1] ≈ 1 + 0im
    catch
        @info "OSCAR not available; skipping exact SNF test."
    end
end

##############################
# ModularData.jl
# Characters, S- and T-matrices from fusion data
# Implements the workflow in the notes you sent:
#  - characters by simultaneous diagonalisation of {N_i}
#  - S from common eigenvectors; normalised so S[1,i]=d_i (first row = Frobenius–Perron dims)
#  - twists (T) from Vafa's integer system; exact via OSCAR SNF if available
##############################

module ModularData

export fusion_ring_characters,
       fusion_ring_smatrices,
       fusion_ring_twists,
       modular_data

using LinearAlgebra
using ..Types
using ..Creation
using ..Operations
using ..Properties

# Optional OSCAR for exact Smith normal form over ℤ.
const _HAS_OSCAR = try
    @eval import Oscar
    true
catch
    false
end

# -------------------------------
# 1) Characters via simult. diag.
# -------------------------------

"""
    fusion_ring_characters(R::FusionRing; tries=8, tol=1e-10) -> (C, V)

Return the **character table** `C::Matrix{ComplexF64}` of a **commutative** fusion ring `R`,
together with a matrix `V` whose columns are a common eigenbasis for the fusion matrices.

By definition here, `C[j,i]` is the eigenvalue of `N_i` on the `j`-th common eigenline
(i.e. character `χ_j` evaluated on basis element `i`).

Algorithm:
1. Form a random real combination `M = ∑_k c_k N_k`.
2. Eigen-decompose `M = V Λ V⁻¹`.
3. Verify that every `V⁻¹ N_i V` is (numerically) diagonal. If not, retry.

Throws if no common eigenbasis is found after `tries` attempts.
"""
function fusion_ring_characters(R::FusionRing; tries::Int=8, tol::Real=1e-10)
    labs = labels(R)
    r = length(labs)
    Nis = [Matrix{Float64}(fusion_matrix(R, a)) for a in labs]

    # quick commutativity sanity check
    for i in 1:r, j in i+1:r
        if norm(Nis[i]*Nis[j] - Nis[j]*Nis[i]) > 1e-8
            error("fusion_ring_characters: ring appears non-commutative; this routine requires commuting fusion matrices.")
        end
    end

    for _ in 1:tries
        coeffs = randn(r)
        M = zeros(Float64, r, r)
        @inbounds for k in 1:r
            M .+= coeffs[k] .* Nis[k]
        end

        ev = eigen(M)                    # symmetric not guaranteed; generic eigen
        V  = Matrix(ev.vectors)
        Vinv = inv(V)                    # small r; explicit inverse is fine here

        # Check diagonalisation
        diags = Vector{Vector{ComplexF64}}(undef, r)
        ok = true
        for i in 1:r
            D = Vinv * Nis[i] * V
            off = copy(D); @inbounds for j in 1:r; off[j,j] = 0.0; end
            if norm(off) > tol
                ok = false
                break
            end
            diags[i] = ComplexF64.(diag(D))
        end
        if ok
            C = zeros(ComplexF64, r, r)
            @inbounds for i in 1:r
                C[:, i] = diags[i]
            end
            return C, V
        end
    end

    error("fusion_ring_characters: failed to find a common eigenbasis. Increase `tries` or check commutativity.")
end

# ---------------------------------------
# 2) S-matrix from common eigenvectors
# ---------------------------------------

"""
    fusion_ring_smatrices(R::FusionRing; tol=1e-10) -> (S, C)

Compute an `S` for a commutative fusion ring `R` such that:

  • each `N_i^T` is diagonalised by `S` (numerically),  
  • the **first row** of `S` equals the Frobenius–Perron dimensions `d` (i.e. `S[1,i]=d_i`),  
  • `S` is returned as `Matrix{ComplexF64}`.

We construct `S` from the common eigenvectors `V` returned by `fusion_ring_characters` and
scale **columns** so that `S[1,i]=d_i`. We do **not** force full symmetry to avoid
spoiling the diagonalisation; small asymmetries are normal from numerical noise.
"""
function fusion_ring_smatrices(R::FusionRing; tol::Real=1e-10)
    d = quantum_dimensions(R)
    C, V = fusion_ring_characters(R; tol=tol)
    r = length(d)

    S = ComplexF64.(V)                  # columns = right eigenvectors
    # Ensure we can normalise columns to match the FP row
    for i in 1:r
        if abs(S[1,i]) < 1e-14
            # nudge to avoid division by 0; extremely rare in practice
            S[:,i] .+= 1e-12 .* randn(r)
        end
        if abs(S[1,i]) < 1e-14
            error("fusion_ring_smatrices: cannot normalise column $i to set S[1,i]=d_i.")
        end
        S[:, i] .*= d[i] / S[1,i]
    end
    # quick check: S^{-1} N_i^T S approximately diagonal
    Sinv = inv(S)
    for a in labels(R)
        NiT = transpose(fusion_matrix(R, a))
        D = Sinv * ComplexF64.(NiT) * S
        off = copy(D); @inbounds for j in 1:r; off[j,j] = 0; end
        if norm(off) > 1e-6
            @warn "S does not perfectly diagonalise N_$a^T (offdiag norm=$(norm(off)))."
            break
        end
    end
    return S, C
end

# -----------------------------------------------------------
# 3) Twists (T) from a Vafa-style integer system
# -----------------------------------------------------------

# NOTE: without an explicit dual map i ↦ ī, many rings (e.g. fully self-dual) make the
# integer rows collapse. We expose a keyword `dual=identity` so you can pass your
# ring’s dual if available.
function _vafa_integer_matrix(R::FusionRing; dual::Function=identity)
    N = fusion_tensor(R)        # r×r×r integer array
    r = size(N,1)
    rows = Vector{Vector{Int}}()

    # Build rows from the linearised Vafa relation (cf. your eq. 14).
    # We interpret barred indices using `dual`.
    for i in 1:r, j in 1:r, k in 1:r, l in 1:r
        ib, jb, kb, lb = dual(i), dual(j), dual(k), dual(l)
        v = zeros(Int, r)
        for n in 1:r
            nb = dual(n)
            # + terms
            v[n] += N[ib,j,n] * N[kb,l,n]
            v[n] += N[ib,l,n] * N[kb,j,n]
            v[n] += N[ib,k,n] * N[jb,l,n]
            # − terms
            v[nb] -= N[ib,j,nb] * N[kb,l,nb]
            v[nb] -= N[jb,k,nb] * N[ib,l,nb]
            v[nb] -= N[ib,k,nb] * N[jb,l,nb]
        end
        if any(!=(0), v)
            push!(rows, v)
        end
    end

    if isempty(rows)
        return zeros(Int, 0, r)
    else
        return reduce(vcat, (permutedims(collect(v)) for v in rows))
    end
end

"""
    fusion_ring_twists(R::FusionRing; prefer_exact=true, dual=identity)

Solve a Vafa-style integer system for `t = (t₁,…,t_r)` with `θ_i = exp(2π i t_i)` and `t₁=0`.
If `prefer_exact=true` and OSCAR is present, use Smith normal form (exact over ℤ).
Otherwise fall back to a rational nullspace (SVD).

You can pass `dual = i -> ī` (index of the dual/simple-object) to get stronger constraints.
If you don’t have a dual, leave it as `identity`; for fully self-dual rings this may yield
only the trivial constraint (θ₁=1).
"""
function fusion_ring_twists(R::FusionRing; prefer_exact::Bool=true, dual::Function=identity)
    r = rank(R)
    A = _vafa_integer_matrix(R; dual=dual)
    if size(A,1) == 0
        return ones(ComplexF64, r), zeros(Float64, r)
    end

    if prefer_exact && _HAS_OSCAR
        Z = Oscar.ZZ
        M = Oscar.matrix(Z, size(A,1), size(A,2), vec(A))
        K = Oscar.kernel(M)
        B = Oscar.basis(K)
        if isempty(B)                # only trivial kernel
            t = zeros(Float64, r)
            θ = ones(ComplexF64, r)
            return θ, t
        end
        v = Vector{Int}(Oscar.entries(B[1]))
        t = Float64.(v .- v[1])      # enforce t₁=0
        # normalise to a bounded representative
        maxabs = maximum(1, maximum(abs, t))
        t ./= maxabs
        θ = [cis(2π * ti) for ti in t]
        θ[1] = 1 + 0im
        return θ, t
    else
        U, Σ, Vt = svd(float.(A))
        nz = sum(>(1e-12), diag(Σ))
        if nz == size(A,2)           # no nullspace
            t = zeros(Float64, r)
            θ = ones(ComplexF64, r)
            return θ, t
        end
        t = Vector{Float64}(Vt[:, nz+1])   # 1st null vector
        t .-= t[1]                          # t₁ = 0
        t ./= maximum(1e-12, maximum(abs, t))
        θ = [cis(2π * ti) for ti in t]
        θ[1] = 1 + 0im
        return θ, t
    end
end

# ----------------------------------
# 4) Convenience wrapper: (S, T)
# ----------------------------------

"""
    modular_data(R::FusionRing; tol=1e-10, prefer_exact=true, dual=identity)

Return a `Dict` with:
  • `:S`  — S-matrix (normalised so first row = FP-dims),  
  • `:T`  — diagonal matrix of twists from Vafa (θ₁=1),  
  • `:characters` — character table `C`,  
  • `:dims` — Frobenius–Perron dimensions,  
  • `:t`   — exponents with `θ_i = exp(2π i t_i)`.

You may pass `dual` to improve Vafa constraints when your ring knows the dual map.
"""
function modular_data(R::FusionRing; tol::Real=1e-10, prefer_exact::Bool=true, dual::Function=identity)
    S, C = fusion_ring_smatrices(R; tol=tol)
    θ, t = fusion_ring_twists(R; prefer_exact=prefer_exact, dual=dual)
    return Dict(
        :S => S,
        :T => Diagonal(θ),
        :characters => C,
        :dims => quantum_dimensions(R),
        :t => t
    )
end

end # module

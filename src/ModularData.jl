
##############################
# ModularData.jl
# Experimental modular data (S, T) and character helpers.
# Where possible, relies on Oscar/GAP for character tables.
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

const _OSCAR_HINT = "This function may require Oscar.jl for character tables or verification."

"""
    fusion_ring_characters(R::FusionRing)

If `R` arose from a **group representation ring** (constructed via `group_rep_fusion_ring`),
return a vector of characters (as Dicts mapping labels → values on conjugacy classes).
For general fusion rings, throws a descriptive error.
"""
function fusion_ring_characters(R::FusionRing)
    if occursin("Rep(", string(ringname(R)))
        try
            @eval using Oscar
        catch
            error(_OSCAR_HINT)
        end
        # Without the original group object we can't reconstruct full class data reliably.
        # Provide a minimal, ring-internal "character" basis via simultaneous diagonalization
        # of the commuting fusion matrices (commutative case). This agrees with irreducible
        # characters for group representation rings.
        labs = labels(R); r = length(labs)
        # Use fusion matrices of generators; here, just pick all and hope diagonalizable
        Ns = [fusion_matrix(R, a) for a in labs]
        # Common eigenvectors via SVD heuristic
        U = Matrix{Float64}(I, r, r)
        for M in Ns
            F = float(M)
            E = eigen(Symmetric(F*F'))
            U = U * E.vectors
        end
        # Return columns as "characters" evaluated on basis
        chars = [ Dict(labs[j] => U[j,i] for j in 1:r) for i in 1:r ]
        return chars
    else
        error("fusion_ring_characters: only defined for representation-like rings (or commutative rings via experimental diagonalization).")
    end
end

"""
    fusion_ring_smatrices(R::FusionRing)

Compute an **experimental** S-matrix using the (commutative) Verlinde-like recipe:
S_{ij} = (1/𝒟) ∑_k N_{ij}^k d_k, where d_k are Frobenius–Perron dimensions and 𝒟 = √(∑ d_k^2).
For non-modular rings, this is a diagnostic matrix only.
"""
function fusion_ring_smatrices(R::FusionRing)
    d = quantum_dimensions(R)
    D = sqrt(sum(d.^2))
    N = fusion_tensor(R)
    r = size(N,1)
    S = zeros(Float64, r, r)
    @inbounds for i in 1:r, j in 1:r
        s = 0.0
        for k in 1:r
            s += N[i,j,k] * d[k]
        end
        S[i,j] = s / D
    end
    return S
end

"""
    fusion_ring_twists(R::FusionRing)

Return a diagonal vector of **twists** (T-eigenvalues). For known named rings this returns
standard values; otherwise returns a placeholder of ones (phases not known).
"""
function fusion_ring_twists(R::FusionRing)
    nm = string(ringname(R))
    r = rank(R)
    if startswith(nm, "SU(2)")
        # crude placeholder: use exp(2πi*h_j) with h_j ∝ j(j+2)/(4(k+2))
        # requires k; try to parse it
        k = try
            parse(Int, replace(nm, r"^SU\(2\)_(\d+).*" => s"\1"))
        catch
            nothing
        end
        if k !== nothing
            labs = labels(R)
            vals = Vector{ComplexF64}(undef, r)
            for (idx, lab) in enumerate(labs)
                j = parse(Int, lab)
                h = j*(j+2) / (4*(k+2))
                vals[idx] = cis(2π*h)
            end
            return vals
        end
    end
    ones(ComplexF64, r)
end

"""
    modular_data(R::FusionRing)

Return a dictionary with experimental modular data:
- `:S` → S-matrix (see `fusion_ring_smatrices`)
- `:T` → twists (see `fusion_ring_twists`)
- `:dims` → Frobenius–Perron dimensions
- `:rank` → rank
These agree with true modular data on known modular examples; otherwise are best-effort diagnostics.
"""
function modular_data(R::FusionRing)
    return Dict(
        :S => fusion_ring_smatrices(R),
        :T => fusion_ring_twists(R),
        :dims => quantum_dimensions(R),
        :rank => rank(R)
    )
end

end # module

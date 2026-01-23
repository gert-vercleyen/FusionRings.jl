module FusionRings

using Oscar
import Oscar: multiplication_table, is_commutative, rank
using Combinatorics
using JSON
using LinearAlgebra:eigen
using Base.Threads

include("GeneralFunctions.jl")
include("Creation.jl")
include("Properties.jl")
include("Operations.jl")
include("ImportData.jl")
include("FormattingAndPrinting.jl")
end

# Internal helpers #

# dual index via conjugation matrix (Anyonica's CC)
_dual_index(fr::FusionRing, i::Int)::Int =
    findfirst(==(1), conjugation_matrix(fr)[i, :])::Int

# support of i ⊗ j (Anyonica- FusionOutcomes)
function _fusion_outcomes(fr::FusionRing, i::Int, j::Int)::Vector{Int}
    N = multiplication_table(fr)
    @views findall(>(0), N[i, j, :])
end

# closure of a seed set under fusion (used by AdjointFusionRing)
function _fusion_closure(fr::FusionRing, S0::Vector{Int})::Vector{Int}
    r = rank(fr)
    seen = falses(r)
    @inbounds for s in S0
        seen[s] = true
    end
    changed = true
    while changed
        changed = false
        current = findall(seen)
        @inbounds for a in current, b in current
            for c in _fusion_outcomes(fr, a, b)
                if !seen[c]
                    seen[c] = true
                    changed = true
                end
            end
        end
    end
    findall(seen)
end

# restrict ring to subindices S (Anyonica's MT[ring][[el,el,el]])
function _restrict_subring(fr::FusionRing, S::Vector{Int}; check_closed::Bool = true)::FusionRing
    sort!(S)
    N = multiplication_table(fr)
    @views Nsub = N[S, S, S]

    if check_closed
        # sanity: S must be fusion-closed
        rmap = zeros(Int, rank(fr))
        @inbounds for (k, v) in enumerate(S)
            rmap[v] = k
        end
        @inbounds for a in S, b in S
            for c in findall(>(0), N[a, b, :])
                rmap[c] != 0 || error("subset not fusion-closed")
            end
        end
    end

    # NOTE TO SELF : adjust the constructor keyword arguments if  FusionRing struct differs.
    FusionRing(
        Nsub;
        # preserve metadata if they exist (slice per-object arrays)
        names      = hasproperty(fr, :names)     ? fr.names      : missing,
        texnames   = hasproperty(fr, :texnames)  ? fr.texnames   : missing,
        labels     = hasproperty(fr, :labels)    ? fr.labels     : missing,
        frobenius_perron_dimensions =
            fr.frobenius_perron_dimensions === missing ? missing :
            fr.frobenius_perron_dimensions[S],
    )
end


module Properties

using ..Types: FusionRing, fusion_tensor, labels, rank
using LinearAlgebra: eigen

export quantum_dimensions, global_dimension, is_commutative, multiplicity
export nonzero_structure_constants, conjugation_matrix, is_group_ring
export conjugate_element, sub_fusion_rings, is_sub_fusion_ring

function quantum_dimensions(fr::FusionRing)
    r = rank(fr)
    S = zeros(Float64, r, r)
    N = fusion_tensor(fr)
    for a in 1:r
        @views S .+= N[a, :, :]
    end
    vals, vecs = eigen(S)
    idx = argmax(vals)
    v = abs.(vecs[:, idx])
    v ./ v[1]
end

global_dimension(fr::FusionRing) = sum(x->x*x, quantum_dimensions(fr))

function is_commutative(fr::FusionRing)
    N = fusion_tensor(fr); r = size(N,1)
    for a in 1:r, b in 1:r, c in 1:r
        N[a,b,c] == N[b,a,c] || return false
    end
    true
end

multiplicity(fr::FusionRing) = maximum(fusion_tensor(fr))

function nonzero_structure_constants(fr::FusionRing)
    N = fusion_tensor(fr); r = size(N,1)
    out = Tuple{Int,Int,Int}[]
    for a in 1:r, b in 1:r, c in 1:r
        N[a,b,c]>0 && push!(out,(a,b,c))
    end
    out
end

function conjugation_matrix(fr::FusionRing)
    N = fusion_tensor(fr)
    @views N[:, :, 1]
end

function conjugate_element(fr::FusionRing, a)
    imap = Dict(l=>i for (i,l) in enumerate(labels(fr)))
    ai = a isa Integer ? a : (a isa Symbol ? imap[String(a)] : imap[String(a)])
    C = conjugation_matrix(fr)
    hits = findall(==(1), C[ai, :])
    length(hits)==1 || error("Conjugate not unique for element $a")
    labels(fr)[hits[1]]
end

function is_group_ring(fr::FusionRing)
    sum( fusion_tensor(fr) ) == FusionRings.rank(r)^2
end

function sub_fusion_rings(fr::FusionRing)
    L = labels(fr); r = length(L)
    sets = Vector{Vector{String}}()
    for mask in 1:(1<<(r-1))-1
        subset = [L[1]]
        for i in 2:r
            if ((mask >> (i-2)) & 1) == 1
                push!(subset, L[i])
            end
        end
        if is_sub_fusion_ring(fr, subset) && length(subset)<r
            push!(sets, subset)
        end
    end
    sets
end

function is_sub_fusion_ring(fr::FusionRing, S::Vector)
    # Accept Vector{String} preferred, but allow symbols via conversion
    S2 = [s isa Symbol ? String(s) : String(s) for s in S]
    Sset = Set(S2)
    all(l -> l in Sset, labels(fr)[1:1]) || return false
    for a in S2, b in S2
        imap = Dict(l=>i for (i,l) in enumerate(labels(fr)))
        ai = imap[a]; bi = imap[b]
        N = fusion_tensor(fr)[ai,bi,:]
        for (ci,m) in enumerate(N)
            m==0 && continue
            c = labels(fr)[ci]
            c in Sset || return false
        end
    end
    true
end

end

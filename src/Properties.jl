
module Properties

using ..Types: FusionRing, fusion_tensor, labels, rank
using ..GeneralFunctions: indexmap
using LinearAlgebra: eigen

export fpdims, fpdim, is_commutative, multiplicity
export nonzero_structure_constants, conjugation_matrix, is_group_ring
export conjugate_element, conjugate_label, sub_fusion_rings, is_sub_fusion_ring

function fpdims(fr::FusionRing)
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

fpdim(fr::FusionRing) = sum(x->x*x, fpdims(fr))

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

"""
    conjugate_element(fr, a) -> Int

Return the integer index of the dual (conjugate) simple object of `a`.
Accepts an integer index, a `String`, or a `Symbol`.

Rationale: internal computations (e.g. composing with other index-based
operations) are simpler when the result is an index rather than a label.
Use `conjugate_label` if you need the string form.
"""
function conjugate_element(fr::FusionRing, a)
    imap = indexmap(fr)
    ai = a isa Integer ? a : imap[String(a)]
    C = conjugation_matrix(fr)
    findfirst(==(1), C[ai, :])::Int
end

"""
    conjugate_label(fr, a) -> String

Return the label (String) of the dual simple object. Thin wrapper over
`conjugate_element`.
"""
function conjugate_label(fr::FusionRing, a)
    labels(fr)[conjugate_element(fr, a)]
end

function is_group_ring(fr::FusionRing)
    sum( fusion_tensor(fr) ) == FusionRings.rank(r)^2
=======
    all(isapprox.(fpdims(fr), 1.0; atol=1e-8)) || return false
    N = fusion_tensor(fr); r = size(N,1)
    for a in 1:r, b in 1:r
        s = sum(N[a,b,:])
        s == 1 || return false
    end
    true
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
    imap = indexmap(fr)
    for a in S2, b in S2
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

"""
    is_sub_fusion_ring(big::FusionRing, small::FusionRing) -> Bool

Return true if `small` is (isomorphic to) a fusion subring of `big` on the
same labeling of simples restricted to a subset.

Criteria:
* All labels of `small` must appear among `big`'s labels (string match).
* The fusion tensors must agree for all triples of indices inside that subset.

Note: This is a structural containment check under the identity embedding
determined by matching labels; it does not attempt to permute labels of the
subring to achieve a match. Use an equivalence/permutation utility first if
you need to test up to relabeling.
"""
function is_sub_fusion_ring(big::FusionRing, small::FusionRing)
    Lbig = labels(big); Lsmall = labels(small)
    # Quick label inclusion
    label_pos = Dict{String,Int}()
    for (i,l) in enumerate(Lbig); label_pos[l] = i; end
    idxs = Int[]
    for l in Lsmall
        i = get(label_pos, l, 0)
        i == 0 && return false
        push!(idxs, i)
    end
    Nbig = fusion_tensor(big)
    Nsmall = fusion_tensor(small)
    r2 = length(Lsmall)
    @inbounds for a in 1:r2, b in 1:r2, c in 1:r2
        if Nsmall[a,b,c] != Nbig[idxs[a], idxs[b], idxs[c]]
            return false
        end
    end
    true
end

end

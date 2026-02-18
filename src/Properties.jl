
import Combinatorics: combinations


function _internal_multiplication(fr::FusionRing, S::Vector{Int})::Bool
    Sset = Set(S)
    @inbounds for i in S, j in S
        for c in fusion_outcomes(fr, i, j)
            c in Sset || return false
        end
    end
    true
end

"""
    _internal_closed_subsets(fr, k) -> Vector{Vector{Int}}

Return all fusion-closed subsets of size `k` containing the unit `1`,
generated as `S = [1; T]` where `T` ranges over (k-1)-subsets of `2:r`.

This avoids ever generating candidates that omit the vacuum.
"""
function _internal_closed_subsets(fr::FusionRing, k::Int)::Vector{Vector{Int}}
    r = rank(fr)
    (k <= 0 || k > r) && return Vector{Vector{Int}}()

    if k == 1
        return _internal_multiplication(fr, [1]) ? [[1]] : Vector{Vector{Int}}()
    end

    candidates = Vector{Vector{Int}}()
    for T in combinations(collect(2:r), k - 1)
        S = vcat(1, collect(T))
        _internal_multiplication(fr, S) && push!(candidates, S)
    end
    candidates
end


"""
    _diag_channel_groups(N) -> Vector{Vector{Int}}

Partition indices by the invariant:
k(i) = |{ c : N[i,i,c] > 0 }|.

We return groups in deterministic order:
- increasing k
- increasing indices within each group
"""
function _diag_channel_groups(N::Array{Int,3})::Vector{Vector{Int}}
    r = size(N, 1)

    k = Vector{Int}(undef, r)
    @inbounds for i in 1:r
        cnt = 0
        for c in 1:r
            (N[i,i,c] > 0) && (cnt += 1)
        end
        k[i] = cnt
    end

    groups = Dict{Int,Vector{Int}}()
    @inbounds for i in 1:r
        push!(get!(groups, k[i], Int[]), i)
    end

    out = Vector{Vector{Int}}()
    for kk in sort!(collect(keys(groups)))
        g = groups[kk]
        sort!(g)
        push!(out, g)
    end
    out
end

# Apply permutation P on all three indices: A'[i,j,k] = A[P[i],P[j],P[k]]
function _permute_multtab(A::Array{Int,3}, P::Vector{Int})::Array{Int,3}
    r = size(A, 1)
    B = similar(A)
    @inbounds for i in 1:r, j in 1:r, k in 1:r
        B[i,j,k] = A[P[i], P[j], P[k]]
    end
    B
end

"""
    _permutation_vector_equiv(A, B) -> Vector{Int} or nothing

Find permutation `perm` with `_permute_multtab(A, perm) == B`, using
diagonal-channel groups for pruning. Returns `nothing` if not found.
"""
function _permutation_vector_equiv(A::Array{Int,3}, B::Array{Int,3})
    r = size(A, 1)
    size(B, 1) == r || return nothing

    grpA = _diag_channel_groups(A)
    grpB = _diag_channel_groups(B)
    sort(map(length, grpA)) == sort(map(length, grpB)) || return nothing

    used = falses(length(grpB))
    cur  = Vector{Int}(undef, r)
    cur[1] = 1  # unit fixed

    function backtrack(gidx::Int)::Bool
        if gidx > length(grpA)
            return _permute_multtab(A, cur) == B
        end
        GA = grpA[gidx]
        for j in eachindex(grpB)
            (used[j] || length(grpB[j]) != length(GA)) && continue
            used[j] = true
            for σ in Base.Iterators.permutations(grpB[j])
                if 1 in GA
                    σ[findfirst(==(1), GA)] == 1 || continue
                end
                for (u, v) in zip(GA, σ)
                    cur[u] = v
                end
                backtrack(gidx + 1) && return true
            end
            used[j] = false
        end
        false
    end

    backtrack(1) ? cur : nothing
end


export which_injection
"""
    which_injection(subring, ring) -> Dict{Int,Int} or nothing

Find an injection of `subring` into `ring` by:
1) enumerating fusion-closed subsets `S ⊂ ring` of size rank(subring) containing 1,
2) comparing multiplication tables up to permutation.
"""
function which_injection(subring::FusionRing, ring::FusionRing)
    rs = rank(subring)
    rr = rank(ring)
    rs > rr && return nothing

    Nbig = multiplication_table(ring)
    Nsub = multiplication_table(subring)

    for S in _internal_closed_subsets(ring, rs)
        Nres = @views Nbig[S, S, S]
        perm = _permutation_vector_equiv(Nsub, Nres)
        perm === nothing && continue

        inj = Dict{Int,Int}()
        @inbounds for i in 1:rs
            inj[i] = S[perm[i]]
        end
        return inj
    end
    nothing
end

export fusion_ring_automorphisms
"""
    fusion_ring_automorphisms(fr) -> Vector{Vector{Int}}

Return all permutations `p` with `_permute_multtab(N,p) == N`.
Uses diagonal-channel pruning (same idea as Anyonica).
"""
function fusion_ring_automorphisms(fr::FusionRing)
    N = multiplication_table(fr)
    r = size(N, 1)

    groups = _diag_channel_groups(N)
    perms  = Vector{Vector{Int}}()
    cur    = Vector{Int}(undef, r)
    cur[1] = 1

    function backtrack(gidx::Int)
        if gidx > length(groups)
            _permute_multtab(N, cur) == N && push!(perms, copy(cur))
            return
        end
        G = groups[gidx]
        for σ in Base.Iterators.permutations(G)
            if 1 in G
                σ[findfirst(==(1), G)] == 1 || continue
            end
            for (u, v) in zip(G, σ)
                cur[u] = v
            end
            backtrack(gidx + 1)
        end
    end

    backtrack(1)
    unique!(perms)
    sort!(perms, by = p -> (sum(p), p))
    perms
end


export commutator
"""
    commutator(fr, sub) -> (els, subring)

commutator (centralizer) of a subring `sub = (S, RS)` inside `fr`:
return all simples `i` with `FusionOutcomes(i ⊗ i*) ⊆ S`.

 the raw set `els0` need not be fusion-closed, so we return the
*smallest fusion subring generated by els0*, i.e. `els = _fusion_closure(fr, els0)`,
and then restrict.
"""
function commutator(fr::FusionRing, sub::Tuple{Vector{Int},FusionRing})
    is_commutative(fr) || error("commutator: ring must be commutative ")

    subEls, _ = sub
    Sset = Set(subEls)

    in_sub(i::Int)::Bool = begin
        di = conjugate_element(fr, i)
        outs = fusion_outcomes(fr, i, di)   
        all(in(Sset), outs)
    end

    els0 = [i for i in 1:rank(fr) if in_sub(i)]
    isempty(els0) && (els0 = [1])

    els = _fusion_closure(fr, els0)  
    return els, _restrict_subring(fr, els; check_closed=true)
end

function commutator(fr::FusionRing)
    return derived_subring_commutator(fr)
end

export universal_grading
"""
    universal_grading(fr) -> (grading, groupRing)


- Compute `irreps = adjoint_irreps(fr)` (a partition of simples).
- Let the universal grading group have `n = length(irreps)` elements.
- `grading` maps each simple `x` to the index `a` of the adjoint-irrep block containing it.
- Define multiplication on group elements a,b by:
    N[a,b,c] = 1  iff  FusionOutcomes(i ⊗ j) ⊆ irreps[c]
for all i ∈ irreps[a], j ∈ irreps[b].
"""
function universal_grading(fr::FusionRing)
    irreps = adjoint_irreps(fr)
    n = length(irreps)

    # grading: Vector{Pair{Int,Int}} mapping (simple -> grade)
    grading = Pair{Int,Int}[]
    @inbounds for a in 1:n
        for x in irreps[a]
            push!(grading, x => a)
        end
    end
    sort!(grading, by = p -> first(p))

    # helper: check cond(l1,l2,l3)
    function _cond(l1::Vector{Int}, l2::Vector{Int}, l3::Vector{Int})::Bool
        S = Set(l3)
        @inbounds for i in l1, j in l2
            outs = fusion_outcomes(fr, i, j)
            for c in outs
                c in S || return false
            end
        end
        true
    end

    mt = zeros(Int, n, n, n)
    @inbounds for a in 1:n, b in 1:n, c in 1:n
        mt[a,b,c] = _cond(irreps[a], irreps[b], irreps[c]) ? 1 : 0
    end

    # build the grading group ring
    groupRing = fusion_ring(mt; labels = string.(1:n))
    return grading, replace_by_known(groupRing)
end

export UG
UG(fr::FusionRing) = universal_grading(fr)

function all_gradings(fr::FusionRing)
    error("all_gradings not implemented yet")
end

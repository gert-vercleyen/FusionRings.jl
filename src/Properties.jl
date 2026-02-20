import Combinatorics: combinations

#
# Changes made 
# -  Replaced custom k_subsets with combinations
#    and forced vacuum inclusion in subring search.
#
# -  Removed redundant `_fusion_outcomes`.
#    
#
# - Clarified diagonal-channel pruning logic
#
# - Fixed commutator to return smallest generated
#    fusion-closed subring instead of raw candidate set.
#
#    - Moved structural logic fully into Properties.jl


#function change_fusion_ring_property(r::FusionRing, dict)
#end

export multiplication_table
function multiplication_table(r::FusionRing)::Array{Int,3}
  return r.multiplication_table
end

# pmt = print_multiplication_table

function element_to_string(mult, elem)::String
  if mult == 0
    return ""
  elseif mult == 1
    return elem
  else
    return string(mult) * " " * elem
  end
end

export rank
function rank(r::FusionRing)::Int
  size(multiplication_table(r))[1]
end

export names
function names(r::FusionRing)::Array{String,1}
  return r.names
end

export tex_names
function tex_names(r::FusionRing)::Array{String,1}
  return r.texnames
end

export labels
function labels(r::FusionRing)::Array{String,1}
  return r.labels
end

export conjugation_matrix
function conjugation_matrix(fr::FusionRing)
    N = fusion_tensor(fr)
    @views N[:, :, 1]
end

export multiplicity
function multiplicity(r::FusionRing)::Int
  maximum(multiplication_table(r))
end

export nonzero_structure_constants
function nonzero_structure_constants(r::FusionRing)::Vector{Tuple{Int64, Int64, Int64}}
  mt = multiplication_table(r)
  Tuple.(findall(x -> x > 0, mt))
end

export nzsc
nzsc = nonzero_structure_constants

function num_nonzero_structure_constants(r::FusionRing)::Int64
  length(nzsc(r))
end
nnzsc = num_nonzero_structure_constants

export frobenius_perron_dimensions
function frobenius_perron_dimensions(r::FusionRing)::Vector{QQBarFieldElem}
  stored_dims = r.frobenius_perron_dimensions
  if stored_dims === missing
    mt = multiplication_table(r)
    multmats = [matrix(ZZ, mt[i, :, :]) for i in 1:rank(r)]
    return [first(eigenvalues(QQBar, A)) for A in multmats]
  else
    return stored_dims
  end
end

export fpdims
fpdims = frobenius_perron_dimensions

export frobenius_perron_dimension
function frobenius_perron_dimension(r::FusionRing)::QQBarFieldElem
  return sum(fpdims(r) .^ 2)
end

export fpdim
fpdim = frobenius_perron_dimension

export num_self_dual_non_self_dual
function num_self_dual_non_self_dual(r::FusionRing)::Array{Int,1}
  sd  = count(x -> x == 1, diag(conjugation_matrix(r)))
  nsd = rank(r) - sd
  return [sd nsd]
end

export nsdnsd
nsdnsd = num_self_dual_non_self_dual

export num_self_dual
function num_self_dual(r::FusionRing)::Int
  first(nsdnsd(r))
end

export nsd
nsd = num_self_dual

export num_non_self_dual
function num_non_self_dual(r::FusionRing)::Int
  last(nsdnsd(r))
end

export nnsd
nnsd = num_non_self_dual

export is_group_ring
function is_group_ring(r::FusionRing)::Bool
  return sum(multiplication_table(r)) == rank(r)^2
end


export conjugate_element

# 1-arg: returns a function i -> dual(i)
function conjugate_element(r::FusionRing)
   return i -> (conjugation_matrix(r) * collect(1:rank(r)))[i]
end

"""
    conjugate_element(fr, a) -> Int

Return the integer index of the dual (conjugate) simple object of `a`.
Accepts an integer index, a `String`, or a `Symbol`.
"""
function conjugate_element(fr::FusionRing, a)
    imap = indexmap(fr)
    ai = a isa Integer ? a : imap[String(a)]
    C = conjugation_matrix(fr)
    findfirst(==(1), C[ai, :])::Int
end

export anyonwiki_code
function anyonwiki_code(r::FusionRing)::Array{Int,1}
  return r.anyonwiki_code
end

export barcode
function barcode(r::FusionRing)::Int
  return r.barcode
end

function mult_tab_code(mat::Array{Int,2}, mult::Int)::Int
end

export sub_fusion_rings
function sub_fusion_rings(r::FusionRing)
  return r.sub_fusion_rings
end

function sub_ring_tables(mat::Array{Int,2})
end

function injection_form(subring::FusionRing, ring::FusionRing)
end

function is_sub_fusion_ring(subring::FusionRing, ring::FusionRing)::Bool
end

function is_equivalent_fusion_ring(ring1::FusionRing, ring2::FusionRing)::Bool
end

function permutation_vector(mt1::Array{Int,3}, mt2::Array{Int,3})::Array{Int,1}
end

function automorphisms(r::FusionRing)::Array{Int,2}
end

export decompositions
function decompositions(r::FusionRing, product="TensorProduct")::Array{FusionRing,1}
  if product == "TensorProduct"
    return r.tensor_product_decompositions
  else
    return error("Only tensor product decompositions are defined at the moment.")
  end
end


export adjoint_fusion_ring
function adjoint_fusion_ring(ring::FusionRing)::Tuple{Vector{Int},FusionRing}
    d(i) = conjugate_element(ring, i)

    el_seen = falses(rank(ring))
    for (i, j, c) in nzsc(ring)
        if j == d(i)
            el_seen[c] = true
        end
    end
    el = findall(el_seen)

    generatedEl = _fusion_closure(ring, el)

    if length(generatedEl) == rank(ring)
        return (generatedEl, ring)
    else
        return (generatedEl, _restrict_subring(ring, generatedEl; check_closed = true))
    end
end


function upper_central_series(r::FusionRing)::Array{FusionRing,1}
end

function is_nilpotent(r::FusionRing)::Bool
end

function adjoint_irreps(r::FusionRing)::Array{Array{Int,1},1}
end


# TODO 
#riginal issue:
# - There was a third function named `_fusion_outcomes`
#   that duplicated functionality already provided elsewhere.



function _internal_multiplication(fr::FusionRing, S::Vector{Int})::Bool
    Sset = Set(S)
    @inbounds for i in S, j in S
        for c in fusion_outcomes(fr, i, j)
            c in Sset || return false
        end
    end
    true
end


# TODO 
#
# Original issue:
# -  previously generated all k-subsets of {1,...,r}.
# - This included subsets not containing the vacuum element 1.
# - Such subsets can never define fusion subrings.
#
# 
#
#   changed:
# - Replaced custom k_subsets logic with Combinatorics.combinations.
# - Only generate subsets S = [1; T], where T ⊂ {2,...,r}.
# -  avoids unnecessary candidates 


"""
    _internal_closed_subsets(fr, k) -> Vector{Vector{Int}}

Return all fusion-closed subsets of size `k` containing the unit `1`,
generated as `S = [1; T]` where `T` ranges over (k-1)-subsets of `2:r`.
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

# TODO

# -  compute the invariant:
#       k(i) = |{ c : N[i,i,c] > 0 }|
# - We partition indices by this invariant.
# - Groups are sorted deterministically (increasing k, then index).

"""
    _diag_channel_groups(N) -> Vector{Vector{Int}}

Partition indices by invariant k(i)=|{c : N[i,i,c]>0}|.

Return groups in deterministic order:
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

Find `perm` such that `_permute_multtab(A, perm) == B`, using diagonal-channel
groups for pruning. Returns `nothing` if not found.
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



function derived_subring_commutator(fr::FusionRing)::FusionRing
    r = rank(fr)
    return derived_subring_commutator(fr, collect(1:r), collect(1:r))
end

function derived_subring_commutator(fr::FusionRing, A::Vector{Int}, B::Vector{Int})::FusionRing
    r = rank(fr)
    all(1 .≤ A .≤ r) || throw(ArgumentError("derived_subring_commutator: A has out-of-bounds indices"))
    all(1 .≤ B .≤ r) || throw(ArgumentError("derived_subring_commutator: B has out-of-bounds indices"))

    seen = falses(r)
    @inbounds for a in A
        aᵗ = conjugate_element(fr, a)
        for b in B
            bᵗ = conjugate_element(fr, b)

            # (a ⊗ b) ⊗ a* ⊗ b*
            for u in fusion_outcomes(fr, a, b)
                for v in fusion_outcomes(fr, u, aᵗ)
                    for w in fusion_outcomes(fr, v, bᵗ)
                        seen[w] = true
                    end
                end
            end
        end
    end

    S0 = findall(seen)
    isempty(S0) && (S0 = [1])

    S = _fusion_closure(fr, S0)
    return _restrict_subring(fr, S; check_closed=true)
end

# TODO 
#
# Original issue:
# - Commutator collected all elements x such that xx* ∈ S.
# -  directly passed that raw set to `_restrict_subring`.
# - This can fail if the set is not fusion-closed.
#

# What  changed:
# - First compute raw candidate set `els0`.
# - Then compute:
#       els = _fusion_closure(fr, els0)
# - then restrict to subring.

export commutator
"""
    commutator(fr, sub) -> (els, subring)

Centralizer-style commutator of subring `sub = (S, RS)` inside `fr`:
return all simples `i` with `FusionOutcomes(i ⊗ i*) ⊆ S`.

Fix per your professor: the raw set `els0` need not be fusion-closed, so we return
the smallest fusion subring generated by `els0`, i.e. `els = _fusion_closure(fr, els0)`,
then restrict.
"""
function commutator(fr::FusionRing, sub::Tuple{Vector{Int},FusionRing})
    is_commutative(fr) || error("commutator: ring must be commutative")

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

# TODO 
#

#   1) Compute adjoint_irreps (partition of simples).
#   2) Define grading map from simples to block indices.
#   3) Construct multiplication table on grading blocks:
#        mt[a,b,c] = 1  iff
#            FusionOutcomes(i ⊗ j) ⊆ irreps[c]
#        for all i ∈ irreps[a], j ∈ irreps[b].
#
#  implements  universal grading group of the fusion ring.

export universal_grading
"""
    universal_grading(fr) -> (grading, groupRing)



- Compute `irreps = adjoint_irreps(fr)` (partition of simples).
- Create group object with `n = length(irreps)` elements.
- `grading` maps each simple `x` to block index `a`.
- Multiplication table on the grading group is:
    mt[a,b,c] = 1  iff  FusionOutcomes(i ⊗ j) ⊆ irreps[c]
for all i ∈ irreps[a], j ∈ irreps[b].
"""
function universal_grading(fr::FusionRing)
    irreps = adjoint_irreps(fr)
    n = length(irreps)

    grading = Pair{Int,Int}[]
    @inbounds for a in 1:n
        for x in irreps[a]
            push!(grading, x => a)
        end
    end
    sort!(grading, by = p -> first(p))

    function _cond(l1::Vector{Int}, l2::Vector{Int}, l3::Vector{Int})::Bool
        S = Set(l3)
        @inbounds for i in l1, j in l2
            for c in fusion_outcomes(fr, i, j)
                c in S || return false
            end
        end
        true
    end

    mt = zeros(Int, n, n, n)
    @inbounds for a in 1:n, b in 1:n, c in 1:n
        mt[a,b,c] = _cond(irreps[a], irreps[b], irreps[c]) ? 1 : 0
    end

    groupRing = fusion_ring(mt; labels = string.(1:n))
    return grading, replace_by_known(groupRing)
end

export UG
UG(fr::FusionRing) = universal_grading(fr)


function all_gradings(fr::FusionRing)
    error("all_gradings not implemented yet")
end


export characters
function characters(ring::FusionRing; use_numerics = true )
  if !(ring.characters === missing)
    return ring.characters
  elseif !FusionRings.is_commutative(ring)
    error("Calculation of characters for non-commutative fusion ring is not implemented yet.")
  elseif rank(ring) == 1
    return [ qqbar(1) ]
  else
    mt   = multiplication_table( ring )
    r    = rank( ring )
    mats = [ mt[i,:,:] for i in 1:r ]
    qqb  = algebraic_closure(QQ)

    sort_mat( mat ) = sortslices( mat, dims = 1, by = char_sort_crit )

    if use_numerics
        V     = ( normalize_first_col ∘ numeric_diagonalizing_matrix )( mats )
        evals = union( vcat( [ eigenvalues( matrix( qqb, m ) ) for m in mats ] ... ) )
        chars = replace_by_known(evals).(V)
    else
        chars = ( normalize_first_col ∘ diagonalizing_matrix )( mats )
    end
    matrix( qqb, sort_mat( chars ) )
  end
end

function numeric_diagonalizing_matrix( mats, tries::Int=64, tol::Real=1e-12 )
    r = mats |> length
    m, n = mats |> first |> size

    function is_diagonalizing_matrix( mat, mats )
        imat = inv(mat)
        for i in 1:r
            D = imat * mats[i] * mat
            @inbounds for j in 1:r; D[j,j] = 0.0; end
            if norm(D) > tol
                return false
            end
        end
        return true
    end

    for _ in 1:tries
        coeffs = rand(Float64, r)
        M = zeros(Float64, r, r)
        @inbounds for k in 1:r
            M .+= coeffs[k] .* mats[k]
        end

        V = Matrix(eigen(M).vectors)
        is_diagonalizing_matrix(V, mats) && return inv(V)
    end

    error("No diagonalizing matrix found in $(tries) tries.")
end

function diagonalizing_matrix( mats )
  qqbmats = [ matrix( algebraic_closure(QQ), m ) for m in mats ]

  function is_diagonalizing_matrix( mat, mats )
    invmat = inv(mat)
    for m in mats
        is_diagonal(mat * m * invmat) || return false
    end
    return true
  end

  proposed_mat = qqbmats[1]
  r = first(size(first(mats)))
  diagq = false; upi = 4; upj = 4

  while !diagq
    upi += 1; upj += 1
    rvec    = rand(unique([ i//j for i ∈ 1:upi, j ∈ 1:upj ]), r)
    sgnvec  = rand([ -1 1 ], r)

    combinedmat = sgnvec[1] * rvec[1] * qqbmats[1]
    for i ∈ 1:r
      combinedmat += sgnvec[i] * rvec[i] * qqbmats[i]
    end

    proposed_mat =
      reduce(vcat, (collect ∘ values ∘ eigenspaces)(combinedmat))

    diagq = is_diagonalizing_matrix(proposed_mat, qqbmats)
  end

  return proposed_mat
end

function char_sort_crit( v )
	RR = ArbField(64)
	CC = AcbField(64)
	conv(x) = convert(Float64,x)
	absval(vec) = conv.( RR.( abs2.( vec ) ) )
	angl(vec) = conv.( real.( log.( CC.( vec ) ) ./ CC( 2 * pi * im ) ) )
	( - Int(all(isreal.(v))), angl(v), absval(v) )
end

function normalize_first_col( mat )
    m, n = size(mat)
    [ mat[i,j] / mat[i,1] for i in 1:m, j in 1:n ]
end

function is_diagonalizing_matrix( mat, ring::FusionRing )
	mt   = FusionRings.multiplication_table( ring )
	r    = FusionRings.rank(ring)
	mats = [ matrix( qqb, mt[ i, :, : ] ) for i ∈ 1:r ]
    mat  = matrix( qqb, mat )
    all( is_diagonal( mat * m * inv(mat) ) for m in mats )
end

export numeric_characters
function numeric_characters( ring, tries::Int = 64, tol=1e-12 )
  if !(ring.numeric_characters === missing)
    return ring.numeric_characters
  elseif !FusionRings.is_commutative(ring)
    error("Calculation of characters for non-commutative fusion ring is not implemented yet.")
  elseif rank(ring) == 1
    return [ qqbar(1) ]
  else
    mt   = multiplication_table( ring )
    r    = rank( ring )
    mats = [ mt[i,:,:] for i in 1:r ]
    V    = ( normalize_first_col ∘ numeric_diagonalizing_matrix )( mats )
    sortslices( V, dims = 1, by = num_char_sort_crit )
  end
end

function num_char_sort_crit( v )
    norm(vec) = sqrt(sum(abs2.(vec)))
    are_real(vec) = (Int ∘ all)(abs(imag(x)) < 1e-14 for x in vec ./ norm(vec))
	absval(vec) = abs2.(vec)
	angl(vec) = angle.(vec)
	( - are_real(v), angl(v), absval(v) )
end


export sl_2_ZZ_reps
function sl_2_ZZ_reps(r::FusionRing)
  return r.modular_data
end

function s_matrices(r::FusionRing)
end

function normalized_s_matrices(r::FusionRing)
end

function twist_factors(r::FusionRing)
end


function numeric_fpdims(fr::FusionRing)
    r = rank(fr)
    S = zeros(Float64, r, r)
    N = multiplication_table(fr)
    for a in 1:r
        @views S .+= N[a, :, :]
    end
    vals, vecs = eigen(S)
    idx = argmax(vals)
    v = abs.(vecs[:, idx])
    v ./ v[1]
end

numeric_fpdim(fr::FusionRing) = sum(x->x*x, numeric_fpdims(fr))

function is_commutative(fr::FusionRing)
    N = multiplication_table(fr)
    r = rank(fr)
    for a in 1:r, b in 1:r, c in 1:r
        N[a,b,c] == N[b,a,c] || return false
    end
    true
end

#function change_fusion_ring_property(r::FusionRing, dict)

#end


export multiplication_table

function multiplication_table(r::FusionRing)::Array{Int,3}
  return r.multiplication_table
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
  Tuple.( findall( x -> x > 0, mt ) ) 
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
  if stored_dims===missing
    mt = multiplication_table(r)
    multmats = [ matrix( ZZ, mt[i,:,:] ) for i in 1:rank(r) ]
    return [ first(eigenvalues(QQBar, A)) for A in multmats ]
  else
    return stored_dims
  end
end

export fpdims 

fpdims = frobenius_perron_dimensions

export frobenius_perron_dimension

function frobenius_perron_dimension(r::FusionRing)::QQBarFieldElem
  return sum( fpdims(r).^2 )
end

export fpdim

fpdim = frobenius_perron_dimension

export num_self_dual_non_self_dual

function num_self_dual_non_self_dual(r::FusionRing)::Array{Int,1}
  sd  = count( x -> x == 1, diag( conjugation_matrix(r) ) )
  nsd = rank(r) - sd
  return [ sd nsd ]
end

export nsdnsd

nsdnsd = num_self_dual_non_self_dual

export num_self_dual

function num_self_dual(r::FusionRing)::Int
  first( nsdnsd(r) )
end

export nsd

nsd = num_self_dual

export num_non_self_dual 

function num_non_self_dual(r::FusionRing)::Int
  last( nsdnsd(r) )
end

export nnsd 

nnsd = num_non_self_dual

export is_group_ring

function is_group_ring(r::FusionRing)::Bool
  return sum( multiplication_table(r) ) == rank(r)^2
end

export conjugate_element

function conjugate_element(r::FusionRing)
   return i -> ( conjugation_matrix(r) * collect( 1:rank(r) ) )[i]
end

export anyonwiki_code

function anyonwiki_code(r::FusionRing)::Array{Int,1}
  return r.anyonwiki_code
end

export barcode 

function barcode(r::FusionRing)
  return r.barcode
end

function mult_tab_code(mat::Array{Int,2},mult::Int)::Int
end

export sub_fusion_rings

function sub_fusion_rings(r::FusionRing)
    dictvec = r.sub_fusion_rings
    if dictvec !== missing
        [
            Dict(
                "injection"   => dict["injection"],
                "fusion_ring" => awc(dict["anyonwiki_code"])
            )
            for dict in dictvec
        ]
    else
        error("Method sub_fusion_rings not full implemented yet")
    end
# TODO: still uses labels and doesn't return injections
#function sub_fusion_rings(fr::FusionRing)
#    L = labels(fr); r = length(L)
#    sets = Vector{Vector{String}}()
#    for mask in 1:(1<<(r-1))-1
#        subset = [L[1]]
#        for i in 2:r
#            if ((mask >> (i-2)) & 1) == 1
#                push!(subset, L[i])
#            end
#        end
#        if is_sub_fusion_ring(fr, subset) && length(subset)<r
#            push!(sets, subset)
#        end
#    end
#    sets
#end
end

function sub_ring_tables(mat::Array{Int,2})

end

# injection_form( subring, ring ) returns the vector of elements of
# ring that form subring in the correct order  
function injection_form( subring::FusionRing, ring::FusionRing )

end

function is_sub_fusion_ring(subring::FusionRing,ring::FusionRing)::Bool
  
end

function is_equivalent_fusion_ring(ring1::FusionRing,ring2::FusionRing)::Bool

end

function permutation_vector( mt1::Array{Int,3}, mt2::Array{Int,3})::Array{Int,1}

end

function automorphisms(r::FusionRing)::Array{Int,2}
  
end

export decompositions

function decompositions( fr::FusionRing, product="TensorProduct" )#::Vector{ Vector{FusionRing} }
    product == "TensorProduct" ||  error("Only tensor product decompositions are defined at the moment.")

    tpd = fr.tensor_product_decompositions
    if tpd !== missing
        [ [ awc( code ) for code in decomp ] for decomp in tpd ]
    else
        tensor_product_decompositions(fr)
    end
end

function tensor_product_decompositions( r::FusionRing )
    error("Not implemented yet.")
end

function adjoint_fusion_ring(r::FusionRing)::FusionRing
  
end

function upper_central_series(r::FusionRing)::Array{FusionRing,1}
  
end

function is_nilpotent(r::FusionRing)::Bool
  
end

function adjoint_irreps(r::FusionRing)::Array{Array{Int,1},1}
  
end

function universal_grading(r::FusionRing)
  
end

function all_gradings(r::FusionRing)

end

function commutator(r::FusionRing)

end

export characters 

function characters(ring::FusionRing; use_numerics = true )
  if !(ring.characters === missing)
    return from_qqb_id( ring.characters )
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
        # find numeric diagonalizing matrix
        V     = ( normalize_first_col ∘ numeric_diagonalizing_matrix )( mats )
        # find symbolic evals
        evals = union( vcat( [ eigenvalues( matrix( qqb, m ) ) for m in mats ] ... ) )
        # identify evals in V, sort mat, and convert to matrix of qqb's
        chars = replace_by_known(evals).(V)
    else
        chars = ( normalize_first_col ∘ diagonalizing_matrix )( mats )
    end
    matrix( qqb, sort_mat( chars ) )
  end
end

function numeric_diagonalizing_matrix( mats, tries::Int=64, tol::Real=1e-12 )
    r = length( mats )

    m, n = size( first( mats ) )

    function is_diagonalizing_matrix( mat, mats )
        imat = inv(mat)
        for i in 1:r
            D = imat * mats[i] * mat
            @inbounds for j in 1:r; D[j,j] = 0.0; end #set diagonal el = 0

            if norm(D) > tol # D should be 0 mat now
                return false
            else
                continue
            end
        end
        return true
    end

    for i in 1:tries
        coeffs = rand( Float64, r )
        M = zeros( Float64, r, r )
        @inbounds for k in 1:r
            M .+= coeffs[k] .* mats[k]
        end


        V  = Matrix(eigen(M).vectors)    # symmetric not guaranteed; generic eigen

        # need to take inv sicne Oscar's (and our) definition of diagonalizing matrix
        # is inverse of that of LinearAlgebra package
        is_diagonalizing_matrix( V, mats ) && return inv( V )

        continue
    end

    error("No diagonalizing matrix found in " * string(tries) *" tries. Try setting tries to a higher value or set a different seed for random generator.")
end

function diagonalizing_matrix( mats )
  qqbmats = [ matrix( algebraic_closure(QQ), m ) for m in mats ]

  function is_diagonalizing_matrix( mat, mats )
    invmat = inv(mat)
    for m in mats
        if !is_diagonal( mat * m * invmat )
            return false
        else
            continue
        end
    end
    return true
  end

  proposed_mat = qqbmats[1]

  r = first( size( first( mats ) ) )
    
  diagq = false; upi = 4; upj = 4

  while !diagq
    upi += 1
    upj += 1

    # Take random linear rational combination of matrices in mats
    rvec    = rand( unique( [ i//j for i ∈ 1:upi, j ∈ 1:upj ] ), r )
    sgnvec  = rand( [ -1 1 ], r )
    combinedmat = sgnvec[1] * rvec[1] * qqbmats[1]
    for i ∈ 1:r
      combinedmat += sgnvec[i] * rvec[i] * qqbmats[i]
    end

    # Find diagonalizing matrix
    proposed_mat = 
      reduce( 
        vcat,
        (collect ∘ values ∘ eigenspaces)( combinedmat )
      )

    # Check whether it works
    diagq = is_diagonalizing_matrix( proposed_mat, qqbmats )
  end

  return proposed_mat 
end

# Sort criterion for characters
function char_sort_crit( v )
	RR = ArbField(64);
	CC = AcbField(64);
	conv(x) = convert(Float64,x)
	# Abs values of elements of v
	absval(vec) = conv.( RR.( abs2.( vec ) ) )
	# Angles of elements of v
	angl(vec) = conv.( real.( log.( CC.( vec ) ) ./ CC( 2 * pi * im ) ) )
			
	( - Int( all( isreal.(v) ) ) , angl(v), absval(v) )
end

function normalize_first_col( mat )
    m, n = size( mat )
    [ mat[i,j] / mat[i,1] for i in 1:m, j in 1:n ]
end

function is_diagonalizing_matrix( mat, ring::FusionRing )
	mt   = FusionRings.multiplication_table( ring )
	r    = FusionRings.rank(ring)
	mats = [ matrix( qqb, mt[ i, :, : ] ) for i ∈ 1:r ]
  mat  = matrix( qqb, mat )

  all( is_diagonal( mat * m * inv(mat) ) for m in mats )
end

#export numeric_characters
# There are some issues with the numeric characters function. 
# Mainly the fact that the rows aren't sorted according to some 
# criterion
export numeric_characters
"""
    numeric_characters(R::FusionRing; tries::Int=64, tol=1e-12) -> (C, V)

Return the **character table** `C::Matrix{ComplexF64}` of a **commutative** fusion ring `R`.

Algorithm:
1. Form a random real combination `M = ∑_k c_k N_k`.
2. Eigen-decompose `M = V Λ V⁻¹`.
3. Verify that every `V⁻¹ N_i V` is (numerically) diagonal. If not, retry.
4. Normalize and sort V

Throws if no common eigenbasis is found after `tries` attempts.
"""
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

# Numeric sort criterion for characters
function num_char_sort_crit( v )
    # for normalizing demands
    norm(vec) = sqrt( sum( abs2.(vec) ) )
    # Check whether all elements are real
    are_real(vec) = ( Int ∘ all )( abs( imag( x ) ) < 1e-14 for x in vec ./ norm(vec) )
	# Abs values of elements of v
	absval(vec) = abs2.( vec )
	# Angles of elements of v
	angl(vec) = angle.( vec )

	( - are_real(v) , angl(v), absval(v) )
end

# finds diagonalizing matrix using floating point arithmetic

export projective_SL_2_ZZ_reps

function projective_SL_2_ZZ_reps( fr::FusionRing )
    md = fr.projective_SL2Z_reps
    if md !== missing 
        return md
    else
        error("No data available and calculation not implemented yet")
    end
end

export numeric_fpdims

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

export numeric_fpdim

numeric_fpdim(fr::FusionRing) = sum(x->x*x, numeric_fpdims(fr))

function is_commutative(fr::FusionRing)
    N = multiplication_table(fr)
    r = rank(fr)
    for a in 1:r, b in 1:r, c in 1:r
        N[a,b,c] == N[b,a,c] || return false
    end
    true
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




# TODO: need to addapt to definition EGNO
# TODO: not sure whether you need a ⊗ b ⊗ a* ⊗ b*, isn't a ⊗ b* enough?
"""
    commutator(fr::FusionRing, A::Vector{Int}, B::Vector{Int}) -> FusionRing

Return  commutator subring `[A,B]`, defined as the smallest fusion-closed subring
containing the support of each product `a ⊗ b ⊗ a* ⊗ b*` with `a ∈ A`, `b ∈ B`.

`A` and `B` are vectors of simple indices (assumed to be subsets of `1:rank(fr)`).
"""
function commutator(fr::FusionRing, A::Vector{Int}, B::Vector{Int})::FusionRing
    r = rank(fr)
    all(1 .≤ A .≤ r) || throw(ArgumentError("commutator: A has out-of-bounds indices"))
    all(1 .≤ B .≤ r) || throw(ArgumentError("commutator: B has out-of-bounds indices"))

    # Seed S0 with the union of supports of a ⊗ b ⊗ a* ⊗ b*
    seen = falses(r)
    @inbounds for a in A
        aᵗ = _dual_index(fr, a)
        for b in B
            bᵗ = _dual_index(fr, b)

            # First multiply a ⊗ b
            for (u, mu) in tensor_product(fr, a, b)
                mu == 0 && continue
                # Then multiply by a* ⊗ b* ; we do it as (u ⊗ a*) ⊗ b*
                for (v, mv) in tensor_product(fr, u, aᵗ)
                    mv == 0 && continue
                    for (w, mw) in tensor_product(fr, v, bᵗ)
                        mw == 0 && continue
                        seen[w] = true
                    end
                end
            end
        end
    end

    S0 = findall(seen)
    isempty(S0) && (S0 = [1])  # at minimum, the unit

    # Close under fusion and build the subring
    S = _fusion_closure(fr, S0)
    _restrict_subring(fr, S; check_closed=true)
end


export categories_with_properties

function categories_with_properties( fr::FusionRing )
    return fr.has_categories_with_props
end

export is_categorifiable

function is_categorifiable( fr::FusionRing )
    return fr.categorifiable
end


# generate all k-subsets of 1:n without external deps - could not find corresponding combinatorics function
# so can be replaced by function once found
function _k_subsets(n::Int, k::Int)
    out = Vector{Vector{Int}}()
    buf = Vector{Int}(undef, k)
    function go(start::Int, depth::Int)
        if depth > k
            push!(out, copy(buf)); return
        end
        # ensure enough remaining
        last = n - (k - depth)
        for v in start:last
            buf[depth] = v
            go(v + 1, depth + 1)
        end
    end
    k == 0 && return [Int[]]
    (k < 0 || k > n) && return out
    go(1,1)
    out
end

#internal-closed subsets of size k (must contain unit=1) 
function _internal_closed_subsets(fr::FusionRing, k::Int)
    r = rank(fr)
    (k == 0 || k > r) && return Vector{Vector{Int}}()
    #  require unit in subset
    candidates = [s for s in _k_subsets(r, k) if first(s) == 1]
    filter(s -> _internal_multiplication(fr, s), candidates)
end

function _internal_multiplication(fr::FusionRing, S::Vector{Int})
    Sset = Set(S)
    @inbounds for i in S, j in S
        for (c, m) in tensor_product(fr, i, j)
            m == 0 && continue
            c in Sset || return false
        end
    end
    true
end

# group by diagonal-channel counts: keys are counts of nonzeros in N[i,i,:] 
function _diag_channel_groups(N::Array{Int,3})
    r = size(N,1)
    key = Vector{Int}(undef, r)
    @inbounds for i in 1:r
        cnt = 0
        for c in 1:r
            N[i,i,c] > 0 && (cnt += 1)
			#can replace with sum 
        end
        key[i] = cnt
    end
    dict = Dict{Int,Vector{Int}}()
    @inbounds for i in 1:r
        push!(get!(dict, key[i], Int[]), i)
    end
    collect(values(dict))
end

# apply permutation P on all three indices: A'[i,j,k] = A[P[i],P[j],P[k]]
function _permute_multtab(A::Array{Int,3}, P::Vector{Int})
    r = size(A,1)
    B = similar(A)
    @inbounds for i in 1:r, j in 1:r, k in 1:r
        B[i,j,k] = A[P[i], P[j], P[k]]
    end
    B
end

function _permutation_vector_equiv(A::Array{Int,3}, B::Array{Int,3})
    r = size(A,1); size(B,1) == r || return nothing
	#think of something else instead of nothing, people might be confused if they run code and nothing appears

    grpA = _diag_channel_groups(A)
    grpB = _diag_channel_groups(B)
    # multiset of group sizes must match
    sort!(map(length, grpA)) == sort!(map(length, grpB)) || return nothing

    used = falses(length(grpB))
    cur  = Vector{Int}(undef, r)

    # enforce unit -> unit if 1 exists (it always should)
    cur[1] = 1

    function backtrack(gidx::Int)
        if gidx > length(grpA)
            return _permute_multtab(A, cur) == B
        end
        GA = grpA[gidx]
        # candidate B-groups of equal size
        for j in eachindex(grpB)
            if used[j] || length(grpB[j]) != length(GA)
                continue
            end
            used[j] = true
            # try all bijections within the group
            for σ in Base.Iterators.permutations(grpB[j])
                # enforce 1 -> 1 if 1 in GA
                if 1 in GA
                    # image of 1 must be 1
                    σ[findfirst(==(1), GA)] == 1 || continue
                end
                # set partial mapping for this group
                for (u, v) in zip(GA, σ)
                    cur[u] = v
                end
                backtrack(gidx + 1) && return true
                # overwrite these entries next loop
            end
            used[j] = false
        end
        return false
    end

    backtrack(1) ? cur : nothing
end



export which_injection
function which_injection(subring::FusionRing, ring::FusionRing)
    rs = rank(subring); rr = rank(ring)
    rs > rr && return nothing
    Nbig = multiplication_table(ring)
    Nsub = multiplication_table(subring)
    for S in _internal_closed_subsets(ring, rs)
        Nres = @views Nbig[S, S, S]
        if (perm = _permutation_vector_equiv(Nsub, Nres)) !== nothing
            inj = Dict{Int,Int}()
            @inbounds for i in 1:rs
                inj[i] = S[perm[i]]
            end
            return inj
        end
    end
    nothing
end

export fusion_ring_automorphisms
function fusion_ring_automorphisms(fr::FusionRing)
    N = multiplication_table(fr)
    r = size(N,1)
    groups = _diag_channel_groups(N)
    perms  = Vector{Vector{Int}}()
    cur    = Vector{Int}(undef, r)

    # enforce unit mapping
    cur[1] = 1

    function backtrack(gidx::Int)
        if gidx > length(groups)
            _permute_multtab(N, cur) == N && push!(perms, copy(cur))
            return
        end
        G = groups[gidx]
        for σ in Base.Iterators.permutations(G)
            # if 1 is in  group,  image must be 1
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








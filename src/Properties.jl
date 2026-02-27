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
    @views multiplication_table(fr)[:, :, 1]
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
                # Then multiply by a* ⊗ b* - do it as (u ⊗ a*) ⊗ b*
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



# TODO 
# 1. Could use unicode to make everything more readable
# 2. Could use dictionaries rather than elseif statements 
# 3. Should put code in Creation.jl

# Anyonica rulesodd[m_]  (metaplectic / SO(m)_2, m odd)
# 
# - Anyonica builds matX as a Table of vectors (length=rank),
#   then does Transpose /@ at the end.
# - i follow that exactly: build "row_table" (rank×rank),
#   then transpose before converting to mt.

# basis vector e_i in ℤ^rank
@inline function _e(i::Int, rank::Int)::Vector{Int}
    v = zeros(Int, rank)
    v[i] = 1
    return v
end

# convert a list of fusion matrices mats[a][b,c] into mt[a,b,c]
function _mats_to_mt(mats::Vector{Matrix{Int}})::Array{Int,3}
    r = length(mats)
    mt = zeros(Int, r, r, r)
    @inbounds for a in 1:r
        mt[a, :, :] .= mats[a]
    end
    return mt
end

function _son2_rules_odd(m::Integer)::Array{Int,3}
    isodd(m) || throw(ArgumentError("_son2_rules_odd expects odd m, got m=$m"))
    m ≥ 5    || throw(ArgumentError("_son2_rules_odd expects m≥5 (odd), got m=$m"))

    r    = (m - 1) ÷ 2
    rank = (m + 7) ÷ 2 

    # convenience
    ar(i) = _e(i, rank)

    # mat1 = IdentityMatrix[rank]
    mat1 = Matrix{Int}(I, rank, rank)

    # matZ = Table[ Which[...], {i,rank} ]
    matZ = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == 1
            ar(2)
        elseif i == 2
            ar(1)
        elseif 3 <= i <= 4
            ar(3 + mod(i, 2))   # i=3 -> 4, i=4 -> 3
        else
            ar(i)
        end
        matZ[i, :] .= v
    end

    # matXe1 = Table[ Which[...], {i,rank} ]
    matXe1 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == 1
            ar(3)
        elseif i == 2
            ar(4)
        elseif i == 3
            # ar[1] + Sum[ar[j], {j, 5, rank}]
            v0 = zeros(Int, rank)
            v0[1] = 1
            for j in 5:rank
                v0[j] += 1
            end
            v0
        elseif i == 4
            # ar[2] + Sum[ar[j], {j, 5, rank}]
            v0 = zeros(Int, rank)
            v0[2] = 1
            for j in 5:rank
                v0[j] += 1
            end
            v0
        else
            ar(3) .+ ar(4)
        end
        matXe1[i, :] .= v
    end

    # matXe2 = Table[ Which[...], {i,rank} ]
    matXe2 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == 1
            ar(4)
        elseif i == 2
            ar(3)
        elseif i == 3
            # ar[2] + Sum[ar[j], {j, 5, rank}]
            v0 = zeros(Int, rank)
            v0[2] = 1
            for j in 5:rank
                v0[j] += 1
            end
            v0
        elseif i == 4
            # ar[1] + Sum[ar[j], {j, 5, rank}]
            v0 = zeros(Int, rank)
            v0[1] = 1
            for j in 5:rank
                v0[j] += 1
            end
            v0
        else
            ar(3) .+ ar(4)
        end
        matXe2[i, :] .= v
    end

    # matY[j_] := Table[ Which[...], {i,rank} ]
    function matY(j::Int)::Matrix{Int}
        M = zeros(Int, rank, rank)
        @inbounds for i in 1:rank
            v = if i == 1
                ar(j + 4)
            elseif i == 2
                ar(j + 4)
            elseif i == 3
                ar(3) .+ ar(4)
            elseif i == 4
                ar(3) .+ ar(4)
            else
                ii = i - 4
                if ii == j
                    # ar[1] + ar[2] + ar[ Min[2j, m-2j] + 4 ]
                    t = min(2*j, m - 2*j) + 4
                    ar(1) .+ ar(2) .+ ar(t)
                else
                    # ar[ Abs[ii-j] + 4 ] + ar[ Min[ii+j, m-ii-j+4] + 4 ]
                    t1 = abs(ii - j) + 4
                    t2 = min(ii + j, m - i - j + 4) + 4
                    ar(t1) .+ ar(t2)
                end
            end
            M[i, :] .= v
        end
        return M
    end

    #   Transpose /@ Join[{mat1, matZ, matXe1, matXe2}, matY /@ Range[r]]
    mats = Matrix{Int}[]
    push!(mats, transpose(mat1))
    push!(mats, transpose(matZ))
    push!(mats, transpose(matXe1))
    push!(mats, transpose(matXe2))
    for j in 1:r
        push!(mats, transpose(matY(j)))
    end

    # Convert fusion matrices to multiplication table
    return _mats_to_mt(mats)
end



#  rulesdiv2[p_] and rulesdiv4[p_]  (SO(m)_2 even cases)
#
# Used by son2_fusion_ring:
#   if m % 4 == 0  => rulesdiv4(m/2)
#   elseif m % 2 == 0 => rulesdiv2(m/2)
#
# 
#   build fusion matrices mats[a] as Integer matrices,
#   apply transpose to match Anyonica's Transpose /@,
#   convert mats -> multiplication table mt[a,b,c].



function _mats_to_mt(mats::Vector{Matrix{Int}})::Array{Int,3}
    r = length(mats)
    mt = zeros(Int, r, r, r)
    @inbounds for a in 1:r
        mt[a, :, :] .= mats[a]
    end
    return mt
end

# index convention (matches  Evaluate[...] = IdentityMatrix[rank]):
# 1: Id
# 2: Θ
# 3: Φ1
# 4: Φ2
# 5: σ1
# 6: σ2
# 7: τ1
# 8: τ2
# 9..rank:  Φ[j] with j = 1..(rank-8)
@inline _Id() = 1
@inline _Th() = 2
@inline _Phi1() = 3
@inline _Phi2() = 4
@inline _Sig1() = 5
@inline _Sig2() = 6
@inline _Tau1() = 7
@inline _Tau2() = 8
@inline _Phi(j::Int) = 8 + j

# rulesdiv2[p_]
function _son2_rules_div2(p::Integer)::Array{Int,3}
    p ≥ 1 || throw(ArgumentError("_son2_rules_div2 expects p≥1, got p=$p"))
    rank = p + 7
    maxphi = rank - 8  # = p-1

    ar(i) = _e(i, rank)

    # sums: sumEvenΛs = Σ_{i=2,4,...,p-1} Φ[i], sumOddΛs = Σ_{i=1,3,...,p-1} Φ[i]
    sumEven = zeros(Int, rank)
    sumOdd  = zeros(Int, rank)
    for i in 1:maxphi
        (isodd(i) ? (sumOdd[_Phi(i)] += 1) : (sumEven[_Phi(i)] += 1))
    end

    mats = Matrix{Int}[]

    # matId = IdentityMatrix[rank]
    matId = Matrix{Int}(I, rank, rank)
    push!(mats, transpose(matId))

    # matΘ
    matTh = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Th())
        elseif i == _Th()
            ar(_Id())
        elseif i == _Phi1()
            ar(_Phi2())
        elseif i == _Phi2()
            ar(_Phi1())
        elseif i == _Sig1()
            ar(_Tau1())
        elseif i == _Sig2()
            ar(_Tau2())
        elseif i == _Tau1()
            ar(_Sig1())
        elseif i == _Tau2()
            ar(_Sig2())
        else
            ar(i) # Φ-lambdas fixed
        end
        matTh[i, :] .= v
    end
    push!(mats, transpose(matTh))

    # matΦ1
    matPhi1 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Phi1())
        elseif i == _Th()
            ar(_Phi2())
        elseif i == _Phi1()
            ar(_Th())
        elseif i == _Phi2()
            ar(_Id())
        elseif i == _Sig1()
            ar(_Sig2())
        elseif i == _Sig2()
            ar(_Tau1())
        elseif i == _Tau1()
            ar(_Tau2())
        elseif i == _Tau2()
            ar(_Sig1())
        else
            # Φ[p - (i-8)]
            j = i - 8
            ar(_Phi(p - j))
        end
        matPhi1[i, :] .= v
    end
    push!(mats, transpose(matPhi1))

    # matΦ2
    matPhi2 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Phi2())
        elseif i == _Th()
            ar(_Phi1())
        elseif i == _Phi1()
            ar(_Id())
        elseif i == _Phi2()
            ar(_Th())
        elseif i == _Sig1()
            ar(_Tau2())
        elseif i == _Sig2()
            ar(_Sig1())
        elseif i == _Tau1()
            ar(_Sig2())
        elseif i == _Tau2()
            ar(_Tau1())
        else
            j = i - 8
            ar(_Phi(p - j))
        end
        matPhi2[i, :] .= v
    end
    push!(mats, transpose(matPhi2))

    # matσ1
    matSig1 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Sig1())
        elseif i == _Th()
            ar(_Tau1())
        elseif i == _Phi1()
            ar(_Sig2())
        elseif i == _Phi2()
            ar(_Tau2())
        elseif i == _Sig1()
            ar(_Phi2()) .+ sumOdd
        elseif i == _Sig2()
            ar(_Id()) .+ sumEven
        elseif i == _Tau1()
            ar(_Phi1()) .+ sumOdd
        elseif i == _Tau2()
            ar(_Th()) .+ sumEven
        else
            # If[OddQ[i], σ2+τ2, σ1+τ1]
            isodd(i) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
        end
        matSig1[i, :] .= v
    end
    push!(mats, transpose(matSig1))

    # matσ2
    matSig2 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Sig2())
        elseif i == _Th()
            ar(_Tau2())
        elseif i == _Phi1()
            ar(_Tau1())
        elseif i == _Phi2()
            ar(_Sig1())
        elseif i == _Sig1()
            ar(_Id()) .+ sumEven
        elseif i == _Sig2()
            ar(_Phi1()) .+ sumOdd
        elseif i == _Tau1()
            ar(_Th()) .+ sumEven
        elseif i == _Tau2()
            ar(_Phi2()) .+ sumOdd
        else
            # If[EvenQ[i], σ2+τ2, σ1+τ1]
            iseven(i) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
        end
        matSig2[i, :] .= v
    end
    push!(mats, transpose(matSig2))

    # matτ1
    matTau1 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Tau1())
        elseif i == _Th()
            ar(_Sig1())
        elseif i == _Phi1()
            ar(_Tau2())
        elseif i == _Phi2()
            ar(_Sig2())
        elseif i == _Sig1()
            ar(_Phi1()) .+ sumOdd
        elseif i == _Sig2()
            ar(_Th()) .+ sumEven
        elseif i == _Tau1()
            ar(_Phi2()) .+ sumOdd
        elseif i == _Tau2()
            ar(_Id()) .+ sumEven
        else
            isodd(i) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
        end
        matTau1[i, :] .= v
    end
    push!(mats, transpose(matTau1))

    # matτ2
    matTau2 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Tau2())
        elseif i == _Th()
            ar(_Sig2())
        elseif i == _Phi1()
            ar(_Sig1())
        elseif i == _Phi2()
            ar(_Tau1())
        elseif i == _Sig1()
            ar(_Th()) .+ sumEven
        elseif i == _Sig2()
            ar(_Phi2()) .+ sumOdd
        elseif i == _Tau1()
            ar(_Id()) .+ sumEven
        elseif i == _Tau2()
            ar(_Phi1()) .+ sumOdd
        else
            iseven(i) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
        end
        matTau2[i, :] .= v
    end
    push!(mats, transpose(matTau2))

    # matΦ[j] for j = 1..(rank-8)
    function matPhi(j::Int)::Matrix{Int}
        M = zeros(Int, rank, rank)
        @inbounds for i in 1:rank
            v = if i == _Id()
                ar(_Phi(j))
            elseif i == _Th()
                ar(_Phi(j))
            elseif i == _Phi1()
                ar(_Phi(p - j))
            elseif i == _Phi2()
                ar(_Phi(p - j))
            elseif i == _Sig1()
                isodd(j) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
            elseif i == _Sig2()
                iseven(j) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
            elseif i == _Tau1()
                isodd(j) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
            elseif i == _Tau2()
                iseven(j) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
            else
                ii = i - 8
                if ii == j && (2*j < p)
                    ar(_Id()) .+ ar(_Th()) .+ ar(_Phi(2*j))
                elseif ii == j && (2*j > p)
                    ar(_Id()) .+ ar(_Th()) .+ ar(_Phi(2*(p - j)))
                elseif ii + j < p
                    ar(_Phi(abs(ii - j))) .+ ar(_Phi(ii + j))
                elseif ii + j > p
                    ar(_Phi(abs(ii - j))) .+ ar(_Phi(2*p - ii - j))
                else
                    # ii == p - j
                    ar(_Phi1()) .+ ar(_Phi2()) .+ ar(_Phi(abs(p - 2*ii)))
                end
            end
            M[i, :] .= v
        end
        return M
    end

    for j in 1:maxphi
        push!(mats, transpose(matPhi(j)))
    end

    return _mats_to_mt(mats)
end


# rulesdiv4[p_]
function _son2_rules_div4(p::Integer)::Array{Int,3}
    p ≥ 1 || throw(ArgumentError("_son2_rules_div4 expects p≥1, got p=$p"))
    rank = p + 7
    maxphi = rank - 8  # = p-1

    ar(i) = _e(i, rank)

    sumEven = zeros(Int, rank)
    sumOdd  = zeros(Int, rank)
    for i in 1:maxphi
        (isodd(i) ? (sumOdd[_Phi(i)] += 1) : (sumEven[_Phi(i)] += 1))
    end

    mats = Matrix{Int}[]

    matId = Matrix{Int}(I, rank, rank)
    push!(mats, transpose(matId))

    # matΘ
    matTh = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Th())
        elseif i == _Th()
            ar(_Id())
        elseif i == _Phi1()
            ar(_Phi2())
        elseif i == _Phi2()
            ar(_Phi1())
        elseif i == _Sig1()
            ar(_Tau1())
        elseif i == _Sig2()
            ar(_Tau2())
        elseif i == _Tau1()
            ar(_Sig1())
        elseif i == _Tau2()
            ar(_Sig2())
        else
            ar(i)
        end
        matTh[i, :] .= v
    end
    push!(mats, transpose(matTh))

    # matΦ1
    matPhi1 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Phi1())
        elseif i == _Th()
            ar(_Phi2())
        elseif i == _Phi1()
            ar(_Id())
        elseif i == _Phi2()
            ar(_Th())
        elseif i == _Sig1()
            ar(_Tau1())
        elseif i == _Sig2()
            ar(_Sig2())
        elseif i == _Tau1()
            ar(_Sig1())
        elseif i == _Tau2()
            ar(_Tau2())
        else
            j = i - 8
            ar(_Phi(p - j))
        end
        matPhi1[i, :] .= v
    end
    push!(mats, transpose(matPhi1))

    # matΦ2
    matPhi2 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Phi2())
        elseif i == _Th()
            ar(_Phi1())
        elseif i == _Phi1()
            ar(_Th())
        elseif i == _Phi2()
            ar(_Id())
        elseif i == _Sig1()
            ar(_Sig1())
        elseif i == _Sig2()
            ar(_Tau2())
        elseif i == _Tau1()
            ar(_Tau1())
        elseif i == _Tau2()
            ar(_Sig2())
        else
            j = i - 8
            ar(_Phi(p - j))
        end
        matPhi2[i, :] .= v
    end
    push!(mats, transpose(matPhi2))

    # matσ1
    matSig1 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Sig1())
        elseif i == _Th()
            ar(_Tau1())
        elseif i == _Phi1()
            ar(_Tau1())
        elseif i == _Phi2()
            ar(_Sig1())
        elseif i == _Sig1()
            ar(_Id()) .+ ar(_Phi2()) .+ sumEven
        elseif i == _Sig2()
            sumOdd
        elseif i == _Tau1()
            ar(_Th()) .+ ar(_Phi1()) .+ sumEven
        elseif i == _Tau2()
            sumOdd
        else
            isodd(i) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
        end
        matSig1[i, :] .= v
    end
    push!(mats, transpose(matSig1))

    # matσ2
    matSig2 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Sig2())
        elseif i == _Th()
            ar(_Tau2())
        elseif i == _Phi1()
            ar(_Sig2())
        elseif i == _Phi2()
            ar(_Tau2())
        elseif i == _Sig1()
            sumOdd
        elseif i == _Sig2()
            ar(_Id()) .+ ar(_Phi1()) .+ sumEven
        elseif i == _Tau1()
            sumOdd
        elseif i == _Tau2()
            ar(_Th()) .+ ar(_Phi2()) .+ sumEven
        else
            iseven(i) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
        end
        matSig2[i, :] .= v
    end
    push!(mats, transpose(matSig2))

    # matτ1
    matTau1 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Tau1())
        elseif i == _Th()
            ar(_Sig1())
        elseif i == _Phi1()
            ar(_Sig1())
        elseif i == _Phi2()
            ar(_Tau1())
        elseif i == _Sig1()
            ar(_Th()) .+ ar(_Phi1()) .+ sumEven
        elseif i == _Sig2()
            sumOdd
        elseif i == _Tau1()
            ar(_Id()) .+ ar(_Phi2()) .+ sumEven
        elseif i == _Tau2()
            sumOdd
        else
            isodd(i) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
        end
        matTau1[i, :] .= v
    end
    push!(mats, transpose(matTau1))

    # matτ2
    matTau2 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Tau2())
        elseif i == _Th()
            ar(_Sig2())
        elseif i == _Phi1()
            ar(_Tau2())
        elseif i == _Phi2()
            ar(_Sig2())
        elseif i == _Sig1()
            sumOdd
        elseif i == _Sig2()
            ar(_Th()) .+ ar(_Phi2()) .+ sumEven
        elseif i == _Tau1()
            sumOdd
        elseif i == _Tau2()
            ar(_Id()) .+ ar(_Phi1()) .+ sumEven
        else
            iseven(i) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
        end
        matTau2[i, :] .= v
    end
    push!(mats, transpose(matTau2))

    # matΦ[j]
    function matPhi(j::Int)::Matrix{Int}
        M = zeros(Int, rank, rank)
        @inbounds for i in 1:rank
            v = if i == _Id()
                ar(_Phi(j))
            elseif i == _Th()
                ar(_Phi(j))
            elseif i == _Phi1()
                ar(_Phi(p - j))
            elseif i == _Phi2()
                ar(_Phi(p - j))
            elseif i == _Sig1()
                isodd(j) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
            elseif i == _Sig2()
                iseven(j) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
            elseif i == _Tau1()
                isodd(j) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
            elseif i == _Tau2()
                iseven(j) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
            else
                ii = i - 8
                if ii == j && (2*j < p)
                    ar(_Id()) .+ ar(_Th()) .+ ar(_Phi(2*j))
                elseif ii == j && (2*j == p)
                    ar(_Id()) .+ ar(_Th()) .+ ar(_Phi1()) .+ ar(_Phi2())
                elseif ii == j && (2*j > p)
                    ar(_Id()) .+ ar(_Th()) .+ ar(_Phi(2*(p - j)))
                elseif ii + j < p
                    ar(_Phi(abs(ii - j))) .+ ar(_Phi(ii + j))
                elseif ii + j > p
                    ar(_Phi(abs(ii - j))) .+ ar(_Phi(2*p - ii - j))
                else
                    # ii == p - j
                    ar(_Phi1()) .+ ar(_Phi2()) .+ ar(_Phi(abs(p - 2*ii)))
                end
            end
            M[i, :] .= v
        end
        return M
    end

    for j in 1:maxphi
        push!(mats, transpose(matPhi(j)))
    end

    return _mats_to_mt(mats)
end




# SO(m)_2 / Metaplectic(m) t)

# Uses:
#   _son2_rules_odd(m)      # for odd m
#   _son2_rules_div2(p)     # for m ≡ 2 (mod 4), with p = m÷2
#   _son2_rules_div4(p)     # for m ≡ 0 (mod 4), with p = m÷2

export son2_fusion_ring


# odd m: rank = (m+7)/2, elements are [1, Z, X_e1, X_e2, Y_1, ..., Y_r], r=(m-1)/2
function _son2_labels_odd(m::Int)::Vector{String}
    r = (m - 1) ÷ 2
    labels = String["1", "Z", "Xₑ₁", "Xₑ₂"]
    for j in 1:r
        push!(labels, "Y_$j")
    end
    return labels
end

# even m: rank = p+7 with p=m/2, elements are [Id, Θ, Φ1, Φ2, σ1, σ2, τ1, τ2, Φ_1..Φ_{p-1}]
function _son2_labels_even(p::Int)::Vector{String}
    labels = String[
        "1", "Θ", "Φ₁", "Φ₂", "σ₁", "σ₂", "τ₁", "τ₂"
    ]
    for j in 1:(p - 1)
        push!(labels, "Φ_$j")
    end
    return labels
end


"""
    son2_fusion_ring(m::Int) -> FusionRing

Return fusion ring SO(N)_2 (metaplectic) .
"""

#- odd `N`: uses `_son2_rules_odd(m)`
#- even `N ≡ 0 (mod 4)`: uses `_son2_rules_div4(N÷2)`
#- even `N ≡ 2 (mod 4)`: uses `_son2_rules_div2(N÷2)`

function son2_fusion_ring(m::Int)::FusionRing
    m ≥ 4 || throw(ArgumentError("son2_fusion_ring(m): requires integer m ≥ 4, got m=$m"))

    mt::Array{Int,3}
    labels::Vector{String}

    if isodd(m)
        mt = _son2_rules_odd(m)
        labels = _son2_labels_odd(m)
    else
        p = m ÷ 2
        if m % 4 == 0
            mt = _son2_rules_div4(p)
        else
            mt = _son2_rules_div2(p)
        end
        labels = _son2_labels_even(p)
    end

    #  label count must match rank
    size(mt, 1) == length(labels) || error("son2_fusion_ring: label length mismatch with mt rank")

    R = fusion_ring(
        mt;
        names  = ["SO($m)_2", "Metaplectic($m)"],
        labels = labels,
    )

    return replace_by_known(R)
end





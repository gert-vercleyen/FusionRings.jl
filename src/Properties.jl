
function change_fusion_ring_property(r::FusionRing, dict)::FusionRing

end

export multiplication_table

function multiplication_table(r::FusionRing)::Array{Int,3}
  return r.multiplication_table
end

export print_multiplication_table

function print_multiplication_table(r::FusionRing)
  rk = rank(r)
  mt = multiplication_table(r)

  tab = fill( "", rk, rk )
  for i in 1:rk, j in 1:rk
    tab[i,j] = row_to_string(r,mt[i,j,:]) 
  end
  tab
end

export row_to_string

function row_to_string(r::FusionRing, row)::String
  n             = length(row)
  el_names      = element_names(r)
  non_zero_ind  = findall(i -> row[i] > 0, 1:n)
  to_string(i)  = element_to_string(row[i], el_names[i])

  join( 
    map(to_string, non_zero_ind), 
    " ⊕ "
  )
end

#check this out
function tensor_product(fr::FusionRing, a, b)
    imap = indexmap(fr)
    normalize(x) = x isa Integer ? x : imap[String(x)]
    ai = normalize(a); bi = normalize(b)
    N = fusion_tensor(fr)[ai,bi,:]
    out = Dict{String,Int}()
    L = labels(fr)
    for (ci,m) in enumerate(N)
        m==0 && continue
        out[L[ci]] = m
    end
    out
end

function element_to_string(mult,elem)::String
  if mult == 0 
    return ""
  elseif mult == 1
    return elem
  else 
    return string(mult) * " " * elem 
  end
end

pmt = print_multiplication_table

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

export element_names

function element_names(r::FusionRing)::Array{String,1}
  return r.element_names
end

export conjugation_matrix

function conjugation_matrix(r::FusionRing)::Array{Int,2}
  return multiplication_table(r)[:,:,1]
end

export is_commutative

function is_commutative(r::FusionRing)::Bool
  mt = multiplication_table(r)
  rk = rank(r)
  all( 
    mat -> mat == mat', 
    [ mt[:,:,i] for i in 1:rk ]
  )
end

export multiplicity

function multiplicity(r::FusionRing)::Int
  maximum(multiplication_table(r))
end

export nonzero_structure_constants

function nonzero_structure_constants(r::FusionRing)::Vector{Tuple{Int64, Int64, Int64}}
  mt = multiplication_table(r)
  map( Tuple, findall( x -> x > 0, mt ) ) 
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
  return r.formal_code
end

export barcode 

function barcode(r::FusionRing)::Int
  return r.barcode
end

function mult_tab_code(mat::Array{Int,2},mult::Int)::Int
end

export sub_fusion_rings

function sub_fusion_rings(r::FusionRing)
  return r.sub_fusion_rings
end

function sub_ring_tables(mat::Array{Int,2})

end

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

function decompositions(r::FusionRing,product="TensorProduct")::Array{FusionRing,1}
  if product == "TensorProduct"
    return r.tensor_product_decompositions
  else 
    return error("Only tensor product decompositions are defined at the moment.")
  end
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

# TODO: Should not be done for non-commutative rings!!!
function characters(ring::FusionRing)
  if !(ring.characters === missing)
    return ring.characters
  elseif !FusionRings.is_commutative(ring) 
    error("Calculation of characters for non-commutative fusion ring is not implemented yet.")
  else
    qqb  = algebraic_closure(QQ) 
    mt   = FusionRings.multiplication_table( ring )
    r    = FusionRings.rank(ring)
    mats = [ matrix( qqb, mt[ i, :, : ] ) for i ∈ 1:r ]

    function is_character_table( mat, mats )
      all( is_diagonal( mat * m * inv(mat) ) for m in mats )
    end
    
    charsq = false
    upi = 9
    upj = 9
    
    proposedchars = mats[1]
    while !charsq
      upi += 1
      upj += 1
      # Take random linear rational combination of fusion mats
      rvec 		= rand( [ i//j for i ∈ 1:upi, j ∈ 1:upj ], r )
      combinedmat = rvec[1] * mats[1]
      for i ∈ 2:r
        combinedmat += rvec[i] * mats[i]
      end

      # Find diagonalizing matrix
      proposedchars = generalized_jordan_form( combinedmat )[2]
      charsq = is_character_table( proposedchars, mats )
    end
    proposedchars
  end 
end

function is_character_table( mat, mats )
  all( is_diagonal( mat * m * inv(mat) ) for m in mats )
end

function is_character_table( mat, ring )
	mt   = FusionRings.multiplication_table( ring )
	r    = FusionRings.rank(ring)
	mats = [ matrix( qqb, mt[ i, :, : ] ) for i ∈ 1:r ]

  all( is_diagonal( mat * m * inv(mat) ) for m in mats )
end

export modular_data

function modular_data(r::FusionRing)
  return r.modular_data
end

function s_matrices(r::FusionRing)

end
  
function normalized_s_matrices(r::FusionRing)
  
end

function twist_factors(r::FusionRing)

end


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
    commutator(fr::FusionRing) -> FusionRing

Return  "derived" fusion subring obtained as  fusion-closure of the supports of
`a ⊗ b ⊗ a* ⊗ b*` for all simples `a,b` in `fr`. This is  commutator of the full ring
with itself
"""
function commutator(fr::FusionRing)::FusionRing
    r = rank(fr)
    return commutator(fr, collect(1:r), collect(1:r))
end

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



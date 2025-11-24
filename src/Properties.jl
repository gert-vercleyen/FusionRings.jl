
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
  el_names      = labels(r)
  non_zero_ind  = findall(i -> row[i] > 0, 1:n)
  to_string(i)  = element_to_string(row[i], el_names[i])

  join( 
    map(to_string, non_zero_ind), 
    " ⊕ "
  )
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
  return r.labels
end

export labels

function labels(r::FusionRing)::Array{String,1}
  return r.labels
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
      rvec = rand( [ i//j for i ∈ 1:upi, j ∈ 1:upj ], r )
      sgnvec = rand( [ -1 1 ], r )
      combinedmat = sgnvec[1] * rvec[1] * mats[1]
      for i ∈ 2:r
        combinedmat += sgnvec[i] * rvec[i] * mats[i]
      end

      # Find diagonalizing matrix
      proposedchars = generalized_jordan_form( combinedmat )[2]
      charsq = is_character_table( proposedchars, mats )
    end
    # TODO: Normalize matrix and set up some convention for order of characters
    [ proposedchars[i,j] for i in 1:r, j in 1:r ]
  end 
end

function is_character_table( mat, ring::FusionRing )
	mt   = FusionRings.multiplication_table( ring )
	r    = FusionRings.rank(ring)
	mats = [ matrix( qqb, mt[ i, :, : ] ) for i ∈ 1:r ]

  all( is_diagonal( mat * m * inv(mat) ) for m in mats )
end

function to_canonical_character_tab()
end

function normalize_characters()
end



"""
    numeric_characters(R::FusionRing; tries=8, tol=1e-10) -> (C, V)

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
function numeric_characters(R::FusionRing; tries::Int=8, tol::Real=1e-10)
    labs = labels(R)
    r = length(labs)
    Nis = [Matrix{Float64}(fusion_matrix(R, a)) for a in labs]

    # quick commutativity sanity check
    if !FusionRings.is_commutative(R) 
        error("fusion_ring_characters: ring appears non-commutative; this routine requires commuting fusion matrices.")
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

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

function characters(r::FusionRing)
  return r.characters
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
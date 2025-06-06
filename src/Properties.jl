
function change_fusion_ring_property(r::FusionRing, dict)::FusionRing

end

export multiplication_table

function multiplication_table(r::FusionRing)::Array{Int,3}
  return r.multiplication_table
end

function print_multiplication_table(r::FusionRing)
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
  rank = rank(r)
  all( 
    mat -> mat == mat', 
    [ mt[:,:,i] for i in 1:rank ]
  )
end

export multiplicity

function multiplicity(r::FusionRing)::Int
  maximum(multiplication_table(r))
end

function nonzero_structure_constants(r::FusionRing)::Array{Int,2}

end

nzsc = nonzero_structure_constants

function num_nonzero_structure_constants(r::FusionRing)::Array{Int,2}

end

nnzsc = num_nonzero_structure_constants

export frobenius_perron_dimensions

function frobenius_perron_dimensions(r::FusionRing)
  r.frobenius_perron_dimensions
end

export fpdims 

fpdims = frobenius_perron_dimensions

export frobenius_perron_dimension

function frobenius_perron_dimension(r::FusionRing)
  return sum( fpdims(r).^2 )
end

export fpdim

fpdim = frobenius_perron_dimension

function num_self_dual_non_self_dual(r::FusionRing)::Array{Int,1}

end

nsdnsd = num_self_dual_non_self_dual

function num_self_dual(r::FusionRing)::Int
  
end

nsd = num_self_dual

function num_non_self_dual(r::FusionRing)::Int
  
end

export is_group_ring

function is_group_ring(r::FusionRing)::Bool
  return sum(multiplication_table(r)) == rank(r)^2
end

function conjugate_element(r::FusionRing)

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

# TODO: This should also work for rings for which no information is 
# known yet! 
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
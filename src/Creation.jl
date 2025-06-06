include("GeneralFunctions.jl")

struct FusionRing 
  multiplication_table::Array{Int,3}
  names
  texnames
  element_names::Array{String,1}
  barcode
  formal_code
  direct_product_decompositions
  sub_fusion_rings
  Frobenius_Perron_dimensions
  modular_data
  characters
end 

export fusion_ring

function fusion_ring( 
  mt::Array{Int,3};
  names = missing,
  texnames = missing,
  element_names = missing,
  barcode = missing,
  formal_code = missing,
  tensor_product_decompositions = missing,
  sub_fusion_rings = missing, 
  frobenius_perron_dimensions = missing,
  modular_data = missing,
  characters = missing
  )
  if !check_struct_const(mt)
    return error("All structure constants should be positive integers")
  elseif !check_mt_dims(mt)
    return error("The multiplication table should be a 3D array with consistent depth.")
  elseif !check_unit(mt)
    return error("The multiplication table should have the unit as its first element.")
  elseif !check_inverse(mt)
    return error("At least one of the elements has either no inverse or multiple inverses.")
  elseif !check_associativity(mt)
    return error("The multiplication table does not correspond to associative multiplication.")
  elseif !(element_names===missing) && !check_element_names(mt,element_names)
    return error("The number of element names provides does not match the number of elements.")
  end
  if element_names===missing
    element_names = 
  end
  FusionRing(
    mt,
    names,
    texnames,
    element_names,
    barcode, 
    formal_code,
    tensor_product_decompositions,
    sub_fusion_rings, 
    frobenius_perron_dimensions,
    modular_data,
    characters
  )
end

function check_struct_const(mt)
  all( x -> x >= 0, mt )
end 

function check_mt_dims(mt)
  dims = size(mt)
  dim1 = dims[1];
  length(dims) === 3 && is_constant_array(dims) 
end

function check_unit(mt)
  mt[1,:,:] == mt[:,1,:] == I
end


function check_inverse(mt)
  rank = size(mt)[1]
  sum(mt[:,:,1]) == rank
end


function check_associativity(mt::Array{Int,3})
  assoc = true
  rank = size(mt)[1]
  sum1(a::Int,b::Int,c::Int,d::Int) = sum([mt[a,b,e]*mt[e,c,d] for e in 1:rank])
  sum2(a::Int,b::Int,c::Int,d::Int) = sum([mt[a,f,d]*mt[b,c,f] for f in 1:rank])
  for a in 1:rank, b in 1:rank, c in 1:rank, d in 1:rank
    if sum1(a,b,c,d) != sum2(a,b,c,d)
      assoc = false
      break 
    else 
      continue
    end
  end
  return assoc
end


function check_element_names(mt,names)
  rank = size(mt)[1]
  length(names) == rank
end

# Generators of special fusion rings 

# PSU(2)ₖ
function psu2k_fusion_ring(k::Int)::FusionRing

end

# SU(2)ₖ
function su2k_fusion_ring(k::Int)::FusionRing

end

function son2_fusion_ring(n::Int)::FusionRing
  
end

function metaplectic_fusion_ring(m)::FusionRing
  
end

# Fusion ring from group
function fusion_ring_from_group(grp)::FusionRing
  
end

# Fusion ring from group mult tab
function fusion_ring_from_group(gmt::Array{Int,2}; skipcheck = false)::FusionRing
  
end

function groupname(grp)::String
  
end

# Zₙ
function zn_fusion_ring(n::Int)::FusionRing
  
end

# Rep(G)
function group_rep_fusion_ring(grp)::FusionRing
  
end

# Haagerup_Izumi
function hi_fusion_ring(grp)::FusionRing
  
end

# Tambara Yamagami
function ty_fusion_ring(grp)::FusionRing
  
end


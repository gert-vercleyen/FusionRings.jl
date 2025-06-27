include("GeneralFunctions.jl")

struct FusionRing 
  multiplication_table
  names
  texnames
  element_names::Array{String,1}
  barcode
  formal_code
  direct_product_decompositions
  sub_fusion_rings
  frobenius_perron_dimensions
  modular_data
  characters
end 

export fusion_ring

function fusion_ring( 
  mt;
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
    element_names = [ bold_integer(i) for i in 1:first(size(mt)) ]
  end
  FR = FusionRing(
    ZZ.(mt),
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
  return FR
end

function check_struct_const(mt)
  all( x -> x >= 0, mt )
end 

function check_mt_dims(mt)
  dims = size(mt)
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

export psu2k_fusion_ring

# PSU(2)ₖ
function psu2k_fusion_ring(k::Int)::FusionRing
  rk = div(k,2) + 1
  mt = fill( 0, rk, rk, rk )
  for a = 0:2:k, b = 0:2:k, c = 0:2:k
    if c ∈ range_psu2k( a, b, k )
      mt[div(a,2)+1,div(b,2)+1,div(c,2)+1] = 1
    else
      continue
    end
  end
  ring = 
    fusion_ring(
      mt, 
      names = [ "PSU(2)" * subscript_integer(k) ],
      element_names = [ to_psu2_elname((i-1)//2) for i in 1:rk ]
      # TODO: use general formula for fpdims
    )

  return ring
end

function to_psu2_elname(a::Rational) 
  if denominator(a) == 1
    return string(numerator(a))
  else 
    return string(numerator(a)) * "/" * string(denominator(a))
  end
end

range_psu2k(i,j,k) = abs(i-j):2:minimum([i+j, 2k - i - j])

export su2k_fusion_ring

# SU(2)ₖ
function su2k_fusion_ring(k::Int)::FusionRing
  rk = k+1
  mt = fill( 0, rk, rk, rk )
  for a = 0:k, b = 0:k, c = 0:k
    if c ∈ range_psu2k( a, b, k )
      mt[a+1,b+1,c+1] = 1
    else
      continue
    end
  end
  ring = 
    fusion_ring(
      mt, 
      names = [ "SU(2)" * subscript_integer(k) ]
      # TODO: use general formula for fpdims
      # TODO: use standard notation for spin reps
    )

  return ring
end

function son2_fusion_ring(n::Int)::FusionRing
  
end

function metaplectic_fusion_ring(m)::FusionRing
  
end

# Fusion ring from group
function fusion_ring_from_group(grp)::FusionRing
  
end

function is_cayley_table(gmt::Array{Int,2})
  # if 
  # "Not every row and column in `1` contains numbers from 1 to `2`.";
  # "Not every element in `1` has a unique inverse.";
  # "The multiplication defined by `1` is not associative.";
  # "Not every element is invertible."; 
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


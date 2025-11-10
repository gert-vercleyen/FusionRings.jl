# ┌────────────────────────────────────────────────────────────────────────────┐
# │                              Pretty printing                               │
# └────────────────────────────────────────────────────────────────────────────┘



function transform_integer( i::Int, dict::Dict ) 
  dgts = digits( i, base = 10 )
  join( reverse( [ dict[j] for j in dgts  ] ) )
end

bold_integer(i::Int)::String = transform_integer(i, bold_digits_dict)

bold_digits_dict = 
  Dict(
    0 => "𝟎", 1 => "𝟏", 2 => "𝟐", 3 => "𝟑", 4 => "𝟒", 
    5 => "𝟓", 6 => "𝟔", 7 => "𝟕", 8 => "𝟖", 9 => "𝟗"
  )

subscript_integer(i::Int)::String = transform_integer(i,subs_digits_dict)

subs_digits_dict = 
  Dict(
    0 => "₀", 1 => "₁", 2 => "₂", 3 => "₃", 4 => "₄", 
    5 => "₅", 6 => "₆", 7 => "₇", 8 => "₈", 9 => "₉"
  )

superscript_integer(i::Int) = transform_integer(i,sup_digits_dict)

sup_digits_dict = 
  Dict(
    0 => "⁰", 1 => "¹", 2 => "²", 3 => "³", 4 => "⁴", 
    5 => "⁵", 6 => "⁶", 7 => "⁷", 8 => "⁸", 9 => "⁹"
  )

function is_constant_array( arr; equalfunc = === ) 
  if isempty(arr)
    return true 
  end
  first = arr[1]
  return all( equalfunc( element, first ) for element in arr )
end

# comap applies each function in an array to a single argument
function comap( arr, arg )
  [ f(arg) for f in arr ]
end

export to_numberfield

function to_composite_field( 
    arr::Array{QQBarFieldElem}; 
    simplify_field = true, 
    canonical_simplification = true
    )

  # we don't store the field since its available in 
  # the parent field of each number. The field is Abstract
  # so we do need to return the embedding
  f = number_field( QQ, unique( arr ) )[2]
   
  if simplify_field 
    g = simplify( numfield; canonical = canonical_simplification )[2]
    to_field_elem = x -> preimage( g, preimage( f, x ) )
    return ( to_field_elem.(arr), ( g, f ) ) #TODO: this should be the composition
  else 
    to_field_elem = x ->  preimage( f, x )
    return ( to_field_elem.(arr), f )
  end
end
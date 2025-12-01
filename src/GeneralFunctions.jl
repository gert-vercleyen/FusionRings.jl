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

export to_composite_field

function to_composite_field( 
  arr::Array{QQBarFieldElem}; 
  simplify_field = false, 
  canonical_simplification = true
  )

  K, f = number_field( QQ, unique( arr ) )
  
  if simplify_field 
    L, g = simplify( K; canonical = canonical_simplification )
    to_field_elem  = x -> preimage( g, preimage( f, x ) )
    fg = hom( L, algebraic_closure(QQ), f(g(gen(L))) )
    return ( to_field_elem.(arr), fg )
  else 
    to_field_elem = x -> preimage( f, x )
    return ( to_field_elem.(arr), f )
  end
end

function to_cyclotomic_field( arr::Array{AbsSimpleNumFieldElem}, emb ) 
	length(arr) === 0 && return ( arr, emb )
	
	# Check parrent field of all fields are equal
	is_constant_array( parent.( arr ) ) || error("Elements of array should belong to same field")

	qqb = algebraic_closure(QQ)
	K   = parent( arr[1] )
	C   = ray_class_field( K ) 
	deg = C |> conductor |> first |> minimum |> Int
	L,  = cyclotomic_field( deg )

	gen_K_as_cyclo = first( roots( L, defining_polynomial(K) ) )
	to_cyclo       = hom( K, L, gen_K_as_cyclo )

	for j in 1:deg
		emb_cyclo = hom( L, qqb, roots( qqb, defining_polynomial(L) )[j] )

		if emb_cyclo(gen_K_as_cyclo) == emb(gen(K)) 	
			return ( to_cyclo.(arr), emb_cyclo )
		else
			continue
		end
	end

	error("Couldn't find embedding from cyclotomics into algebraic_closure(QQ)")
end
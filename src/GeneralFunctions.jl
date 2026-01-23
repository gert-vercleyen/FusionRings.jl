
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
    x::QQBarFieldElem,
    simplify_field = false, 
    canonical_simplification = true
    )
    arr, emb = to_composite_field(
        [x],
        simplify_field,
        canonical_simplification
    )
    ( arr[1], emb )
end

function to_composite_field( 
  arr::Array{QQBarFieldElem}; 
  simplify_field = false, 
  canonical_simplification = true
  )

  K, f = number_field( QQ, unique( arr ) )
  
  if simplify_field 
    L, g = simplify( K; canonical = canonical_simplification )
    to_field_elem  = x -> preimage( g, preimage( f, x ) )
    fg = hom( L, algebraic_closure(QQ), (f ∘ g ∘ gen)(L) )
    return ( to_field_elem.(arr), fg )
  else 
    to_field_elem = x -> preimage( f, x )
    return ( to_field_elem.(arr), f )
  end
end

export to_cyclotomic_field


function to_cyclotomic_field(
    x::QQBarFieldElem,
    simplify_field = false, 
    canonical_simplification = true
    )
    cfx, emb =
        to_composite_field(
            x,
            simplify_field,
            canonical_simplification
        )

    to_cyclotomic_field( el, emb )
end

function to_cyclotomic_field( x::AbsSimpleNumFieldElem, emb )
    arr, emb, deg = to_cyclotomic_field( [x], emb )

    ( arr[1], emb, deg )
end

function to_cyclotomic_field( arr::Array{AbsSimpleNumFieldElem}, emb ) 
	length(arr) === 0 && return (arr, emb)
	
	# Check parrent field of all fields are equal
	is_constant_array( parent.(arr) ) || error("Elements of array should belong to same field")

	qqb = algebraic_closure(QQ)
	K   = parent( arr[1] )
	C   = ray_class_field(K) 
	deg = (Int ∘ minimum ∘ first ∘ conductor)(C)
	L,  = cyclotomic_field(deg)

	gen_K_as_cyclo = first( roots( L, defining_polynomial(K) ) )
	to_cyclo       = hom( K, L, gen_K_as_cyclo )

    rts = roots( qqb, defining_polynomial(L) )
	for j in 1:deg
		emb_cyclo = hom( L, qqb, rts[j] )

		if emb_cyclo(gen_K_as_cyclo) == emb( gen(K) ) 	
			return ( to_cyclo.(arr), emb_cyclo, deg )
		else
			continue
		end
	end

	error("Couldn't find embedding from cyclotomics into algebraic_closure(QQ)")
end

# Returns element of v who's value equals x
# Super inneficient implementation at the moment since
# we just loop over the list v
function replace_by_known( v; tol=1e-10 )
    function (x)
        CC = AcbField(64);
    	conv(z) = convert(ComplexF64,z)
        for y in v
            abs(conv(x) - conv(CC(y))) < tol && return y
        end
        error("No matching value found")
    end
end

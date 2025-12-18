### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 4f8788dc-b431-11f0-01df-fd75003adb9b
# ╠═╡ show_logs = false
begin
	using Revise
	using Pkg;
	Pkg.develop(path="/home/gert/Projects/FusionRings.jl/")
	#Pkg.develop(path="/Users/gertvercleyen/Projects/FusionRings.jl/")
	using FusionRings
	using Oscar
	using JSON
	using Base.Threads
end

# ╔═╡ 95eeb20c-3efd-4285-b6a3-8d8d37aaee05
md"""
# Initialize package
"""

# ╔═╡ 5b6a60f2-d925-49ed-a563-f8bc630f43c0
md"""
# Load fusion ring list
"""

# ╔═╡ cbaa67dc-35db-4759-aaab-09e461581c16
function frl(i::Int) 
  js = JSON.parsefile("/home/gert/Tests/JSONExport/ring_"*string(i)*".json")
  fc = [ js["formal_code"][i] for i in 1:4 ]
  r = fc[1]
  mt = zeros(Int, r, r, r)
  for i in 1:r, j in 1:r, k in 1:r 
      mt[i,j,k] = Int.(js["mt"][i][j][k])
  end
  FusionRings.fusion_ring( mt, formal_code = fc)
end

# ╔═╡ ac3382fa-1341-4017-ae7a-b96bdfa0ea67
@time FRL = [ frl(i) for i in 1:28451 ];

# ╔═╡ e6ef9a18-94a6-44e5-b590-f08a12b2f4cf
md"""
# Initialize ``\overline{\mathbb{Q}}``
"""

# ╔═╡ 69d5a5bf-6952-4870-bbfb-39170b17940c
qqb = algebraic_closure(QQ)

# ╔═╡ 2921df86-567e-4149-9677-c6fc821cb210
md"""
# Some useful definitions 
"""

# ╔═╡ 93079240-ee00-4e8b-8bcd-ca98d1c7f016
mt = FusionRings.multiplication_table

# ╔═╡ 9e9d0d16-5907-4601-b59a-f03556541567
function is_constant_array( arr; equalfunc = === ) 
  if isempty(arr)
    return true 
  end
  first = arr[1]
  return all( equalfunc( element, first ) for element in arr )
end

# ╔═╡ 3d817d5e-728a-45b3-bdba-7270a72a7c9b
md"""
# Time to test 
"""

# ╔═╡ b25af0b1-010b-4a9d-8698-dab097f90ed7
# ╠═╡ disabled = true
#=╠═╡
function to_combined_numberfield( 
    arr::Array{QQBarFieldElem}; 
    simplify_field = true, 
    canonical_simplification = true
    )

  K, f = number_field( QQ, unique( arr ), cached = false )
   
  if simplify_field 
    L, g = simplify( K; canonical = canonical_simplification )
    to_field_elem  = x -> preimage( g, preimage( f, x ) )
	fg = hom( L, algebraic_closure(QQ), f(g(gen(L))) )
    return ( to_field_elem.(arr), fg)
  else 
    to_field_elem = x -> preimage( f, x )
    return ( to_field_elem.(arr), f )
  end
end
  ╠═╡ =#

# ╔═╡ 1c6e0524-0b46-4a59-b4a9-8fd15936045d
# ╠═╡ disabled = true
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ e45d89f4-f033-4d9d-8930-02b58606873e
#=╠═╡
function cyclotomic_characters( r )
	mat, emb = to_combined_numberfield(characters(r),simplify_field=true)
	to_cyclotomic_field( mat, emb )
end
  ╠═╡ =#

# ╔═╡ bfbea5b1-cbc2-4bed-946e-2a2066c317e1
function samefield_characters( r )
	to_composite_field( characters( r ), simplify_field = true )
end

# ╔═╡ 1e211d49-5c37-4357-9f0c-b44b97d8bc2b
function export_characters( i::Int ) 
	fn1 = "/home/gert/Tests/characters/chars_"* string(i) *".mrdi"
	fn2 = "/home/gert/Tests/characters/chars_injection_"* string(i) *".mrdi"

	function export_new_chars(i) 
		try 
			chars, f  = samefield_characters(FRL[i])
			generator = gen( parent( chars[1] ) )
			
			Oscar.save( fn1, chars )
			Oscar.save( fn2, f(generator) )
		catch e2
			Oscar.save( fn1, ZZ.( [ 0 ] ) )
			Oscar.save( fn2, ZZ.( [ 0 ] ) )
		end	
	end
	
	try
		chars = Oscar.load(fn1)
		f     = Oscar.load(fn2)
		if chars == ZZ.([0])
			export_new_chars(i)
		end
	catch e
		export_new_chars(i)
	end
end

# ╔═╡ 881071c8-58dc-4f5f-99d3-7e031a3bfcf0
function create_ind()
	indices = []
	fn(i) = "/home/gert/Tests/characters/chars_injection_"* string(i) *".mrdi"
	for j in 1:352
		try 
			f = Oscar.load(fn(j))
			if f == ZZ.([0])
				push!( indices, j )
			end
			continue
		catch e
			push!( indices, j )
		end
	end
	indices
end

# ╔═╡ d656c78d-2df6-4497-82e8-0525a6195ea4


# ╔═╡ 76730eeb-53fa-4d07-b0be-cc14660e5011
@threads for i in create_ind()
	export_characters( i ) 
end

# ╔═╡ aeb7bcf2-fcc8-4e47-8457-0a8ddc59243c
# ╠═╡ disabled = true
#=╠═╡
function char_sort_crit( v )
	RR = ArbField(64);
	CC = AcbField(64);
	conv(x) = convert(Float64,x)
	# Abs values of elements of v
	absval(vec) = conv.( RR.(abs2.(vec)) )
	# Angles of elements of v
	angl(vec) = conv.( real.( log.( CC.( vec) ) ./ CC( 2 * pi * im ) ) )
			
	( Int( all(isreal.(v)) ), absval(v), angl(v) )
end
  ╠═╡ =#

# ╔═╡ 0c123f38-1d93-47da-bef7-35a45c444b6d
# ╠═╡ disabled = true
#=╠═╡
function characters(ring)
  if !(ring.characters === missing)
    return ring.characters
  elseif !FusionRings.is_commutative(ring) 
    error("Calculation of characters for non-commutative fusion ring is not implemented yet.")
  else
    qqb  = algebraic_closure(QQ) 
    mt   = FusionRings.multiplication_table( ring )
    r    = FusionRings.rank(ring)
    mats = [ matrix( qqb, mt[ i, :, : ] ) for i ∈ 1:r ]

    function is_diagonalizing_matrix( mat, mats )
      all( is_diagonal( mat * m * inv(mat) ) for m in mats )
    end
    
    diagq = false
    upi = 9
    upj = 9
    
    proposedchars = mats[1]
    while !diagq
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
      diagq = is_diagonalizing_matrix( proposedchars, mats )
    end

    normalize( mat ) = mat./mat[:,1]
    sort_mat( mat )  = sortslices( mat, dims = 1, by = char_sort_crit )

    sort_mat( normalize( [ proposedchars[i,j] for i in 1:r, j in 1:r ] ) )
  end 
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─95eeb20c-3efd-4285-b6a3-8d8d37aaee05
# ╠═4f8788dc-b431-11f0-01df-fd75003adb9b
# ╟─5b6a60f2-d925-49ed-a563-f8bc630f43c0
# ╠═cbaa67dc-35db-4759-aaab-09e461581c16
# ╠═ac3382fa-1341-4017-ae7a-b96bdfa0ea67
# ╟─e6ef9a18-94a6-44e5-b590-f08a12b2f4cf
# ╠═69d5a5bf-6952-4870-bbfb-39170b17940c
# ╟─2921df86-567e-4149-9677-c6fc821cb210
# ╠═93079240-ee00-4e8b-8bcd-ca98d1c7f016
# ╠═9e9d0d16-5907-4601-b59a-f03556541567
# ╟─3d817d5e-728a-45b3-bdba-7270a72a7c9b
# ╠═b25af0b1-010b-4a9d-8698-dab097f90ed7
# ╠═1c6e0524-0b46-4a59-b4a9-8fd15936045d
# ╠═e45d89f4-f033-4d9d-8930-02b58606873e
# ╠═bfbea5b1-cbc2-4bed-946e-2a2066c317e1
# ╠═1e211d49-5c37-4357-9f0c-b44b97d8bc2b
# ╠═881071c8-58dc-4f5f-99d3-7e031a3bfcf0
# ╠═d656c78d-2df6-4497-82e8-0525a6195ea4
# ╠═76730eeb-53fa-4d07-b0be-cc14660e5011
# ╠═aeb7bcf2-fcc8-4e47-8457-0a8ddc59243c
# ╠═0c123f38-1d93-47da-bef7-35a45c444b6d

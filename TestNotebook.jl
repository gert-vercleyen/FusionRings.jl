### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ 4f8788dc-b431-11f0-01df-fd75003adb9b
# ╠═╡ show_logs = false
begin
	using Revise
	using Pkg;
	#Pkg.develop(path="/home/gert/Projects/FusionRings.jl/")
	Pkg.develop(path="/Users/gertvercleyen/Projects/FusionRings.jl/")
	using FusionRings
	using Oscar
end

# ╔═╡ 95eeb20c-3efd-4285-b6a3-8d8d37aaee05
md"""
# Initialize package
"""

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

# ╔═╡ e1481c23-1f4e-4001-b457-337e42d684b3
r = psu2k_fusion_ring(7)

# ╔═╡ e63f950e-cf36-4199-a043-3238359d1be2


# ╔═╡ 835a0145-d7cf-41ce-9739-40bb0e05d89f
# ╠═╡ disabled = true
#=╠═╡
m = matrix( ZZ, mt(r)[2,:,:] )
  ╠═╡ =#

# ╔═╡ 324dc15c-256f-4b72-b351-ec5cd1c21e77
ch = characters(r)

# ╔═╡ b25af0b1-010b-4a9d-8698-dab097f90ed7
function to_combined_numberfield( 
    arr::Array{QQBarFieldElem}; 
    simplify_field = false, 
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

# ╔═╡ 42dd1286-4197-442d-93f6-0a66b8de2dc6

numch, emb = to_combined_numberfield(ch;simplify_field = true)

# ╔═╡ 6971f1d5-a26b-4a47-98a9-8bcf24ec8e12
emb( numch[1,1] ) == ch[1,1]

# ╔═╡ c136a86f-fe28-440d-b6a3-e2fa22340537
begin
	#qqb = algebraic_closure(QQ)
	arr = [ sqrt(qqb(2)), sqrt(qqb(13) + sqrt(qqb(13)) ) ]
	numarr, m = to_combined_numberfield( arr, simplify_field = true )
end

# ╔═╡ a46a09bd-48a4-4f4f-b94e-f6904b15e841
m(numarr[2])
>>>>>>> Stashed changes

# ╔═╡ 071269d5-9701-4fa4-b4fa-cb7c4dc2120c
emb( numch[1,1] )

# ╔═╡ b775919c-11bf-4d60-bcb1-14994c3513fb
ch[1,1]

# ╔═╡ 1c6e0524-0b46-4a59-b4a9-8fd15936045d
function to_cyclotomic_field( arr::Array{AbsSimpleNumFieldElem}, emb ) 
	length(arr) === 0 && return ( arr, emb )
	
	# Check parrent field of all fields are equal
	is_constant_array( parent.( arr ) ) || error("Elements of array should belong to same field")

	𝕂  = parent( arr[1] )
	C  = ray_class_field(simplify(𝕂)[1]) 
	𝕃, = C |> conductor |> first |> minimum |> Int |> cyclotomic_field
		
	i  = 
		hom( 
			𝕂, 
			𝕃, 
			first( roots( 𝕃, defining_polynomial(𝕂)/leading_coefficient(defining_polynomial(𝕂)) ) )
		)
	
	return ( i.(arr), ( emb, i ) )
end

# ╔═╡ 66779bbe-4aa9-4e5a-a99c-0b9c80f923b1
to_cyclotomic_field( numch, emb )

# ╔═╡ 14cf6184-8740-4e4f-8933-7e117f32d52e
methods(hom)

# ╔═╡ Cell order:
# ╟─95eeb20c-3efd-4285-b6a3-8d8d37aaee05
# ╠═4f8788dc-b431-11f0-01df-fd75003adb9b
# ╟─e6ef9a18-94a6-44e5-b590-f08a12b2f4cf
# ╠═69d5a5bf-6952-4870-bbfb-39170b17940c
# ╟─2921df86-567e-4149-9677-c6fc821cb210
# ╠═93079240-ee00-4e8b-8bcd-ca98d1c7f016
# ╠═9e9d0d16-5907-4601-b59a-f03556541567
# ╟─3d817d5e-728a-45b3-bdba-7270a72a7c9b
# ╠═e1481c23-1f4e-4001-b457-337e42d684b3
# ╠═e63f950e-cf36-4199-a043-3238359d1be2
# ╠═835a0145-d7cf-41ce-9739-40bb0e05d89f
# ╠═324dc15c-256f-4b72-b351-ec5cd1c21e77
# ╠═b25af0b1-010b-4a9d-8698-dab097f90ed7
# ╠═42dd1286-4197-442d-93f6-0a66b8de2dc6
# ╠═6971f1d5-a26b-4a47-98a9-8bcf24ec8e12
# ╠═a46a09bd-48a4-4f4f-b94e-f6904b15e841
# ╠═071269d5-9701-4fa4-b4fa-cb7c4dc2120c
# ╠═b775919c-11bf-4d60-bcb1-14994c3513fb
# ╠═1c6e0524-0b46-4a59-b4a9-8fd15936045d
# ╠═66779bbe-4aa9-4e5a-a99c-0b9c80f923b1
# ╠═14cf6184-8740-4e4f-8933-7e117f32d52e

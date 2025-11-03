### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 4f8788dc-b431-11f0-01df-fd75003adb9b
# ╠═╡ show_logs = false
begin
	using Pkg;
	Pkg.develop(path="/home/gert/Projects/FusionRings.jl/")
	using FusionRings
	using Oscar
end

# ╔═╡ 95eeb20c-3efd-4285-b6a3-8d8d37aaee05
md"""
# Initialize package
"""

# ╔═╡ e6ef9a18-94a6-44e5-b590-f08a12b2f4cf
md"""
# Initialize field
"""

# ╔═╡ 69d5a5bf-6952-4870-bbfb-39170b17940c
qqb = algebraic_closure(QQ)

# ╔═╡ 2921df86-567e-4149-9677-c6fc821cb210
md"""
# Some useful definitions 
"""

# ╔═╡ 93079240-ee00-4e8b-8bcd-ca98d1c7f016
mt = FusionRings.multiplication_table

# ╔═╡ e1481c23-1f4e-4001-b457-337e42d684b3
r = psu2k_fusion_ring(7)

# ╔═╡ 835a0145-d7cf-41ce-9739-40bb0e05d89f
# ╠═╡ disabled = true
#=╠═╡
m = matrix( ZZ, mt(r)[2,:,:] )
  ╠═╡ =#

# ╔═╡ 324dc15c-256f-4b72-b351-ec5cd1c21e77
v = [ 2 3 5 7 ]

# ╔═╡ d0613f9e-d892-4dd9-bdcb-a41f52c8397e
m2 = v[1] * matrix(qqb, mt(r)[2,:,:])

# ╔═╡ ae07fa0e-e5b6-4fc2-920b-f61e635cffc0
typeof(m2)

# ╔═╡ 3fa439a9-b9de-42a5-8cad-06a03cc0685c
generalized_jordan_form(m2)

# ╔═╡ 0454c6ff-a6a0-4582-8343-b72e971914ec
function eigenvectors( field, mat ) 
	zzmat = ZZMatrix( mat )
	evals = eigenvalues( field, mat )
end

# ╔═╡ 29897985-2b97-460e-86c0-9376db4ae8f9
function is_character_table( mat, ring::FusionRings.FusionRing )
	mt   = FusionRings.multiplication_table( ring )
	r    = FusionRings.rank(ring)
	mats = [ matrix( qqb, r, r, mt[ i, :, : ] ) for i ∈ 1:r ]
	all( is_diagonal.( mat * m * inv(mat) for m in mats ) )
end

# ╔═╡ 045e816d-ac5a-4b30-83c7-26424444708c
function is_character_table( mat, mats )
	all( is_diagonal( mat * m * inv(mat) ) for m in mats )
end

# ╔═╡ 0b893b9d-d238-4b2f-afdc-6bdec959c08e
function characters( ring )
	mt   = FusionRings.multiplication_table( ring )
	r    = FusionRings.rank(ring)
	mats = [ matrix( qqb, r, r, mt[ i, :, : ] ) for i ∈ 1:r ]

	function is_character_table( mat, mats )
		!(false ∈ map( is_diagonal, map( x -> mat * x * inv(mat), mats ) ))
	end
	
	charsq = false
	upi = 9
	upj = 9
	
	proposedchars = mats[1]
	while !charsq
		upi += 1
		upj += 1
		# Take random linear rational combination of fusion mats
		rvec 		= rand( [ i//j for i ∈ 1:upi, j ∈ 1:upj ], r )
		combinedmat = rvec[1] * mats[1]
		for i ∈ 2:r
			combinedmat += rvec[i] * mats[i]
		end

		# Find diagonalizing matrix
		proposedchars = generalized_jordan_form( combinedmat )[2]
		# Check whether procedure succeeded
		charsq = is_character_table( proposedchars, mats )
	end
	proposedchars
end

# ╔═╡ 1142ef18-2eb8-4846-895a-fb510d4fd546
x = characters(r)

# ╔═╡ de8ff568-b818-4d9e-92e0-3e27337d64ca
is_character_table( x, r )

# ╔═╡ Cell order:
# ╟─95eeb20c-3efd-4285-b6a3-8d8d37aaee05
# ╠═4f8788dc-b431-11f0-01df-fd75003adb9b
# ╟─e6ef9a18-94a6-44e5-b590-f08a12b2f4cf
# ╠═69d5a5bf-6952-4870-bbfb-39170b17940c
# ╟─2921df86-567e-4149-9677-c6fc821cb210
# ╠═93079240-ee00-4e8b-8bcd-ca98d1c7f016
# ╠═e1481c23-1f4e-4001-b457-337e42d684b3
# ╠═835a0145-d7cf-41ce-9739-40bb0e05d89f
# ╠═324dc15c-256f-4b72-b351-ec5cd1c21e77
# ╠═d0613f9e-d892-4dd9-bdcb-a41f52c8397e
# ╠═ae07fa0e-e5b6-4fc2-920b-f61e635cffc0
# ╠═3fa439a9-b9de-42a5-8cad-06a03cc0685c
# ╠═0454c6ff-a6a0-4582-8343-b72e971914ec
# ╠═29897985-2b97-460e-86c0-9376db4ae8f9
# ╠═045e816d-ac5a-4b30-83c7-26424444708c
# ╠═0b893b9d-d238-4b2f-afdc-6bdec959c08e
# ╠═1142ef18-2eb8-4846-895a-fb510d4fd546
# ╠═de8ff568-b818-4d9e-92e0-3e27337d64ca

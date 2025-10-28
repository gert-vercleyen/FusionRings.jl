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

# ╔═╡ d0613f9e-d892-4dd9-bdcb-a41f52c8397e
m2 = matrix(QQ, mt(r)[2,:,:])

# ╔═╡ ae07fa0e-e5b6-4fc2-920b-f61e635cffc0
typeof(m2)

# ╔═╡ 3fa439a9-b9de-42a5-8cad-06a03cc0685c
eigenspaces(m2)

# ╔═╡ 0454c6ff-a6a0-4582-8343-b72e971914ec
eigenvalues( qqb, ZZMatrix(FusionRings.multiplication_table(r)[2,:,:]) )

# ╔═╡ 29897985-2b97-460e-86c0-9376db4ae8f9


# ╔═╡ 0b893b9d-d238-4b2f-afdc-6bdec959c08e
function characters( ring )
	mt   = FusionRings.multiplication_table( ring )
	r    = FusionRings.rank(ring)
	mats = [ matrix( QQ, mt[ i, :, : ] ) for i ∈ 1:r ]
	rand = [ 2 3 5 7 ]
	combinedmat = rand[1] * mats[1,:,:]
	for i ∈ 2:r
		combinedmat += rand[i] * mats[i,:,:]
	end
	combinedmat
end

# ╔═╡ 1142ef18-2eb8-4846-895a-fb510d4fd546
x = characters(r)

# ╔═╡ 5b2935ff-aa98-4200-a5be-785b261a9aba
2*x[1,:,:]+3*x[2,:,:]

# ╔═╡ 29afb3e0-f5bd-4bb8-96e3-6f346b41708e


# ╔═╡ 835a0145-d7cf-41ce-9739-40bb0e05d89f
# ╠═╡ disabled = true
#=╠═╡
m = matrix( ZZ, mt(r)[2,:,:] )
  ╠═╡ =#

# ╔═╡ 371b6add-7709-4156-a362-4f1d834c060e
begin
	m = QQ[ 1 1 ; 0 1 ]
	print(jordan_decomposition(m))
	jordan_normal_form(m)
	generalized_jordan_form(m)
end

# ╔═╡ Cell order:
# ╟─95eeb20c-3efd-4285-b6a3-8d8d37aaee05
# ╠═4f8788dc-b431-11f0-01df-fd75003adb9b
# ╟─e6ef9a18-94a6-44e5-b590-f08a12b2f4cf
# ╠═69d5a5bf-6952-4870-bbfb-39170b17940c
# ╟─2921df86-567e-4149-9677-c6fc821cb210
# ╠═93079240-ee00-4e8b-8bcd-ca98d1c7f016
# ╠═e1481c23-1f4e-4001-b457-337e42d684b3
# ╠═835a0145-d7cf-41ce-9739-40bb0e05d89f
# ╠═d0613f9e-d892-4dd9-bdcb-a41f52c8397e
# ╠═ae07fa0e-e5b6-4fc2-920b-f61e635cffc0
# ╠═3fa439a9-b9de-42a5-8cad-06a03cc0685c
# ╠═0454c6ff-a6a0-4582-8343-b72e971914ec
# ╠═29897985-2b97-460e-86c0-9376db4ae8f9
# ╠═0b893b9d-d238-4b2f-afdc-6bdec959c08e
# ╠═1142ef18-2eb8-4846-895a-fb510d4fd546
# ╠═5b2935ff-aa98-4200-a5be-785b261a9aba
# ╠═371b6add-7709-4156-a362-4f1d834c060e
# ╠═29afb3e0-f5bd-4bb8-96e3-6f346b41708e

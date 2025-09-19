
# FusionRings.jl

A lightweight toolbox for constructing and working with fusion rings, with optional OSCAR/GAP extensions and a simple JSON registry.

## Quick start

```julia
using Pkg
Pkg.activate(; temp=true)
Pkg.develop(path="C:/path/to/FusionRingsPackage")
Pkg.instantiate()

using FusionRings

F = fibonacci_ring()
print_tensor_table(F)
tensor_product(F, :τ, :τ)    # Dict(:"1"=>1, :τ=>1)
decompose(F, :τ, :τ)         # [(:"1",1), (:τ,1)]
product_string(F, :τ, :τ)    # "τ ⊗ τ = 1 ⊕ τ"
```

## SONG (D₃)

```julia
R = song_extension_D3(:rotations; tildeg="s", n=1, A=:id)
print_tensor_table(R)
exts = enumerate_song_extensions_D3()
```

## JSON registry

```julia
id = save_ring_json(F; extra=Dict("source"=>"demo"))
hits = query_registry_by_fpdims(quantum_dimensions(F))
```


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
tensor_product(F, "τ", "τ")   # Dict("1"=>1, "τ"=>1)
tensor_product(F, :τ, :τ)      # also works (symbols auto-convert)
decompose(F, "τ", "τ")        # [("1",1), ("τ",1)]
product_string(F, "τ", "τ")   # "τ ⊗ τ = 1 ⊕ τ"
```

### Note on label types

As of version NEXT, fusion ring simple object labels are stored internally as `String` instead of `Symbol`.

Rationale:
* Strings are more natural for user-provided names (UTF-8, no interning concerns).
* Symbol semantics (identity & interning) were not exploited in arithmetic; labels serve as display identifiers.
* Avoids accidental global symbol growth in long-running sessions.

Backward compatibility:
* All public APIs still accept `Symbol` inputs; they are converted with `String(sym)`.
* Returned collections (e.g. `labels(fr)`, keys of `tensor_product`) now contain `String` values.
* Update downstream code that relied on `Symbol` pattern matching to use strings (or call `Symbol(x)` if necessary).

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

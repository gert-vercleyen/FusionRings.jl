
# FusionRings.jl

A lightweight toolbox for constructing and working with fusion rings, with optional OSCAR/GAP extensions and a simple JSON registry.

## Quick start

```julia
using Pkg
Pkg.activate(; temp=true)
Pkg.develop(path="C:/path/to/FusionRings.jl")  # adjust path
Pkg.instantiate()

using FusionRings

F = fibonacci_ring()
print_multiplication_table(F)
tensor_product(F, "τ", "τ")    # Dict("1"=>1, "τ"=>1)
tensor_product(F, :τ, :τ)       # also works (Symbols auto-convert)
decompose(F, "τ", "τ")         # [("1",1), ("τ",1)]
product_string(F, "τ", "τ")    # "τ ⊗ τ = 1 ⊕ τ"
fpdims(F)                       # Frobenius–Perron dims vector
fpdim(F)                        # global FP dimension
```

### Conjugation API (since 0.5.0 pending)

`conjugate_element(fr, a)` now returns the integer index of the dual object.
To obtain the display label, call `conjugate_label(fr, a)`.

Motivation: keeping computations (composition, iteration, comparisons) in
index space avoids repeated label lookups and allocations; labels remain for
printing and external serialization only.

### Subring checks

`is_sub_fusion_ring(fr, subset::Vector)` checks closure of a chosen label
subset (must contain the unit). New: `is_sub_fusion_ring(big, small)` checks
that the fusion data of `small` is literally the restriction of `big` to the
matching labels (no relabeling search is attempted).

### Named rings module rename (internal)

The internal source file defining predefined example rings was renamed from
`FusionRingGenerators.jl` to `NamedFusionRings.jl`. A deprecated compatibility
stub remains so user code continues to work unchanged.

### Note on label types (since 0.3.0)

As of version 0.3.0, fusion ring simple object labels are stored internally as `String` instead of `Symbol`.

Rationale:
* Strings are more natural for user-provided names (UTF-8, no interning concerns).
* Symbol semantics (identity & interning) were not exploited in arithmetic; labels serve as display identifiers.
* Avoids accidental global symbol growth in long-running sessions.

Backward compatibility:
* All public APIs still accept `Symbol` inputs; they are converted with `String(sym)`.
* Returned collections (e.g. `labels(fr)`, keys of `tensor_product`) now contain `String` values.
* Update downstream code that relied on `Symbol` pattern matching to use strings (or call `Symbol(x)` if necessary).

### API renames in 0.4.0

Version 0.4.0 introduces clearer, domain-standard names:

* `quantum_dimensions` → `fpdims`
* `global_dimension` → `fpdim`
* `tensor_table` → `multiplication_table`
* `print_tensor_table` → `print_multiplication_table`
* `fusion_ring_from_group_oscar` → `group_fusion_ring`
* `near_group_ring` → `near_group_fusion_ring`

The old names are deprecated but still callable for now; please migrate to the new names.

## SONG (D₃)

```julia
R = song_extension_D3(:rotations; tildeg="s", n=1, A=:id)
print_multiplication_table(R)
exts = enumerate_song_extensions_D3()
```

## JSON registry

```julia
id = save_ring_json(F; extra=Dict("source"=>"demo"))
hits = query_registry_by_fpdims(fpdims(F))
```


## Quick examples
```julia
using FusionRings
F = fibonacci_ring()
print_multiplication_table(F)
save_ring_json(F)
```

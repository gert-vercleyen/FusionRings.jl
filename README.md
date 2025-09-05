# FusionRings (clean, dependency-light)

This package provides a validated `FusionRing` type and a small suite of generators
(Fibonacci, Ising, Semion, Zₙ, SU(2)_k, PSU(2)_k, Tambara–Yamagami, near-group)
plus utilities (quantum dimensions, sub-rings, equivalence).

## Install (Windows-friendly)

```julia
using Pkg
Pkg.develop(path="C:/path/to/FusionRingsPackage")
using FusionRings
```

If you're inside the folder, simply: `Pkg.develop(path=".")`.

## Quickstart

```julia
using FusionRings

F = fibonacci_ring()
tensor_product(F, :τ, :τ)            # Dict(:1=>1, :τ=>1)

G = zn_fusion_ring(5)                 # Z5
is_commutative(G)                     # true

H = su2_fusion_ring(4)
quantum_dimensions(H)

TY = ty_fusion_ring(["0","1","2"])
NG = near_group_ring(["0","1"]; n=1)  # near-group on Z2 with n=1
```

## Tests
```
] test FusionRings
```

## Notes
- No external registry packages are required.
- SONG-specific enumeration is not included in this light build; a general
  near-group constructor is provided.

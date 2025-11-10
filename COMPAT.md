
# Compatibility map (Anyonica ↔ FusionRings.jl)

This package exports a thin shim `src/CompatAnyonica.jl` with names used in the Anyonica codebase.
Where functionality exists in this package, wrappers delegate to the implemented methods.
Otherwise, wrappers throw a clear `Not implemented yet` error so downstream code fails loudly.

**Aliases implemented now**
- `adjoint_irreps` → `GroupTheory.adjoint_irreps`
- `universal_grading` → `GroupTheory.universal_grading`
- `fusion_ring_automorphisms` → `GroupTheory.fusion_ring_automorphisms`
- `num_self_dual`/`nsd` and `num_non_self_dual`/`nnsd` → `Properties`

**Placeholders (PRs welcome)**
- `adjoint_fusion_ring`, `which_decompositions`, `is_nilpotent_fusion_ring`, `commutator`, `upper_central_series`,
  `identify_fusion_ring`, `group_fusion_ring`, `son2_fusion_ring`, `hi_fusion_ring`

The main codebase keeps **labels as strings** to avoid `Symbol` conversion issues and to match typical usage in fusion ring computations.

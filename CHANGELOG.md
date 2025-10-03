# Changelog

## 0.5.0 - 2025-10-02
### Breaking
- `conjugate_element` now returns an `Int` (index) instead of a `String` label.

### Added
- `conjugate_label(fr,a)` returns the label of the dual object.
- `is_sub_fusion_ring(big::FusionRing, small::FusionRing)` overload for structural subring checks.

### Changed
- Subring logic clarified: vector form (`is_sub_fusion_ring(fr, [labels...])`) unchanged; new ring-based form matches labels directly (no permutation search).

### Internal
- Renamed source file `FusionRingGenerators.jl` → `NamedFusionRings.jl` (stub module keeps old name with a warning).

### Notes
- Minor version bump due to return type break in `conjugate_element`.

## 0.4.0 - 2025-10-02
### Changed
- Renamed `quantum_dimensions` to `fpdims` (Frobenius–Perron dimensions).
- Renamed `global_dimension` to `fpdim` (global Frobenius–Perron dimension).
- Renamed `tensor_table` to `multiplication_table`.
- Renamed `print_tensor_table` to `print_multiplication_table`.
- Renamed `fusion_ring_from_group_oscar` to `group_fusion_ring`.
- Renamed `near_group_ring` to `near_group_fusion_ring` for clarity/consistency.

### Deprecated
- Old names listed above now issue deprecation warnings but still function; they will be removed in a future minor release (≥ 0.6.0).

### Added
- Documentation updates reflecting new API names.

### Internal
- Consolidated registry queries and docs to use `fpdims` / `fpdim` naming.

### Notes
- This is a minor version bump due to public API renames (non-breaking via deprecations). Update downstream code to silence warnings.

## 0.3.0 - 2025-09-26
### Changed
- Element labels internally and in all public-facing APIs are now `String` instead of `Symbol`.
- All generators (`fibonacci_ring`, `ising_ring`, etc.) now return rings whose `labels` field is `Vector{String}`.
- Operations (`tensor_product`, `fusion_matrix`, `decompose`, etc.) return dictionaries keyed by `String`.

### Added
- Backward-compatible input handling: functions still accept `Symbol` or `String` (and integer indices where previously supported); inputs are normalized to `String`.
- Added test target configuration in `Project.toml`.

### Fixed
- Removed parse error in `SongsD3.parse_d3` by replacing chained ternary with explicit conditionals.
- Eliminated `LinearAlgebra.rank` name conflict warnings by narrowing imports.
- Ensured top-level module re-exports generator and operation symbols so `using FusionRings` exposes the full API.

### Internal
- Added shared `indexmap(fr)` helper with simple caching to reduce repeated Dict constructions.
- Simplified normalization logic removing redundant Symbol vs String branches.
- Removed uniqueness assertion in `conjugate_element` (trust constructor invariants).
- Normalized extension (OSCAR) label generation to produce `String` labels directly.

### Notes
- This is a minor version bump because although inputs remain backward compatible, return types involving labels changed from `Symbol` to `String`.
- If external code relied on `keys(dict)::Symbol` it should be updated to compare with strings or call `String(sym)`.

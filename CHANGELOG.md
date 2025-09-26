# Changelog

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

### Notes
- This is a minor version bump because although inputs remain backward compatible, return types involving labels changed from `Symbol` to `String`.
- If external code relied on `keys(dict)::Symbol` it should be updated to compare with strings or call `String(sym)`.

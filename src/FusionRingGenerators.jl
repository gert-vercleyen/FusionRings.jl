module FusionRingGenerators

"""Deprecated module. Use `NamedFusionRings` instead. This stub re-exports all
named ring constructors and will be removed in a future release."""
const _deprecated_warned = Ref(false)
function __init__()
    if !_deprecated_warned[]
        @warn "Module FusionRingGenerators is deprecated; use NamedFusionRings instead" maxlog=1
        _deprecated_warned[] = true
    end
end
using ..NamedFusionRings
export fibonacci_ring, ising_ring, semion_ring, su2_fusion_ring, psu2_fusion_ring,
       zn_fusion_ring, ty_fusion_ring, near_group_fusion_ring

end

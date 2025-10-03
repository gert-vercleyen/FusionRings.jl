module GeneralFunctions

using ..Types: FusionRing, labels

"""
    indexmap(fr::FusionRing) -> Dict{String,Int}

Return (and cache per-session) a dictionary mapping each simple object label
to its index. This centralizes repeated constructions that previously
occurred inline in multiple operations.

The mapping is inexpensive to build for small ranks, but many functions call
it repeatedly; having a single helper makes future memoization trivial if
benchmarks suggest it matters.
"""
module_indexmaps = IdDict{FusionRing,Dict{String,Int}}()
function indexmap(fr::FusionRing)
    get!(module_indexmaps, fr) do
        Dict(l=>i for (i,l) in enumerate(labels(fr)))
    end
end

export indexmap

end # module GeneralFunctions


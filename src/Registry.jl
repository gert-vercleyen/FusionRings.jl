
module Registry

using ..Types: FusionRing, fusion_tensor, labels, rank, ringname
using ..Properties: quantum_dimensions, global_dimension, multiplicity, is_commutative
using JSON3, Dates, UUIDs

export save_ring_json, load_registry_json, query_registry_by_fpdims

const DEFAULT_REGISTRY_PATH = normpath(joinpath(@__DIR__, "..", "data-cache", "registry.json"))

function save_ring_json(fr::FusionRing; path::AbstractString=DEFAULT_REGISTRY_PATH,
                        extra::Dict=Dict{String,Any}(), overwrite::Bool=false)
    mkpath(dirname(path))
    rec = Dict{String,Any}(
        "id"            => string(uuid4()),
        "timestamp"     => Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sZ"),
        "name"          => ringname(fr),
        "rank"          => rank(fr),
        "labels"        => string.(labels(fr)),
        "fpdims"        => quantum_dimensions(fr),
        "fpdims_sorted" => sort(quantum_dimensions(fr)),
        "global_dim"    => global_dimension(fr),
        "multiplicity"  => multiplicity(fr),
        "commutative"   => is_commutative(fr),
        "tensor"        => fusion_tensor(fr),
    )
    for (k,v) in extra
        rec[k] = v
    end
    data = if isfile(path)
        JSON3.read(read(path, String))
    else
        JSON3.Array()
    end
    push!(data, rec)
    open(path, "w") do io
        JSON3.write(io, data; indent=2)
    end
    rec["id"]
end

function load_registry_json(; path::AbstractString=DEFAULT_REGISTRY_PATH)
    isfile(path) || return JSON3.Array()
    JSON3.read(read(path, String))
end

function query_registry_by_fpdims(dims::AbstractVector{<:Real};
                                  path::AbstractString=DEFAULT_REGISTRY_PATH,
                                  tol::Real=1e-8)
    want = sort(collect(dims))
    reg = load_registry_json(; path)
    hits = Any[]
    for rec in reg
        v = get(rec, "fpdims_sorted", nothing)
        v===nothing && continue
        length(v) == length(want) || continue
        ok = all(abs(v[i] - want[i]) ≤ tol for i in eachindex(v))
        ok && push!(hits, rec)
    end
    hits
end

end

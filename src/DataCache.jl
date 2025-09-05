module DataCache
using ..Types, ..Creation
using Serialization
export optimized_import, FusionRingByCode, FRBC, FusionRingList, FRL, all_fusion_rings_q, AFRQ

const _data_dir = joinpath(@__DIR__, "..", "data")
isdir(_data_dir) || mkpath(_data_dir)

function optimized_import(base::AbstractString)
    path = joinpath(_data_dir, base * ".bson")
    isfile(path) || error("Missing data file: " * path * " (provide your cached dataset).")
    return Serialization.deserialize(path)
end

const FusionRingList = isfile(joinpath(_data_dir,"FusionRingList.bson")) ? optimized_import("FusionRingList") : Types.FusionRing[]
const FRL = FusionRingList

function FusionRingByCode(t::NTuple{N,Int}) where N
    dictpath = joinpath(_data_dir, "FusionRingAssociation.bson")
    isfile(dictpath) || error("Missing data file: " * dictpath)
    db = Serialization.deserialize(dictpath)  # expect Dict{NTuple{4,Int}, Any}
    if N==3
        return get(db, (t[1],1,t[2],t[3])) do
            error("Code $(t) not found")
        end
    elseif N==4
        return get(db, t) do
            error("Code $(t) not found")
        end
    else
        error("Code must be a 3- or 4-tuple")
    end
end
const FRBC = FusionRingByCode

function all_fusion_rings_q(r::Int, m::Int)
    r == 1 ||
    (r in (2,3,4) && m <= 16) ||
    (r == 5 && m <= 12) ||
    (r == 6 && m <= 4) ||
    (r == 7 && m <= 2) ||
    (r == 8 && m == 1) ||
    (r == 9 && m == 1)
end
const AFRQ = all_fusion_rings_q

end # module


module DataCache

using ..Types: FusionRing
using ..FusionRingGenerators: fibonacci_ring, ising_ring
using Serialization

export optimized_import, FusionRingList, FRL, FusionRingByCode, FRBC, AllFusionRingsQ, AFRQ

const DATADIR = normpath(joinpath(@__DIR__, "..", "data-cache"))

function optimized_import(basename::String)
    mkpath(DATADIR)
    jls = joinpath(DATADIR, basename*".jls")
    if isfile(jls)
        return deserialize(jls)
    end
    if isdefined(@__MODULE__, :_bson_read)
        return _bson_read(joinpath(DATADIR, basename*".bson"))
    end
    nothing
end

const FusionRingList = let data = optimized_import("FusionRingList")
    data === nothing ? [fibonacci_ring(), ising_ring()] : data
end
const FRL = FusionRingList

const _FRBC_ASSOC = let data = optimized_import("FusionRingAssociation")
    data === nothing ? Dict((2,1,1,1)=>fibonacci_ring()) : data
end

function FusionRingByCode(code::NTuple{4,Int})
    get(_FRBC_ASSOC, code, error("Code not found in cache"))
end
function FusionRingByCode(code::NTuple{3,Int})
    FusionRingByCode((code[1], 1, code[2], code[3]))
end
const FRBC = FusionRingByCode

function AllFusionRingsQ(r::Int, m::Int)
    if r==1; return true; end
    if r in (2,3,4) && m≤16; return true; end
    if r==5 && m≤12; return true; end
    if r==6 && m≤4; return true; end
    if r==7 && m≤2; return true; end
    if r==8 && m==1; return true; end
    if r==9 && m==1; return true; end
    false
end
const AFRQ = AllFusionRingsQ

end

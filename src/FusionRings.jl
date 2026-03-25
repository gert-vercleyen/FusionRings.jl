module FusionRings

using Oscar, Combinatorics, JSON, Base.Threads, Accessors
using LinearAlgebra: eigen, eigvals, diag
import Oscar: multiplication_table, is_commutative, rank, multiplicity
import Base.names

include("GeneralFunctions.jl")
include("Structs.jl")
include("FormattingAndPrinting.jl")
include("Creation.jl")
include("Properties.jl")
include("Operations.jl")
include("ImportData.jl")

export qqb_dict, fusion_ring_list, frl, fusion_ring_dict, frd, from_anyonwiki_code, fawc

"""Internal helper: ensure ring lookup dictionary has been initialized."""
function _ensure_frd_initialized()
    isdefined(@__MODULE__, :frd) || error("FusionRings data not initialized (frd is undefined).")
    frd isa Dict || error("frd is defined but not a Dict).")
    nothing
end

"""
    from_anyonwiki_code(r, m, nnsd, i) -> FusionRing
    from_anyonwiki_code(v::AbstractVector{<:Integer}) -> FusionRing

Lookup a stored fusion ring by its AnyonWiki code `[r, m, nnsd, i]`.
This requires package data to be loaded (normally happens during `__init__`).
"""
function from_anyonwiki_code(r::Integer, m::Integer, nnsd::Integer, i::Integer)
    _ensure_frd_initialized()
    frd[Int[r, m, nnsd, i]]
end

function from_anyonwiki_code(v::AbstractVector{<:Integer})
    _ensure_frd_initialized()
    length(v) == 4 || error("anyonwiki_code expects a vector of 4 integers.")
    frd[Int.(collect(v))]
end

const fawc = from_anyonwiki_code
const awc = from_anyonwiki_code

function __init__()
    # GLOBAL VARIABLES
    global QQb     = algebraic_closure(QQ)
    global QQab, ζ = abelian_closure(QQ)

    datadir = joinpath( @__DIR__, "data" )

    # IMPORT DICTIONARY OF QQB ELEMENTS
    global qqb_dict = begin
        ids     = Oscar.load( joinpath( datadir, "qqb_ids.mrdi") )
        nums    = Oscar.load( joinpath( datadir, "qqb_vals.mrdi") )

        Dict( ids[i] => nums[i] for i in 1:length(ids) )
    end

    # IMPORT FUSION RINGS
    global fusion_ring_list =
        sort( # Stored list is unsorted so we still need to sort
            import_rings( joinpath( datadir, "fusionrings.json" ) ),
            by = ( x -> (x.anyonwiki_code)[ [ 2, 1, 3, 4 ] ] )
        )
    global frl = fusion_ring_list

    global fusion_ring_dict = Dict( anyonwiki_code(r) => r for r in frl )
    global frd = fusion_ring_dict
end

end

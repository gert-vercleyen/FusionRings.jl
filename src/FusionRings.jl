module FusionRings

using Oscar, Combinatorics, JSON, LinearAlgebra:eigen, Base.Threads, Accessors
import Oscar: multiplication_table, is_commutative, rank, multiplicity
import Base.names

include("GeneralFunctions.jl")
include("Creation.jl")
include("Properties.jl")
include("Operations.jl")
include("ImportData.jl")
include("FormattingAndPrinting.jl")

export qqb_dict, fusion_ring_list, frl, fusion_ring_dict, frd, from_anyonwiki_code, fawc

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

    function from_anyonwiki_code( r, m, nnsd, i )
         frd[ [ r, m, nnsd, i ] ]
    end
    
    function from_anyonwiki_code( v::Vector{Int64} )
        if length( v ) == 4
            return frd[ v ]
        else
            error( "anyonwiki_code expects a vector of 4 Int64's." )
        end
    end
    global fawc = from_anyonwiki_code
end

end

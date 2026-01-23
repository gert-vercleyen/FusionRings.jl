module FusionRings

using Oscar
import Oscar: multiplication_table, is_commutative, rank, multiplicity
using Combinatorics
using JSON
using LinearAlgebra:eigen
using Base.Threads
using Accessors

include("GeneralFunctions.jl")
include("Creation.jl")
include("Properties.jl")
include("Operations.jl")
include("ImportData.jl")
include("FormattingAndPrinting.jl")

export qqb_dict
export fusion_ring_list
export frl
export fusion_ring_dict
export frd
export from_anyonwiki_code
export fawc
function __init__()
    global QQb     = algebraic_closure(QQ)
    global QQab, ζ = abelian_closure(QQ)

    
    global qqb_dict = begin
        datadir = joinpath( @__DIR__, "data", "Numbers", "QQBarFieldElems" )
        ids     = Oscar.load( joinpath( datadir, "qqb_ids.mrdi") )
        nums    = Oscar.load( joinpath( datadir, "qqb_vals.mrdi") )

        Dict( ids[i] => nums[i] for i in 1:length(ids) )
    end


    local function load_frl()
        path = joinpath( @__DIR__, "data", "FusionRingsJSON" )

        prep_path( fn ) = joinpath( path, fn )

        filenames = prep_path.( readdir(path) )

        [ import_ring( fn ) for fn in filenames ]
    end

    # We want the mult-free rings to be first
    function frlsortcrit( c )
        (c.anyonwiki_code)[ [ 2, 1, 3, 4 ] ]
    end

    global fusion_ring_list = sort( load_frl(), by = frlsortcrit )
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

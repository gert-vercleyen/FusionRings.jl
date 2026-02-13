using Pkg; Pkg.activate("/home/gert/Projects/FusionRings.jl")
using Revise, Oscar, Accessors, JSON, FusionRings 

to_qqb_arr( mat ) = [ mat[i] for i in eachindex(mat) ] 

function adddata( fr )
    code = string.( anyonwiki_code( fr ) )
    charfn =
        "chars" *
        code[1] * "_" *
        code[2] * "_" *
        code[3] * "_" *
        code[4] * ".mrdi"
    chrs =
        try 
            qqb_id(
                to_qqb_arr(
                    Oscar.load(
                        joinpath(
                            "/home/gert/Projects/FusionRings.jl/src/data/Characters/",
                            charfn
                        )
                    )
                )
            )
        catch e
            missing
        end

    ncr = is_categorifiable(fr) === true ? "none" : nothing
   
    newfr = @set fr.characters = chrs 
    if chrs === missing
        newfr = @set fr.numeric_characters = chrs
    end
    @set newfr.non_cat_reason = ncr
end

function new_qqb_id( x::QQBarFieldElem ) 
    mp = minimal_polynomial(x)
    coeffs = string.( collect( coefficients(mp) ) )
    us = fill( "_", degree(mp) + 1 )
    
    numstring = string( rootnum( x ) )

    stringriffle( coeffs, us ) *  "_" * numstring
end

function fr_fn(fr::FusionRing)
    c  = string.( anyonwiki_code( fr ) )
    us = fill( "_", 3 )
    stringriffle( c, us ) * ".json"
end

function export_new_rings(i::Int64)
    datadir = "/home/gert/Projects/FusionRings.jl/src/data/NewRings/"
    probrings = []
    for fr in frl[1:i]
        try
            export_ring(
                joinpath( datadir, fr_fn(fr) ),
                adddata(fr)
            )
        catch error
            push!(probrings, anyonwiki_code(fr) )
            continue
        end
    end
    probrings
end

function export_new_rings(v::Vector{Int64})
    datadir = "/home/gert/Projects/FusionRings.jl/src/data/NewRings/"
    probrings = []
    for fr in frl[v]
        #try
            export_ring(
                joinpath( datadir, fr_fn(fr) ),
                adddata(fr)
            )
        #catch error
        #    print(error)
        #    push!(probrings, anyonwiki_code(fr) )
        #    continue
        #end
    end
    probrings
end
export_new_rings(20)

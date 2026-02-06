using Pkg; Pkg.activate("/home/gert/Projects/FusionRings.jl")
using Revise, Oscar, Accessors, JSON, FusionRings 

to_qqb_arr( mat ) = [ mat[i] for i in eachindex(mat) ] 

function import_and_fix( fr )
    path = "/home/gert/Projects/FusionRings.jl/src/data/FusionRings/"
    
    js = JSON.parsefile( joinpath( path, fr_fn( fr ) ) )

    # Mult tables were incorrect
    mt = multiplication_table( fr )

    # Song names contain \\| which should be |
    tnames = replace( tnfromjs( js ), "\\|" => "|" )
    
    iscategorifiable = cfromjs( js ) 
    if iscategorifiable === true # could be missing
         
        prpdict = ctpfromjs( js )
        uq = prpdict["Unitary"]
        sq = prpdict["Spherical"]
        bq = prpdict["Braided"]
        rq = prpdict["Ribbon"]
        mq = prpdict["Modular"]

        function reasontext( string, bool )
            bool === missing && [  ]
            
            if bool
                [ "Example", "See categorifications" ]
            else
                [ "Computer", "No "* string *" category found using Anyonica v < 0.9.7.0"]
            end
        end
        
        catswithprops =
            [
                 [ "Fusion"    , true, reasontext( "fusion", true ) ]
                ,[ "Pivotal"   , true, reasontext( "pivotal", true ) ]
                ,[ "Unitary"   , uq, reasontext( "unitary", uq ) ]
                ,[ "Spherical" , sq, reasontext( "spherical", sq ) ]
                ,[ "BraidedQ"  , bq, reasontext( "braided", bq )]
                ,[ "Ribbon"    , rq, reasontext( "ribbon", rq ) ]
                ,[ "Modular"   , mq, reasontext( "modular", mq ) ]
            ]
    elseif iscategorifiable === false
        catswithprops =
            [
                 [ "Fusion"    , false, [ "Computer", "No category found using Anyonica v < 0.9.7.0" ] ]
                ,[ "Pivotal"   , false, [ "Computer", "See fusion category" ] ]
                ,[ "Unitary"   , false, [ "Computer", "See fusion category" ] ]
                ,[ "Spherical" , false, [ "Computer", "See fusion category" ] ]
                ,[ "Braided"   , false, [ "Computer", "See fusion category" ] ]
                ,[ "Ribbon"    , false, [ "Computer", "See fusion category" ] ]
                ,[ "Modular"   , false, [ "Computer", "See fusion category" ] ]
            ]
    else 
        catswithprops =
            [
                 [ "Fusion"    , missing, [ "", "Unknown to AnyonWiki" ] ]
                ,[ "Pivotal"   , missing, [ "", "Unknown to AnyonWiki" ] ]
                ,[ "Unitary"   , missing, [ "", "Unknown to AnyonWiki" ] ]
                ,[ "Spherical" , missing, [ "", "Unknown to AnyonWiki" ] ]
                ,[ "Braided"   , missing, [ "", "Unknown to AnyonWiki" ] ]
                ,[ "Ribbon"    , missing, [ "", "Unknown to AnyonWiki" ] ]
                ,[ "Modular"   , missing, [ "", "Unknown to AnyonWiki" ] ]
            ]
    end
    
    fusion_ring(
        mt,
        names                               = nfromjs( js ),
        texnames                            = tnames,
        barcode                             = bcfromjs( js ),
        anyonwiki_code                      = fcfromjs( js ),
        characters                          = chfromjs( js ),
        sub_fusion_rings                    = sfrfromjs( js ),
        projective_SL2Z_reps                = psrfromjs( js ),
        frobenius_perron_dimension          = fpdfromjs( js ),
        frobenius_perron_dimensions         = fpdsfromjs( js ),
        tensor_product_decompositions       = tpdfromjs( js ),
        numeric_characters                  = nchfromjs( js ),
        numeric_projective_SL2Z_reps        = npsrfromjs( js ),
        numeric_frobenius_perron_dimension  = nfpdfromjs( js ),
        numeric_frobenius_perron_dimensions = nfpdsfromjs( js ),
        has_categories_with_props           = catswithprops,
        categorifications                   = ctsfromjs( js ),
        references                          = js["references"],
        software                            = js["software"],
        comments                            = js["comments"],
  )
end

function fix_characters(fr::FusionRing)
    ch =  fr.characters
    ch === missing && return fr
    
    @set fr.characters = string.(permutedims(ch) )
end

function fr_fn(fr::FusionRing)
    c  = string.( anyonwiki_code( fr ) )
    us = fill( "_", 3 )
    stringriffle( c, us ) * ".json"
end

function export_new_rings(v::Vector{Int64})
    datadir = "/home/gert/Projects/FusionRings.jl/src/data/FusionRings/"
    export_rings( joinpath(datadir,"newrings.json"), import_and_fix.(frl[v]) )
end


macro timeout(seconds, expr_to_run, expr_when_fails=nothing)
    quote
        timer = Channel{Timer}(1)
        tsk = @task begin
            x = $(esc(expr_to_run))
            close(take!(timer))
            x
        end
        schedule(tsk)

        put!(timer, Timer($(esc(seconds))) do timer
            istaskdone(tsk) || schedule(tsk, InterruptException(); error=true)
        end)

        try
            fetch(tsk)
        catch _
            $(esc(expr_when_fails))
        end
    end
end

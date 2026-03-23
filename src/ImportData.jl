############################################################
# Exporting and importing QQBarFieldElems
############################################################

# The QQBarFieldElem objects are quite heavy to load so we will load them all in a dictionary 
# and provide functions to convert qqbar elems to keys and vice versa

# Generate unique ID for a QQBarFieldElem
export qqb_id

function qqb_id( x::QQBarFieldElem ) 
    mp = minimal_polynomial(x)
    coeffs = string.( collect( coefficients(mp) ) )
    us = fill( "_", degree(mp) + 1 )
    
    numstring = string( rootnum( x ) )

    stringriffle( coeffs, us ) *  "_" * numstring
end

qqb_id( arr::Array{QQBarFieldElem} ) = qqb_id.(arr)

function rootnum( x::QQBarFieldElem )
    p   = minimal_polynomial( x ) 
    rts = roots( QQBar, p )
    sr  = sort( rts, by = root_sort_crit )
    findfirst( y -> y == x, sr )
end

# This is the sort criterion for roots used by mathematica 
# and by the anyonwiki on 28/12/2025
function root_sort_crit( x )
    ( - Int( is_real( x ) ), real(x), imag(x) )
end


#TODO:implement
function save_qqb_num( dir::String, x::QQBarFieldElem )
    return nothing
end

#TODO:implement
function save_qqb_num( x::QQBarFieldElem )
    return nothing
end

# Load the dictionary of qqbar elems

# Get from dict
export from_qqb_id

from_qqb_id( s::String ) = qqb_dict[s]
from_qqb_id( a::Array{Any} ) = from_qqb_id.(a)
from_qqb_id( a::Matrix{String} ) = from_qqb_id.(a)

############################################################
# Importing fusion rings
############################################################
# The fusion rings are stored as json files. Not all data 
# types (e.g. complex numbers) are supported by JSON so we 
# had to store those using a variety of hacks. 
# The following functions convert the stored data back 
# to their proper types.
#
# TODO: some of the if clauses below are necessary for legacy
# compatibility. Once all json files have the correct format
# we should remove it since it slows down the import

# formal code
function fcfromjs( js )::Vector{Int64}
    k = keys( js )

    if "formal_code" ∈ k
        fc = js["formal_code"]
    elseif "anyonwiki_code" ∈ k
        fc = js["anyonwiki_code"]
    else
        return missing
    end

    if length(fc) == 0
        missing
    else
        [ fc[i] for i in 1:4 ]  
    end
end

# mult tab
function mtfromjs( js )::Array{Int64, 3}
  jsmt = js["mult_tab"]
  r = length(jsmt)
  mt = zeros(Int, r, r, r)
  for i in 1:r, j in 1:r, k in 1:r 
      mt[i,j,k] = Int.(jsmt[i][j][k])
  end
  mt
end

# barcode
function bcfromjs(js)::ZZRingElem
  ZZ( parse( BigInt, js["barcode"] ) )
end

# tensor product decompositions
function tpdfromjs(js)
  tps = js["tensor_product_decompositions"]

    if length(tps) == 0
        []
    elseif typeof(tps) === Vector{Any}
        [
            [ Int.( code ) for code in decomp ]
            for decomp in tps
        ]
    else
        tps = tps["value"]
        [
            [ Int.( code ) for code in decomp ]
            for decomp in tps
        ]
    end
end

# sub-fusion rings
function sfrfromjs(js)
    srs = js["non_trivial_sub_fusion_rings"]

    # if length(srs) == 0
    #     return []
    # end

    # intData =
    #     [
    #         [ Int.( sr["injection"] ), Int.( sr["anyonwiki_code"] ) ]
    #         for sr in srs
    #     ]

    # [
    #     Dict(
    #         "injection"       => data[1],
    #         "anyonwiki_code"  => data[2]
    #     )
    #     for data ∈ intData
    # ]
end

function vec_to_cflt( v::Vector{Any} )::ComplexF64
    v[1] + v[2]*1im
end

# numeric characters
function nchfromjs(js)::Union{Missing,Matrix{ComplexF64}}
    ncvecs = js["numeric_characters"]
    if ncvecs === nothing
        return missing
    else
        r = length(ncvecs)
        [ vec_to_cflt( ncvecs[i][j] ) for i in 1:r, j in 1:r ]
    end
end

# characters
function chfromjs(js)
    try
        vecs = js["characters"]
        string.(mapreduce( permutedims, vcat, vecs))
    catch e
        return missing
    end
end

# fpdims
function nfpdsfromjs(js)::Vector{ComplexF64}
    nfpdims = js["numeric_frobenius_perron_dimensions"]
    vec_to_cflt.( nfpdims )
end

# fpdim
function nfpdfromjs(js)::ComplexF64
    vec_to_cflt( js["numeric_frobenius_perron_dimension"] )
end

function cfromjs(js)
  # Known to be non categorifiable
  if js["categorifiable"] === false
    return false 
  end

  # Nothing known about categorifiability
  if js["categorifiable"] === nothing 
    return missing
  end

  # Has fusion categories
  return true

end

# TODO: only works for cats given by anyonwiki_code
# categorifications
function ctsfromjs(js)
    cats = js["categorifications"]

    if cats === nothing # Nothing known about categorifiability
        return missing
    elseif length(cats) === 0 # Known to have no cats
        return Vector{Int64}[]
    else 
        [ Int.(code) for code in cats ]
    end
end

function ctpfromjs(js)
    props = js["has_categories_with_props"]
end

# TODO: it should be possible to add type to output but I get the following error when importing FR^{2,10,0}_{1}:
# MethodError: Cannot `convert` an object of type Vector{Dict{String, Array}} to an object of type Dict{String, Array}
# The error is not reproducible when using the REPL
function npsrfromjs(js)#::Vector{Dict{String, Array}}
    try
    npsr = js["numeric_projective_SL2Z_reps"]
    if npsr == Any[]
        return Dict{String, Array}[]
    else
        dicts = Dict{String, Array}[]
        for rep in eachindex( npsr )
            sm = npsr[rep]["SMatrix"];
            tf = npsr[rep]["TwistFactors"];
            r  = size(sm,1);
            push!(
                dicts,
                Dict(
                "S_matrix"      =>
                    [ vec_to_cflt( sm[i][j]  ) for i in 1:r, j in 1:r ],
                "twist_factors" =>
                    [ [ vec_to_cflt( tf[i][j] ) for j in 1:r ] for i in 1:length(tf) ]
                )
            )
        end
        return dicts
    end
    catch e
        return missing
    end
end

# import names. Might fail
function nfromjs(js)
    try
        Vector{String}( js["names"] )
    catch e
        String[]
    end
end

# import texnames. Might fail
function tnfromjs(js)
    try
        Vector{String}( js["texnames"] )
    catch e
        String[]
    end
end

# import projective SL2Z reps
function psrfromjs(js)
    k = keys( js )
    if "projective_SL2Z_reps" ∈ k
        psr = js["projective_SL2Z_reps"]
    else
        return missing
    end

    if psr == "NotImplementedYet"
        return missing
    end

    if psr == Any[] 
        return Dict{String, Array}[]
    else
        dicts = Dict{String, Array}[]
        for rep in eachindex( psr )
            sm = psr[rep]["SMatrix"];
            tf = psr[rep]["TwistFactors"];
            r  = size(sm,1);
            push!(
                dicts,
                Dict(
                    "S_matrix"      =>
                        from_qqb_id(
                            [ sm[i][j] for i in 1:r, j in 1:r ]
                        ),
                    "twist_factors" =>
                        from_qqb_id(
                            [
                                [  vec[j] for j in 1:r ]
                                for vec in tf
                            ]
                        )
                )
            )
        end
        return dicts
    end
end


# import fpdim. Might fail
function fpdfromjs(js)
    try
        from_qqb_id(js["frobenius_perron_dimension"])
    catch e
        missing
    end
end

function fpdsfromjs(js)
    try
        from_qqb_id(js["frobenius_perron_dimensions"])
    catch e
        missing
    end
end

function ncrfromjs(js)
    k = keys( js )
    if "non_cat_reasons" ∈ k
        return js["non_cat_reasons"]
    else
        return Dict(
            "Fusion"    => "missing",
            "Pivotal"   => "missing",
            "Spherical" => "missing",
            "Unitary"   => "missing",
            "Braided"   => "missing",
            "Ribbon"    => "missing",
            "Modular"   => "missing"
        )
    end
end

export import_ring

function import_ring( filename::String )
    js = JSON.parsefile( filename );

    fusion_ring(
        mtfromjs( js ),
        names                               = nfromjs( js ),
        texnames                            = tnfromjs( js ),
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
        has_categories_with_props           = ctpfromjs( js ),
        #categorifiable                      = cfromjs( js ),
        categorifications                   = ctsfromjs( js ),
        references                          = js["references"],
        software                            = js["software"],
        comments                            = js["comments"],
        non_cat_reasons                     = ncrfromjs( js )
  )
end

export import_rings

function import_rings( filename::String )
    jsdict = JSON.parsefile( filename );

    frlist = FusionRing[]
    
    for ind in eachindex( jsdict["data"] )
        js = jsdict["data"][ind]
        fr = 
        fusion_ring(
            mtfromjs( js ),
            names                               = nfromjs( js ),
            texnames                            = tnfromjs( js ),
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
            has_categories_with_props           = ctpfromjs( js ),
            categorifications                   = ctsfromjs( js ),
            references                          = js["references"],
            software                            = js["software"],
            comments                            = js["comments"],
            non_cat_reasons                     = ncrfromjs( js )
        )
        push!( frlist, fr )
    end

    frlist
end



############################################################
# Exporting Fusion rings
############################################################

function missing_to_nothing(x)
    if x !== missing
        return x
    else
        return nothing
    end
end

function mttojs( fr::FusionRing )::Vector{Vector{Vector{Int64}}}
    mt = multiplication_table( fr )
    r  = rank( fr )
    [ [ [ mt[i,j,k] for k in 1:r ] for j in 1:r ] for i in 1:r ]
end

function chtojs( fr::FusionRing )
    ch = fr.characters
    if ch !== missing 
        r  = rank( fr )
        return [ [ ch[i,j] for j in 1:r ] for i in 1:r ]
    else
        return nothing
    end 
end

function sfrtojs( fr::FusionRing )
    missing_to_nothing( fr.sub_fusion_rings )
end

function psrtojs( fr::FusionRing )
    "NotImplementedYet" 
end

function fpdtojs( fr::FusionRing )::String
    qqb_id( fpdim( fr ) )
end

function fpdstojs( fr::FusionRing )::Vector{String}
    [ qqb_id( d ) for d in fpdims( fr ) ]
end

function tpdtojs( fr::FusionRing )
    missing_to_nothing( fr.tensor_product_decompositions )
end

function reim( x::ComplexF64 )::Vector{Float64}
    [ real(x), imag(x) ]
end

function reim( x::Float64 )::Vector{Float64}
    [ x, 0.0 ]
end

function reim( mat::Matrix{ComplexF64} )::Vector{Vector{Vector{Float64}}}
    [
        [
            reim( r[i] )
            for i in eachindex( r )
        ]
        for r in eachrow( mat )
    ]
end

function reim( vv::Vector{Vector{ComplexF64}} )::Vector{Vector{Vector{Float64}}}
    [ [ reim( coef ) for coef in vec ] for vec in vv ]
end

function nchtojs( fr::FusionRing )::Union{Vector{Vector{Vector{Float64}}},Nothing}
    if fr.numeric_characters === missing
        return nothing 
    else
        r = rank( fr )
        splitchars = reim.( numeric_characters( fr ) )
        return [ [ splitchars[i,j] for j in 1:r ] for i in 1:r ]
    end
end

function npsrtojs( fr::FusionRing )
    npsr = fr.numeric_projective_SL2Z_reps
    if npsr !== missing
        dicts = []
        for rep in npsr
            tf = reim(rep["twist_factors"])
            sm = reim(rep["S_matrix"])
            push!(
                dicts,
                Dict(
                    "twist_factors" => tf,
                    "S_matrix"      => sm
                )
            )
        end
        return dicts
    else
        return nothing
    end
end


function cpropstojs( fr::FusionRing )
    props = fr.has_categories_with_props
    function missing_to_nothing( v )
        if v[2] === missing
            [ v[1], nothing, v[3]  ]
        else
            v
        end
    end
    
    missing_to_nothing.(props)
end

function ctojs( fr::FusionRing )
    missing_to_nothing( is_categorifiable( fr ) )
end

function ctstojs( fr::FusionRing )
    missing_to_nothing( fr.categorifications )
end

function nfpdtojs( fr::FusionRing )
    reim( numeric_fpdim(fr) )
end

function nfpdstojs( fr::FusionRing )
    reim.( numeric_fpdims( fr ) )
end

function ncrtojs( fr::FusionRing )
    missing_to_nothing( fr.non_cat_reasons )
end

function write_json( filename::String, data::Dict )
    open( filename, "w" ) do f
        JSON.json( f, data, pretty = true, inline_limit = 10 )
    end 
end

function ring_to_dict( fr )
    infostring = "Fusion ring. mult_tab: structure constants. barcode & formal_code: unique identifiers see (DOI: 10.1063/5.0148848). non_trivial_sub_fusion_rings: tuples where the first element = elements of ring that form subring isomorphic to subring identified by second element of the tuple. software: doi of original software used to represent fusion ring. references: doi of paper from which data was obtained. categorifiable: false=not categorifiable, null= unknown. categorifications: if categorifiable then anyonwiki codes of pivotal (braided) fusion cats that categorify ring. numeric_projective_SL2Z_reps: each rep consists of a generalized S-matrix and a vector of vectors representing the ln(diag(T))/(2 pi i) of a generalized T-matrix. Algebraic numbers are encoded as a0_..._an__m where ai are polynomial coefficients and m is root number, ordered via Mathematica's convention."
    Dict(
        "mult_tab"                            => mttojs( fr )
       ,"names"                               => names( fr )
       ,"texnames"                            => tex_names( fr )
       ,"barcode"                             => string( barcode( fr ) )
       ,"anyonwiki_code"                      => anyonwiki_code( fr )
       ,"characters"                          => chtojs( fr )
       ,"non_trivial_sub_fusion_rings"        => sfrtojs(fr)
       ,"projective_SL2Z_reps"                => psrtojs( fr )
       ,"frobenius_perron_dimension"          => fpdtojs( fr )
       ,"frobenius_perron_dimensions"         => fpdstojs( fr )
       ,"tensor_product_decompositions"       => tpdtojs( fr )
       ,"numeric_characters"                  => nchtojs( fr )
       ,"numeric_projective_SL2Z_reps"        => npsrtojs( fr )
       ,"numeric_frobenius_perron_dimension"  => nfpdtojs( fr )
       ,"numeric_frobenius_perron_dimensions" => nfpdstojs( fr )
       ,"has_categories_with_props"           => cpropstojs( fr )
       ,"categorifications"                   => ctstojs( fr )
       ,"references"                          => fr.references
       ,"software"                            => fr.software
       ,"comments"                            => fr.comments
       ,"info"                                => infostring
    )
end

export rings_to_dict

function rings_to_dict( frs::Vector{FusionRing} )

    infostring = "Fusion ring. mult_tab: structure constants. barcode & formal_code: unique identifiers see (DOI: 10.1063/5.0148848). non_trivial_sub_fusion_rings: tuples (els,sr) with els = elements of ring that form subring isomorphic to ring sr. software: doi of original software used to represent fusion ring. references: doi of paper from which data was obtained. categorifications: if categorifiable then anyonwiki codes of pivotal (braided) fusion cats. numeric_projective_SL2Z_reps: each rep consists of a generalized S-matrix and a vector of vectors representing the ln(diag(T))/(2 pi i) of a generalized T-matrix. Algebraic numbers are encoded as a0_..._an__m where ai are polynomial coefficients and m is root number, ordered via Mathematica's convention. has_categories_with_props: triples [ prop, bool, reason ] where prop is the property, bool is true when its known at least one cat with prop exists, false when its known no cat with prop exists and null when no information is known. reason is a tuple [ method, str ] where method could be computer or theory and str gives more info."

    # We don't want to copy the infostring for each ring
    function ringtodict( fr )
        Dict(
         "mult_tab"                            => mttojs( fr )
        ,"names"                               => names( fr )
        ,"texnames"                            => tex_names( fr )
        ,"barcode"                             => string( barcode( fr ) )
        ,"anyonwiki_code"                      => anyonwiki_code( fr )
        ,"characters"                          => chtojs( fr )
        ,"non_trivial_sub_fusion_rings"        => sfrtojs(fr)
        ,"projective_SL2Z_reps"                => psrtojs( fr )
        ,"frobenius_perron_dimension"          => fpdtojs( fr )
        ,"frobenius_perron_dimensions"         => fpdstojs( fr )
        ,"tensor_product_decompositions"       => tpdtojs( fr )
        ,"numeric_characters"                  => nchtojs( fr )
        ,"numeric_projective_SL2Z_reps"        => npsrtojs( fr )
        ,"numeric_frobenius_perron_dimension"  => nfpdtojs( fr )
        ,"numeric_frobenius_perron_dimensions" => nfpdstojs( fr )
        ,"has_categories_with_props"           => cpropstojs( fr )
        ,"categorifications"                   => ctstojs( fr )
        ,"references"                          => fr.references
        ,"software"                            => fr.software
        ,"comments"                            => fr.comments
        )
    end
    
    Dict(
        "data" => Dict( fusion_ring_string(fr) => ringtodict(fr) for fr in frs ),
        "info" => infostring
    )
end

function fusion_ring_string(fr::FusionRing)
    c  = string.( anyonwiki_code( fr ) )
    us = fill( "_", 3 )
    stringriffle( c, us )
end

function fusion_ring_file_name(fr::FusionRing)
    fusion_ring_string(fr) * ".json"
end

export export_ring

function export_ring( filename::String,  fr::FusionRing )
    write_json( filename, ring_to_dict( fr ) )
end

export export_rings

function export_rings( filename::String, frs::Vector{FusionRing})
    write_json( filename, rings_to_dict( frs ) ) 
end

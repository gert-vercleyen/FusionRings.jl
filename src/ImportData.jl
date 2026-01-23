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
from_qqb_id( a::Array{String} ) = from_qqb_id.(a)

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
function fcfromjs( js::JSON.Object{String, Any} )::Vector{Int64}
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
function mtfromjs( js::JSON.Object{String, Any} )::Array{Int64, 3}
  jsmt = js["mult_tab"]
  r = length(jsmt)
  mt = zeros(Int, r, r, r)
  for i in 1:r, j in 1:r, k in 1:r 
      mt[i,j,k] = Int.(jsmt[i][j][k])
  end
  mt
end

# barcode
function bcfromjs(js::JSON.Object{String, Any})::ZZRingElem
  ZZ( parse( BigInt, js["barcode"] ) )
end

# tensor product decompositions
function tpdfromjs(js::JSON.Object{String, Any})
  tps = js["tensor_product_decompositions"]

    if length(tps) == 0
        []
    else
        tps = tps["value"]
        [
            [ Int.( code ) for code in decomp ]
            for decomp in tps
        ]
    end
end

# sub-fusion rings
function sfrfromjs(js::JSON.Object{String, Any})
    srs = js["non_trivial_sub_fusion_rings"]
    if length(srs) == 0
        return []
    end

    srs = srs["value"]

    intData =
        [
            [ Int.( vec ) for vec in sr ]
            for sr in srs
        ]

    [
        Dict(
            "injection"       => data[1],
            "anyonwiki_code"  => data[2]
        )
        for data ∈ intData
    ]
end

function vec_to_cflt( v::Vector{Any} )::ComplexF64
    v[1] + v[2]*1im
end

# numeric characters
function nchfromjs(js::JSON.Object{String, Any})::Matrix{ComplexF64}
    ncvecs = js["numeric_characters"]
    r = length(ncvecs)
    [ vec_to_cflt( ncvecs[i][j] ) for i in 1:r, j in 1:r ]
end

# characters
function chfromjs(js::JSON.Object{String, Any})
    try
        vecs = js["characters"]
        getfromqqbdict([ vec[i] for vec in vecs, i in eachindex(vec) ])
    catch e
        return missing
    end
end

# fpdims
function nfpdsfromjs(js::JSON.Object{String, Any})::Vector{ComplexF64}
    nfpdims = js["numeric_frobenius_perron_dimensions"]
    vec_to_cflt.( nfpdims )
end

# fpdim
function nfpdfromjs(js::JSON.Object{String, Any})::ComplexF64
    vec_to_cflt( js["numeric_frobenius_perron_dimension"] )
end

function cfromjs(js::JSON.Object{String, Any})
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
function ctsfromjs(js::JSON.Object{String, Any})
    # Known to be non categorifiable
    if js["categorifiable"] === false
        return Vector{Int64}[]
    end

    # Nothing known about categorifiability
    if js["categorifiable"] === nothing 
        return missing
    end

    # Has fusion categories
    cats = js["categorifications"]
    k    = keys( cats )

    # Legacy compatibility
    if "categories" ∈ k
        [ Int.(code) for code in cats["categories"] ]
    else
        [ Int.(code) for code in cats ]
    end
end

function ctpfromjs(js::JSON.Object{String, Any})
  if cfromjs(js) === missing
    return missing
  end

  if cfromjs(js) 
    props = js["has_categories_with_props"]

    Dict( 
      props[i][1] => props[i][2] 
      for i in 1:length(props) 
    )
  else
    Dict( "Fusion" => false )
  end
end

# TODO: it should be possible to add type to output but I get the following error when importing FR^{2,10,0}_{1}:
# MethodError: Cannot `convert` an object of type Vector{Dict{String, Array}} to an object of type Dict{String, Array}
# The error is not reproducible when using the REPL
function npsrfromjs(js::JSON.Object{String, Any})#::Vector{Dict{String, Array}}
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
end

# import names. Might fail
function nfromjs(js::JSON.Object{String, Any})
    try
        Vector{String}( js["names"] )
    catch e
        String[]
    end
end

# import texnames. Might fail
function tnfromjs(js::JSON.Object{String, Any})
    try
        Vector{String}( js["texnames"] )
    catch e
        String[]
    end
end

# import projective SL2Z reps
function psrfromjs(js::JSON.Object{String, Any})
    try
        psr = js["projective_SL2Z_reps"]
    catch e
        return missing
    end

    if psr == Any[]
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
                        getfromqqbdict(
                            [ sm[i][j] for i in 1:r, j in 1:r ]
                        ),
                    "twist_factors" =>
                        getfromqqbdict(
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
function fpdfromjs(js::JSON.Object{String, Any})
    try
        getfromqqbdict(js["frobenius_perron_dimension"])
    catch e
        missing
    end
end

function fpdsfromjs(js::JSON.Object{String, Any})
    try
        getfromqqbdict(js["frobenius_perron_dimensions"])
    catch e
        missing
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
        categorifiable                      = cfromjs( js ),
        categorifications                   = ctsfromjs( js ),
        references                          = js["references"],
        software                            = js["software"],
        comments                            = js["comments"]
  )
end


function export_ring( dir,  fr::FusionRing )

end

export load_frl

function load_frl()
    path = joinpath( @__DIR__, "data", "FusionRingsJSON" )

    prep_path( fn ) = joinpath( path, fn )

    filenames = prep_path.( readdir(path) )

    [ import_ring( fn ) for fn in filenames ]
end


# We want the mult-free rings to be first
frlsortcrit( c ) = anyonwiki_code(c)[ [ 2, 1, 3, 4 ] ]

export fusion_ring_list

fusion_ring_list = sort( load_frl(), by = frlsortcrit )


export frl

frl = fusion_ring_list


export fusion_ring_dict

fusion_ring_dict = Dict( anyonwiki_code(r) => r for r in frl )


export frd

frd = fusion_ring_dict


export anyonwiki_code

anyonwiki_code( r, m, nnsd, i ) = frd[ [ r, m, nnsd, i ] ]

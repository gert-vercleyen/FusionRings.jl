function import_ring(i::Int) 
  js = JSON.parsefile( joinpath(@__DIR__, "data","FusionRingsJSON", "ring"*string(i)*".json") )
  fc = [ js["formal_code"][i] for i in 1:4 ]
  r = fc[1]
  mt = zeros(Int, r, r, r)
  for i in 1:r, j in 1:r, k in 1:r 
      mt[i,j,k] = Int.(js["mult_tab"][i][j][k])
  end
  FusionRings.fusion_ring( mt, formal_code = fc)
end

export fusion_ring_list

fusion_ring_list =  [ import_ring(i) for i in 1:28451 ]

export frl

frl = fusion_ring_list

export fusion_ring_dict

fusion_ring_dict = Dict( anyonwiki_code(r) => r for r in frl )

export frd

frd = fusion_ring_dict

# The QQBarFieldElem objects are quite heavy to load so we will load them all in a dictionary 
# and provide functions to convert qqbar elems to keys and vice versa

# Generate unique ID for a QQBarFieldElem
function qqb_id( x::QQBarFieldElem ) 
    mp = minimal_polynomial(x)
    degstring = string( degree( mp ) )
    polstring = 
        replace( 
            string(mp),  
            "*" => "", " " => ""  
        )
    numstring = string( rootnum( x ) )

    degstring * "_" * polstring * "_" * numstring
    
end

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

function save_qqb_num( dir::String, x::QQBarFieldElem )
    idstring = qqb_id(x)
    Oscar.save( joinpath( dir, idstring * ".mrdi" ), idstring => x )
end

function save_qqb_num( x::QQBarFieldElem )
    save_qqb_num( joinpath(@__DIR__, "data", "Numbers", "QQBarFieldElems"), x )
end

# Load the dictionary of qqbar elems 


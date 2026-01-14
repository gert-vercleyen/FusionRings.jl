
# Formatting of fusion rings
function Base.show( io::IO, ring::FusionRing )
    p(str) = print( io, str );
    if ring.names != []
        p( "FR(" * names(ring)[1] * ")" )
    elseif ring.anyonwiki_code != []
        p( "FR(" * string(ring.anyonwiki_code)[2:end-1] * ")" )
    else
        props = map( string, comap( [ rank, multiplicity, nnsd ], ring ) )
        p( "FR(" * join( props, ", "  ) * ")" )
    end
end

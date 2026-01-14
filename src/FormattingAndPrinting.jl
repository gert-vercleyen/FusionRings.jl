
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

export print_multiplication_table

#TODO: 
# * no nlonger use names - should be labels
# * labels should be printed as bold integers
"""
Pretty prints the multiplication table as strings (no mutation).
"""
function print_multiplication_table(fr::FusionRing; include_zeros::Bool=false)
    N = multiplication_table(fr)
    names = fr.element_names
    r = length(names)
    head = "× │ " * join(names, " │ ")
    sep  = "──┼" * "───┼"^(r-1) * "──"
    println(head); println(sep)
    for i in 1:r
        rowcells = String[]
        for j in 1:r
            d = fusion_product(fr, i, j)
            if include_zeros
                parts = String[]
                for c in 1:r
                    m = get(d,c,0)
                    if m==0; push!(parts, "0 "*names[c])
                    elseif m==1; push!(parts, names[c])
                    else; push!(parts, string(m," ",names[c]))
                    end
                end
                push!(rowcells, join(parts, " + "))
            else
                isempty(d) && push!(rowcells, "0") && continue
                push!(rowcells,
                    join([ m==1 ? names[c] : string(m," ",names[c]) for (c,m) in d ], " + "))
            end
        end
        println(names[i], " │ ", join(rowcells, " │ "))
    end
    nothing
end

"Pretty one-liner: `a ⊗ b = ...` using printed names; `a,b` are indices."
function product_string(fr::FusionRing, a::Int, b::Int)
    rhs = let d = fusion_product(fr,a,b), names = fr.element_names
        isempty(d) ? "0" :
            join([ m==1 ? names[c] : string(m," ",names[c]) for (c,m) in d ], " ⊕ ")
    end
    string(fr.element_names[a], " × ", fr.element_names[b], " = ", rhs)
end
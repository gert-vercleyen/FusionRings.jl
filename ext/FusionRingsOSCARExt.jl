
module FusionRingsOSCARExt

using FusionRings
using Oscar
const GAP = Oscar.GAP

import FusionRings: fusion_ring
export group_fusion_ring, character_table_ring, song_extension_ring

function group_fusion_ring(G)
    elts = GAP.Globals.Elements(G)
    n = GAP.Globals.Length(elts)
    labstr = [String(GAP.gap_to_julia_string(GAP.Globals.String(e))) for e in elts]
    idx = Dict{Any,Int}(elts[i]=>i for i in 1:n)
    N = fill(0, n,n,n)
    for i in 1:n, j in 1:n
        prod = GAP.Globals.\*(elts[i], elts[j])
        k = idx[prod]
        N[i,j,k] = 1
    end
    fusion_ring(N; labels=labstr, name=string("GroupRing(", GAP.Globals.String(G), ")"))
end

function character_table_ring(G)
    T = GAP.Globals.CharacterTable(G)
    irrs = GAP.Globals.Irr(T)
    m = GAP.Globals.Length(irrs)
    N = fill(0, m,m,m)
    for i in 1:m, j in 1:m
        prod = GAP.Globals.\*(irrs[i], irrs[j])
        for k in 1:m
            mult = GAP.Globals.ScalarProduct(T, prod, irrs[k])
            N[i,j,k] = Int(mult)
        end
    end
    labs = ["χ"*string(i) for i in 1:m]
    fusion_ring(N; labels=labs, name="RepRing("*String(GAP.Globals.String(G))*")")
end

function song_extension_ring(G; n::Int=1)
    GAP.Globals.IsAbelian(G) == GAP.True || error("SONG defined here for abelian groups only")
    elts = GAP.Globals.Elements(G)
    nG = GAP.Globals.Length(elts)
    labs = [String(GAP.gap_to_julia_string(GAP.Globals.String(e))) for e in elts]
    labs = vcat(labs, ["X"])
    r = nG + 1
    N = fill(0, r,r,r)
    idx = Dict{Any,Int}(elts[i]=>i for i in 1:nG)
    for i in 1:nG, j in 1:nG
        prod = GAP.Globals.\*(elts[i], elts[j])
        k = idx[prod]
        N[i,j,k] = 1
    end
    X = r
    for i in 1:nG
        N[i,X,X] = 1
        N[X,i,X] = 1
    end
    for i in 1:nG
        N[X,X,i] = 1
    end
    if n>0
        N[X,X,X] = n
    end
    fusion_ring(N; labels=labs, name="SONG(" * String(GAP.Globals.String(G)) * "; n=$(n))")
end

end

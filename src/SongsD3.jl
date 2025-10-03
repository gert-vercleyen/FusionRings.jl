
module SongsD3

using ..Types: FusionRing
using ..Creation: fusion_ring
using ..Operations: multiplication_table, print_multiplication_table
using LinearAlgebra

export song_extension_D3, enumerate_song_extensions_D3, demo_song_D3

const D3 = ["1","r","r^2","s","sr","sr^2"]

"""
    parse_d3(x::String) -> (i::Int, flip::Bool)

Robust parser for D3 elements. Replaces previous chained ternary expression
to avoid parsing / style issues on some Julia versions and improve
readability.
"""
function parse_d3(x::String)::Tuple{Int,Bool}
    if x == "1";      return (0,false)
    elseif x == "r";  return (1,false)
    elseif x == "r^2";return (2,false)
    elseif x == "s";  return (0,true)
    elseif x == "sr"; return (1,true)
    elseif x == "sr^2"; return (2,true)
    else
        error("Unknown element $x")
    end
end

function unparse_d3(i::Int, flip::Bool)::String
    j = mod(i,3)
    if !flip
        return j==0 ? "1" : j==1 ? "r" : "r^2"
    else
        return j==0 ? "s" : j==1 ? "sr" : "sr^2"
    end
end

function d3mul(a::String, b::String)::String
    (i,fa) = parse_d3(a)
    (j,fb) = parse_d3(b)
    unparse_d3(i + (fa ? -j :  j), xor(fa,fb))
end

function inv_d3(x::String)::String
    for c in D3
        d3mul(x,c)=="1" && return c
    end
    error("No inverse for $x")
end

const H_trivial   = ["1"]
const H_rotations = ["1","r","r^2"]
const H_whole     = D3

same_coset(a::String, b::String, H::Vector{String}) = d3mul(inv_d3(a), b) in H

function coset_reps(G::Vector{String}, H::Vector{String})
    reps = String[]
    for g in G
        any(same_coset(g, r, H) for r in reps) || push!(reps, g)
    end
    reps
end

function automorphisms_on_factor(G::Vector{String}, H::Vector{String})
    idx = length(G) ÷ length(H)
    if idx == 1
        return [identity]
    elseif idx == 2
        A_id(x::String)    = x
        A_swap(x::String)  = x
        return [A_id, A_swap]
    else
        return [identity]
    end
end

build_orbit_set(G::Vector{String}, H::Vector{String}) = [r*"*t" for r in coset_reps(G,H)]
build_lift(T::Vector{String}) = Dict(x => replace(x,"*t"=>"") for x in T)
one_hot(x::String) = Dict(x=>1)
group_sum(H::Vector{String}) = Dict(h=>1 for h in H)

function conj_in_factor(tg::String, c::String, H::Vector{String}, reps::Vector{String})
    z = d3mul(inv_d3(tg), d3mul(c, tg))
    for r in reps
        same_coset(z,r,H) && return r
    end
    z
end

function check_A2_equals_conjugation(A::Function, tg::String, H::Vector{String}, reps::Vector{String})
    for c in reps
        lhs = A(A(c))
        lhsrep = begin
            found = c
            for r in reps; same_coset(lhs,r,H) && (found=r; break); end
            found
        end
        rhsrep = conj_in_factor(tg,c,H,reps)
        lhsrep == rhsrep || return false
    end
    true
end

function song_product(x::String, y::String;
    G::Vector{String}, H::Vector{String}, T::Vector{String},
    tildeg::String, n::Int, lift::Dict{String,String}, A::Function
)::Dict{String,Int}
    inGx = x in G; inGy = y in G
    inTx = x in T; inTy = y in T

    if inGx && inGy
        return one_hot(d3mul(x,y))
    elseif inGx && inTy
        return one_hot(d3mul(x, lift[y]) * "*t")
    elseif inTx && inGy
        return one_hot(d3mul(lift[x], y) * "*t")
    elseif inTx && inTy
        g1 = lift[x]; g2 = lift[y]
        rep   = d3mul(g1, A(g2))
        rep_t = d3mul(rep, inv_d3(tildeg))
        new_sum = Dict{String,Int}()
        for (h,coeff) in group_sum(H)
            p = d3mul(rep_t, h)
            new_sum[p] = get(new_sum,p,0) + coeff
        end
        if n != 0
            for t in T
                new_sum[t] = get(new_sum,t,0) + n
            end
        end
        return new_sum
    else
        error("Invalid SONG basis elements: $x, $y")
    end
end

function ringmul_to_fusionring(basis::Vector{String},
                               ringmul::Dict{Tuple{String,String},Dict{String,Int}};
                               name::String)
    labels = basis  # now Vector{String}
    r = length(labels)
    idx  = Dict{String,Int}(basis[i]=>i for i in 1:r)
    N = fill(0, r,r,r)
    for ((x,y), d) in ringmul
        a = idx[x]; b = idx[y]
        for (z,m) in d
            c = idx[z]
            N[a,b,c] += m
        end
    end
    fusion_ring(N; labels=labels, name=name)
end

function song_extension_D3(H; tildeg::String="s", n::Int=1, A=:id)
    Hs = H isa Symbol ? H : Symbol(H)
    As = A isa Symbol ? A : Symbol(A)
    Hset = Hs === :trivial   ? H_trivial   :
           Hs === :rotations ? H_rotations :
           Hs === :whole     ? H_whole     :
           error("H must be :trivial, :rotations, or :whole")
    G = D3
    T = build_orbit_set(G,Hset)
    L = build_lift(T)
    reps = coset_reps(G,Hset)
    auts = automorphisms_on_factor(G,Hset)
    Afunc = As === :id   ? auts[1] :
        As === :swap ? auts[min(2,length(auts))] :
            error("A must be :id or :swap")
    check_A2_equals_conjugation(Afunc, tildeg, Hset, reps) || @warn "A^2 != conj_{g̃} (continuing)"
    basis = vcat(G,T)
    ringmul = Dict{Tuple{String,String},Dict{String,Int}}()
    for x in basis, y in basis
        ringmul[(x,y)] = song_product(x,y; G=G,H=Hset,T=T,tildeg=tildeg,n=n,lift=L,A=Afunc)
    end
    nm = "SONG(D3; H=$(String(Hs)), g̃=$tildeg, n=$n, A=$(String(As)))"
    ringmul_to_fusionring(basis, ringmul; name=nm)
end

function enumerate_song_extensions_D3(; nset=(0,1), Aset=(:id,:swap))
    out = NamedTuples = NamedTuple[]
    G = D3
    results = NamedTuple[]
    for (Hsym,Hset) in ((:trivial,H_trivial), (:rotations,H_rotations), (:whole,H_whole))
        reps = coset_reps(G,Hset)
        auts = automorphisms_on_factor(G,Hset)
        T = build_orbit_set(G,Hset); L = build_lift(T)
        for tg in G
            for Asym in Aset
                Afunc = Asym===:id ? auts[1] : auts[min(2,length(auts))]
                check_A2_equals_conjugation(Afunc, tg, Hset, reps) || continue
                for n in nset
                    basis = vcat(G,T)
                    ringmul = Dict{Tuple{String,String},Dict{String,Int}}()
                    for x in basis, y in basis
                        ringmul[(x,y)] = song_product(x,y; G=G,H=Hset,T=T,tildeg=tg,n=n,lift=L,A=Afunc)
                    end
                    nm = "SONG(D3; H=$(String(Hsym)), g̃=$tg, n=$n, A=$(String(Asym)))"
                    push!(results, (H=Hsym, tildeg=tg, n=n, A=Asym,
                                    ring=ringmul_to_fusionring(basis, ringmul; name=nm)))
                end
            end
        end
    end
    results
end

function demo_song_D3(; n=1)
    R = song_extension_D3(:rotations; tildeg="s", n=n, A=:id)
    println("Built: ", R)
    print_tensor_table(R)
    R
end

end

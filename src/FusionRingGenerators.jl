module FusionRingGenerators
using ..Types, ..Creation, ..GeneralFunctions

export fibonacci_ring, ising_ring, semion_ring, zn_fusion_ring,
       su2_fusion_ring, psu2_fusion_ring, ty_fusion_ring, near_group_ring

# Z_n
function zn_fusion_ring(n::Int)
    mt = fill(0, n, n, n)
    for i in 0:n-1, j in 0:n-1
        k = mod(i+j, n)
        mt[i+1, j+1, k+1] = 1
    end
    labs = Symbol.(string.(0:n-1))
    Creation.fusion_ring(mt; labels=labs, names=["Z_"*string(n)], element_names=string.(labs))
end

# Fibonacci
function fibonacci_ring()
    labels = Symbol.(["1","τ"])
    N = fill(0,2,2,2)
    N[1,1,1]=1; N[1,2,2]=1; N[2,1,2]=1; N[2,2,1]=1; N[2,2,2]=1
    Creation.fusion_ring(N; labels=labels, names=["Fibonacci"], element_names=string.(labels))
end

# Ising
function ising_ring()
    labels = Symbol.(["1","σ","ψ"])
    N = fill(0,3,3,3)
    for i in 1:3; N[1,i,i]=1; N[i,1,i]=1; end
    N[2,2,1]=1; N[2,2,3]=1
    N[2,3,2]=1; N[3,2,2]=1
    N[3,3,1]=1
    Creation.fusion_ring(N; labels=labels, names=["Ising"], element_names=string.(labels))
end

# Semion
function semion_ring()
    labels = [:1,:s]
    N = fill(0,2,2,2)
    for i in 1:2; N[1,i,i]=1; N[i,1,i]=1; end
    N[2,2,1]=1
    Creation.fusion_ring(N; labels=labels, names=["Semion"], element_names=string.(labels))
end

# SU(2)_k
function su2_fusion_ring(k::Int)
    r = k+1
    N = fill(0,r,r,r)
    for a in 0:k, b in 0:k, c in 0:k
        if abs(a-b) <= c <= min(a+b, 2k - a - b) && iseven(a+b+c)
            N[a+1,b+1,c+1] = 1
        end
    end
    labs = Symbol.(string.(0:k))
    Creation.fusion_ring(N; labels=labs, names=["SU(2)"*GeneralFunctions.subscript_integer(k)], element_names=string.(labs))
end

# PSU(2)_k (even k)
range_psu2k(i,j,k) = abs(i-j):2:min(i+j, 2k - i - j)
function psu2_fusion_ring(k::Int)
    @assert iseven(k) "PSU(2)_k requires even k"
    rk = k ÷ 2 + 1
    mt = fill(0, rk, rk, rk)
    for a in 0:2:k, b in 0:2:k, c in 0:2:k
        if c in range_psu2k(a,b,k)
            mt[a÷2+1, b÷2+1, c÷2+1] = 1
        end
    end
    labs = Symbol.(string.(0:rk-1))
    Creation.fusion_ring(mt; labels=labs, names=["PSU(2)"*GeneralFunctions.subscript_integer(k)], element_names=string.(labs))
end

# Tambara–Yamagami for abelian group given as a list G
function ty_fusion_ring(G::AbstractVector)
    n = length(G)
    rank = n+1
    mt = fill(0, rank, rank, rank)
    for i in 1:n, j in 1:n
        k = (i + j - 2) % n + 1
        mt[i,j,k] = 1
    end
    m = rank
    for i in 1:n
        mt[i,m,m]=1; mt[m,i,m]=1
    end
    for i in 1:n
        mt[m,m,i]=1
    end
    labs = Symbol.(vcat(string.(G), ["m"]))
    Creation.fusion_ring(mt; labels=labs, names=["TY(" * join(string.(G),",") * ")"], element_names=string.(labs))
end

# Near-group ring over abelian group with one extra object X: X⊗X = ⊕g g ⊕ n X
function near_group_ring(G::AbstractVector{<:AbstractString}; n::Int=0)
    nG = length(G)
    rank = nG + 1
    mt = fill(0, rank, rank, rank)
    for i in 1:nG, j in 1:nG
        k = (i + j - 2) % nG + 1
        mt[i,j,k] = 1
    end
    X = rank
    for i in 1:nG
        mt[i,X,X]=1; mt[X,i,X]=1
    end
    for i in 1:nG
        mt[X,X,i]=1
    end
    mt[X,X,X]=n
    labs = Symbol.(vcat(G, ["X"]))
    Creation.fusion_ring(mt; labels=labs, names=["NearGroup(" * join(G,",") * "; n=$(n))"], element_names=string.(labs))
end

end # module

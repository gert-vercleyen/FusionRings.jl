
module FusionRingGenerators

using ..Types: FusionRing
using ..Creation: fusion_ring

export fibonacci_ring, ising_ring, semion_ring
export su2_fusion_ring, psu2_fusion_ring, zn_fusion_ring
export ty_fusion_ring, near_group_ring

_sym(x) = x isa Symbol ? x : Symbol(x)

function fibonacci_ring()
    labels = Symbol["1", :τ]
    N1 = [1 0; 0 1]
    Nτ = [0 1; 1 1]
    N = Array{Int,3}(undef, 2,2,2)
    N[1,:,:] = N1
    N[2,:,:] = Nτ
    fusion_ring(N; labels=labels, name="Fibonacci")
end

function ising_ring()
    labels = Symbol["1", :σ, :ψ]
    N1 = [1 0 0; 0 1 0; 0 0 1]
    Nσ = [0 1 0; 1 0 1; 0 1 0]
    Nψ = [0 0 1; 0 1 0; 1 0 0]
    N = Array{Int,3}(undef, 3,3,3)
    N[1,:,:] = N1
    N[2,:,:] = Nσ
    N[3,:,:] = Nψ
    fusion_ring(N; labels=labels, name="Ising")
end

function semion_ring()
    labels = Symbol["1", :s]
    N1 = [1 0; 0 1]
    Ns = [0 1; 1 0]
    N = Array{Int,3}(undef, 2,2,2)
    N[1,:,:] = N1
    N[2,:,:] = Ns
    fusion_ring(N; labels=labels, name="Semion")
end

function su2_fusion_ring(k::Int)
    r = k+1
    labs = Symbol.(string.(0:k))
    N = fill(0, r,r,r)
    for a in 0:k, b in 0:k, c in 0:k
        if abs(a-b) ≤ c ≤ min(a+b, 2k - a - b) && ((a+b+c) % 2 == 0)
            N[a+1, b+1, c+1] = 1
        end
    end
    fusion_ring(N; labels=labs, name="SU(2)_$(k)")
end

function psu2_fusion_ring(k::Int)
    k % 2 == 0 || error("PSU(2)_k only defined for even k")
    su = su2_fusion_ring(k)
    even = collect(1:2:k+1)
    r2 = length(even)
    N2 = fill(0, r2,r2,r2)
    for i in 1:r2, j in 1:r2, l in 1:r2
        N2[i,j,l] = fusion_tensor(su)[even[i], even[j], even[l]]
    end
    labs = Symbol.(string.(0:r2-1))
    fusion_ring(N2; labels=labs, name="PSU(2)_$(k)")
end

function zn_fusion_ring(n::Int)
    labs = Symbol.(string.(0:n-1))
    N = fill(0, n,n,n)
    for i in 0:n-1, j in 0:n-1
        k = (i + j) % n
        N[i+1, j+1, k+1] = 1
    end
    fusion_ring(N; labels=labs, name="Z_$(n)")
end

function ty_fusion_ring(G::Vector{String})
    n = length(G); r = n+1
    labs = Symbol.(vcat(G, ["m"]))
    N = fill(0, r,r,r)
    for i in 1:n, j in 1:n
        k = ((i-1) + (j-1)) % n + 1
        N[i, j, k] = 1
    end
    m = r
    for i in 1:n
        N[i, m, m] = 1
        N[m, i, m] = 1
    end
    for i in 1:n
        N[m, m, i] = 1
    end
    fusion_ring(N; labels=labs, name="TY(" * join(G,",") * ")")
end

function near_group_ring(G::Vector{String}; n::Int=1)
    n ≥ 0 || error("n must be ≥ 0")
    nG = length(G)
    labs = Symbol.(vcat(G, ["X"]))
    r = nG+1
    N = fill(0, r,r,r)
    for i in 1:nG, j in 1:nG
        k = ((i-1) + (j-1)) % nG + 1
        N[i, j, k] = 1
    end
    X = r
    for i in 1:nG
        N[i, X, X] = 1
        N[X, i, X] = 1
    end
    for i in 1:nG
        N[X, X, i] = 1
    end
    if n>0
        N[X, X, X] = n
    end
    fusion_ring(N; labels=labs, name="NearGroup(" * join(G,",") * "; n=$(n))")
end

end

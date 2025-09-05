module Properties
using LinearAlgebra
using ..Types, ..Operations

export quantum_dimensions, global_dimension, is_commutative,
       is_equivalent, is_sub_fusion_ring, sub_fusion_rings

# Compute Frobenius–Perron dimensions via common positive eigenvector
function quantum_dimensions(fr::Types.FusionRing)
    r = size(fr.N,1)
    # positive combination of fusion matrices
    M = zeros(Float64, r, r)
    for a in 1:r
        M .+= Operations.fusion_matrix(fr,a)
    end
    vals, vecs = eigen(Symmetric(M))
    idx = argmax(vals)
    v = abs.(vecs[:,idx]) .+ eps()
    d = zeros(Float64,r)
    for a in 1:r
        w = Operations.fusion_matrix(fr,a) * v
        d[a] = w[1]/v[1]
    end
    return d
end

global_dimension(fr::Types.FusionRing) = sum(quantum_dimensions(fr).^2)

function is_commutative(fr::Types.FusionRing)
    r = size(fr.N,1)
    @inbounds for a in 1:r, b in 1:r, c in 1:r
        fr.N[a,b,c] == fr.N[b,a,c] || return false
    end
    true
end

# Simple permutations generator to avoid external deps
function _permutations(v::Vector{Int})
    if isempty(v); return Vector{Vector{Int}}([Int[]]); end
    out = Vector{Vector{Int}}()
    function backtrack(a, rest)
        if isempty(rest); push!(out, copy(a)); return; end
        for i in eachindex(rest)
            push!(a, rest[i])
            backtrack(a, vcat(rest[1:i-1], rest[i+1:end]))
            pop!(a)
        end
    end
    backtrack(Int[], v)
    return out
end

function is_equivalent(fr1::Types.FusionRing, fr2::Types.FusionRing)
    size(fr1.N) == size(fr2.N) || return false
    fr1.N == fr2.N && return true
    r = size(fr1.N,1)
    for p in _permutations(collect(1:r))
        ok = true
        @inbounds for a in 1:r, b in 1:r, c in 1:r
            if fr1.N[a,b,c] != fr2.N[p[a], p[b], p[c]]
                ok = false; break
            end
        end
        ok && return true
    end
    return false
end

function is_sub_fusion_ring(fr::Types.FusionRing, subset)
    labs = Symbol[]
    for s in subset
        push!(labs, s isa Int ? fr.labels[s] : Symbol(s))
    end
    inds = [findfirst(==(l), fr.labels) for l in labs]
    S = Set(inds)
    r = size(fr.N,1)
    for a in S, b in S, c in 1:r
        (fr.N[a,b,c]>0 && !(c in S)) && return false
    end
    true
end

function sub_fusion_rings(fr::Types.FusionRing; proper::Bool=true)
    r = size(fr.N,1)
    out = Vector{Vector{Symbol}}()
    for mask in 1:(2^r-1)
        proper && mask==(2^r-1) && continue
        subset = [i for i in 1:r if (mask >> (i-1)) & 1 == 1]
        is_sub_fusion_ring(fr, subset) || continue
        push!(out, fr.labels[subset])
    end
    out
end

end # module

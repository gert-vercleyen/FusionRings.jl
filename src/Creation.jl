
module Creation

using ..Types: FusionRing
using LinearAlgebra

export fusion_ring

_nonnegints(A) = all(x -> x isa Integer && x ≥ 0, A)
_cubic3(A) = ndims(A)==3 && size(A,1)==size(A,2)==size(A,3)

function _unit_ok(N::Array{Int,3})
    r = size(N,1)
    for a in 1:r, b in 1:r
        if N[1,a,b] != (a==b ? 1 : 0); return false; end
        if N[a,1,b] != (a==b ? 1 : 0); return false; end
    end
    return true
end

function _associative(N::Array{Int,3})
    r = size(N,1)
    for a in 1:r, b in 1:r, c in 1:r, d in 1:r
        lhs = 0; rhs = 0
        for e in 1:r
            lhs += N[a,b,e]*N[e,c,d]
            rhs += N[a,e,d]*N[b,c,e]
        end
        lhs==rhs || return false
    end
    return true
end

function fusion_ring(N::Array{Int,3}; labels::Union{Nothing,Vector{Symbol},Vector{String}}=nothing, name::String="FusionRing")
    !_nonnegints(N) && error("Structure constants must be nonnegative integers.")
    !_cubic3(N) && error("Tensor must be 3-dimensional and cubic.")
    !_unit_ok(N) && error("First simple must act as the tensor unit.")
    !_associative(N) && error("Tensor fails associativity.")
    if labels === nothing
        labels = [i==1 ? "1" : string(i-1) for i in 1:size(N,1)]
    elseif eltype(labels) <: Symbol
        labels = String.(labels)
    end
    length(labels)==size(N,1) || error("labels length must equal rank.")
    return FusionRing(N, labels, name)
end

end

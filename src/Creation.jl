\
module Creation
using LinearAlgebra
using ..Types
using ..GeneralFunctions

export fusion_ring, is_fusion_ring, _idx

@inline function _idx(fr::Types.FusionRing, x)
    x isa Int && return x
    i = get(fr.dict, Symbol(x), 0)
    i == 0 && throw(ArgumentError("Unknown label $(x)"))
    return i
end

# Validation helpers
_nonneg_ints(N) = all(x->x>=0, N)
_cubic(N) = (ndims(N)==3) && (size(N,1)==size(N,2)==size(N,3))

function _unit_ok(N::Array{Int,3})
    r = size(N,1)
    for b in 1:r, c in 1:r
        if N[1,b,c] != (b==c ? 1 : 0); return false; end
        if N[b,1,c] != (b==c ? 1 : 0); return false; end
    end
    true
end

function _duals_exist(N::Array{Int,3})
    r = size(N,1)
    for a in 1:r
        found = false
        for b in 1:r
            if N[a,b,1]==1 && N[b,a,1]==1
                found = true; break
            end
        end
        found || return false
    end
    true
end

function _associative(N::Array{Int,3})
    r = size(N,1)
    @inbounds for a in 1:r, b in 1:r, c in 1:r, d in 1:r
        lhs = 0
        rhs = 0
        for e in 1:r
            lhs += N[a,b,e]*N[e,c,d]
            rhs += N[a,e,d]*N[b,c,e]
        end
        lhs == rhs || return false
    end
    true
end

"""
    fusion_ring(N; labels, names, element_names)

Build a validated `FusionRing` from a 3-tensor `N` of nonnegative integers.
The first label is the tensor unit.
"""
function fusion_ring(N::Array{<:Integer,3};
    labels::Union{Nothing,Vector{Symbol}}=nothing,
    names::Union{Nothing,Vector{String}}=nothing,
    element_names::Union{Nothing,Vector{String}}=nothing)

    Nint = Int.(N)
    !_nonneg_ints(Nint)  && error("All structure constants must be nonnegative integers")
    !_cubic(Nint)        && error("Tensor must be rank-3 with equal side lengths")
    !_unit_ok(Nint)      && error("First label must be the tensor unit")
    !_duals_exist(Nint)  && error("Some objects are missing two-sided duals")
    !_associative(Nint)  && error("Fusion tensor is not associative")

    r = size(Nint,1)
    labs = labels === nothing ? Symbol.(string.(1:r)) : labels
    eln  = element_names === nothing ? string.(labs) : element_names
    nms  = names === nothing ? String[] : names
    fr = Types.FusionRing(labs, Nint, nms, eln, GeneralFunctions.indexdict(labs))
    return fr
end

is_fusion_ring(N::Array{<:Integer,3}) = try
    fusion_ring(N); true
catch
    false
end

end # module

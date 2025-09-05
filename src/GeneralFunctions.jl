\
module GeneralFunctions
export bold_integer, subscript_integer, indexdict, permute_tensor, tensor_equal

# Simple helpers for pretty names
const _bold_digits = Dict('0'=>'𝟎','1'=>'𝟏','2'=>'𝟐','3'=>'𝟑','4'=>'𝟒','5'=>'𝟓','6'=>'𝟔','7'=>'𝟕','8'=>'𝟖','9'=>'𝟗')
bold_integer(n::Integer) = join(get(_bold_digits, c, c) for c in string(n))

const _sub_digits = Dict('0'=>'₀','1'=>'₁','2'=>'₂','3'=>'₃','4'=>'₄','5'=>'₅','6'=>'₆','7'=>'₇','8'=>'₈','9'=>'₉')
subscript_integer(n::Integer) = join(get(_sub_digits, c, c) for c in string(n))

indexdict(labels::Vector{Symbol}) = Dict(l=>i for (i,l) in enumerate(labels))

function permute_tensor(N::Array{<:Integer,3}, p::AbstractVector{<:Integer})
    r = length(p)
    M = similar(N)
    @inbounds for a in 1:r, b in 1:r, c in 1:r
        M[p[a],p[b],p[c]] = N[a,b,c]
    end
    return M
end

tensor_equal(A::Array{<:Integer,3}, B::Array{<:Integer,3}) = all(A .== B)

end # module

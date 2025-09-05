module Types
export FusionRing, labels, rank, fusion_tensor

mutable struct FusionRing
    labels::Vector{Symbol}
    N::Array{Int,3}
    names::Vector{String}
    element_names::Vector{String}
    dict::Dict{Symbol,Int}
end

labels(fr::FusionRing) = fr.labels
rank(fr::FusionRing) = length(fr.labels)
fusion_tensor(fr::FusionRing) = fr.N

end # module

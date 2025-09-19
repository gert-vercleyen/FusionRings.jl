
module Types

export FusionRing, labels, rank, fusion_tensor, ringname

Base.@kwdef struct FusionRing
    N::Array{Int,3}
    labels::Vector{Symbol}
    name::String = "FusionRing"
end

labels(fr::FusionRing) = fr.labels
rank(fr::FusionRing) = size(fr.N, 1)
fusion_tensor(fr::FusionRing) = fr.N
ringname(fr::FusionRing) = fr.name

end

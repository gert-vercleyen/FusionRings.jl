
module Types

export FusionRing, labels, rank, fusion_tensor, ringname

Base.@kwdef struct FusionRing
    N::Array{Int,3}
    labels::Vector{String}
    name::String = "FusionRing"
end

# Backward compatibility accessor allowing Symbols in legacy code paths
_sym2str(x) = x isa Symbol ? String(x) : String(x)

labels(fr::FusionRing) = fr.labels
rank(fr::FusionRing) = size(fr.N, 1)
fusion_tensor(fr::FusionRing) = fr.N
ringname(fr::FusionRing) = fr.name

end


module PrettyPrint

using ..Types: FusionRing, labels, rank, ringname

function Base.show(io::IO, fr::FusionRing)
    print(io, string(ringname(fr), " with ", rank(fr), " simples: ", collect(labels(fr))))
end

end

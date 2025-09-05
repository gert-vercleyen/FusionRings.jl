module PrettyPrint
using ..Types
import Base: show, getindex

getindex(R::Types.FusionRing, i::Int) = R.element_names[i]

function show(io::IO, ::MIME"text/plain", R::Types.FusionRing)
    nm = isempty(R.names) ? "FR(" * string(length(R.labels)) * ")" : R.names[1]
    print(io, nm, " with ", length(R.labels), " simples: ", R.labels)
end

end # module

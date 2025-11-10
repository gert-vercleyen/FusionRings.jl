
module FusionRingsBSONExt

using FusionRings
using BSON

function FusionRings.DataCache._bson_read(path::AbstractString)
    isfile(path) || return nothing
    d = BSON.load(path)
    haskey(d, :data) ? d[:data] : d
end

end

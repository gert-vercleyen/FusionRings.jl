
module FusionRings

using Oscar
import Oscar: multiplication_table, is_commutative
using Combinatorics
using JSON

include("GeneralFunctions.jl")
include("Creation.jl")
include("Properties.jl")
include("Operations.jl")
include("ImportData.jl")
end

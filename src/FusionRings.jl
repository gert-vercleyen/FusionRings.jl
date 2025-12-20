module FusionRings

using Oscar
import Oscar: multiplication_table, is_commutative,rank
using Combinatorics
using JSON
using LinearAlgebra:eigen

include("GeneralFunctions.jl")
include("Creation.jl")
include("Operations.jl")
include("ImportData.jl")
end

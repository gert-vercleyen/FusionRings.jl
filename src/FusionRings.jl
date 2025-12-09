module FusionRings

include("Types.jl")
include("Creation.jl")
include("Operations.jl")
include("Properties.jl")
# include other files as needed (Generators, ModularData, …)

using .Types
using .Creation
using .Operations
using .Properties

export FusionRing, fusion_ring
export multiplication_table
export fusion_matrix, fusion_coeff, fusion_product, decompose, decompose_all
export product_string, print_multiplication_table, is_equivalent, permute
export numeric_frobenius_perron_dimensions, nfpdims
export is_commutative, conjugate_element
export is_sub_fusion_ring, sub_fusion_rings, sub_fusion_rings_elements
export commutator

end # module

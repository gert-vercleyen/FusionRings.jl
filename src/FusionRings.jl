module FusionRings

include("GeneralFunctions.jl")
include("Types.jl")
include("Creation.jl")
include("Operations.jl")
include("Properties.jl")
include("FusionRingGenerators.jl")
include("PrettyPrint.jl")
include("DataCache.jl")

# Bring submodule APIs into the top-level namespace and re-export
using .Types: FusionRing, labels, rank, fusion_tensor
using .Creation: fusion_ring, is_fusion_ring
using .Operations: fusion_matrix, fusion_coeff, tensor_product, decompose, decompose_all, tensor_table, print_tensor_table
using .Properties: quantum_dimensions, global_dimension, is_commutative, is_equivalent, is_sub_fusion_ring, sub_fusion_rings
using .FusionRingGenerators: fibonacci_ring, ising_ring, semion_ring, zn_fusion_ring, su2_fusion_ring, psu2_fusion_ring, ty_fusion_ring, near_group_ring
using .DataCache: FusionRingByCode, FRBC, FusionRingList, FRL, all_fusion_rings_q, AFRQ

export FusionRing, fusion_ring, is_fusion_ring, labels, rank
export fusion_tensor, fusion_matrix, fusion_coeff, tensor_product, decompose, decompose_all, tensor_table, print_tensor_table
export quantum_dimensions, global_dimension, is_commutative
export is_equivalent, is_sub_fusion_ring, sub_fusion_rings
export fibonacci_ring, ising_ring, semion_ring, zn_fusion_ring
export su2_fusion_ring, psu2_fusion_ring, ty_fusion_ring, near_group_ring
export FusionRingByCode, FRBC, FusionRingList, FRL, all_fusion_rings_q, AFRQ

end # module

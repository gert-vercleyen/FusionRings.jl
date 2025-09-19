module FusionRings

using LinearAlgebra
using Serialization

include("GeneralFunctions.jl")
include("Types.jl")
include("PrettyPrint.jl")

include("Creation.jl")
include("Operations.jl")
include("Properties.jl")

include("FusionRingGenerators.jl")

include("DataCache.jl")
include("Registry.jl")

include("SongsD3.jl")

export FusionRing, labels, rank, fusion_tensor, ringname

export fusion_ring

export fusion_matrix, fusion_coeff, tensor_product, decompose, decompose_all
export tensor_table, print_tensor_table, product_string
export permute_mult_tab, permute, is_equivalent

export quantum_dimensions, global_dimension, is_commutative, multiplicity
export nonzero_structure_constants, conjugation_matrix, is_group_ring
export conjugate_element, sub_fusion_rings, is_sub_fusion_ring

export fibonacci_ring, ising_ring, semion_ring, su2_fusion_ring, psu2_fusion_ring
export zn_fusion_ring, ty_fusion_ring, near_group_ring

export optimized_import, FusionRingList, FRL, FusionRingByCode, FRBC, AllFusionRingsQ, AFRQ

export save_ring_json, load_registry_json, query_registry_by_fpdims

export song_extension_D3, enumerate_song_extensions_D3, demo_song_D3

end

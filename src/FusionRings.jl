module FusionRings

using Serialization

include("Types.jl")
include("GeneralFunctions.jl")
include("PrettyPrint.jl")

include("Creation.jl")
include("Operations.jl")
include("Properties.jl")

include("NamedFusionRings.jl")
include("FusionRingGenerators.jl") # deprecated stub

include("DataCache.jl")
include("Registry.jl")

include("SongsD3.jl")

# Re-export selected symbols from internal submodules so that
# `using FusionRings` provides the public API directly. Each included file
# defines a submodule, therefore we explicitly import the desired names.
using .Types: FusionRing, labels, rank, fusion_tensor, ringname
using .GeneralFunctions: indexmap
using .Creation: fusion_ring
using .Operations: fusion_matrix, fusion_coeff, tensor_product, decompose, decompose_all,
				   multiplication_table, print_multiplication_table, product_string,
				   permute_mult_tab, permute, is_equivalent
using .Properties: fpdims, fpdim, is_commutative, multiplicity,
				   nonzero_structure_constants, conjugation_matrix, is_group_ring,
				   conjugate_element, conjugate_label, sub_fusion_rings, is_sub_fusion_ring
using .NamedFusionRings: fibonacci_ring, ising_ring, semion_ring, su2_fusion_ring,
				 psu2_fusion_ring, zn_fusion_ring, ty_fusion_ring, near_group_fusion_ring
using .Registry: save_ring_json, load_registry_json, query_registry_by_fpdims
using .SongsD3: song_extension_D3, enumerate_song_extensions_D3, demo_song_D3

export FusionRing, labels, rank, fusion_tensor, ringname

export fusion_ring

export fusion_matrix, fusion_coeff, tensor_product, decompose, decompose_all
export multiplication_table, print_multiplication_table, product_string
export permute_mult_tab, permute, is_equivalent

export fpdims, fpdim, is_commutative, multiplicity
export nonzero_structure_constants, conjugation_matrix, is_group_ring
export conjugate_element, conjugate_label, sub_fusion_rings, is_sub_fusion_ring

export fibonacci_ring, ising_ring, semion_ring, su2_fusion_ring, psu2_fusion_ring
export zn_fusion_ring, ty_fusion_ring, near_group_fusion_ring

export optimized_import, FusionRingList, FRL, FusionRingByCode, FRBC, AllFusionRingsQ, AFRQ

export save_ring_json, load_registry_json, query_registry_by_fpdims

export song_extension_D3, enumerate_song_extensions_D3, demo_song_D3

# Deprecations (soft aliases) ----------------------------------------------
@deprecate tensor_table multiplication_table
@deprecate print_tensor_table print_multiplication_table
@deprecate quantum_dimensions fpdims
@deprecate global_dimension fpdim
@deprecate near_group_ring near_group_fusion_ring

end

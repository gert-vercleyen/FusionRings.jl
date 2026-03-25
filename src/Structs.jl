export FusionRing

# TODO:
# * add upper central series
# * make non_cat_reason a list with reasons for each prop
#   so
#   [
#     [ non_pivotal, reason_np ],
#     [ non_unitary, reason_nu ],
#     ...
#   ]
# * grading groups 

struct FusionRing
    multiplication_table::Array{Int,3}
    names::Array{String,1}
    texnames::Array{String,1}
    labels::Array{String,1}
    barcode
    anyonwiki_code::Union{Array{Int,1},Missing}
    characters
    sub_fusion_rings
    projective_SL2Z_reps
    frobenius_perron_dimension
    frobenius_perron_dimensions
    tensor_product_decompositions
    numeric_characters
    numeric_projective_SL2Z_reps
    numeric_frobenius_perron_dimension
    numeric_frobenius_perron_dimensions
    has_categories_with_props
    categorifiable
    categorifications
    references
    software
    comments
    non_cat_reasons
end

export fusion_ring

check_struct_const(mt) = all(x -> x isa Integer && x >= 0, mt)

function check_mt_dims(mt)
    dims = size(mt)
    length(dims) == 3 && is_constant_array(dims)
end

function check_unit(mt)
    r = size(mt)[1]
    δ(i,j) = i == j ? 1 : 0
    for i in 1:r, j in 1:r
        if !(mt[1, i, j] == mt[i, 1, j] == δ(i,j))
            return false
        end
        continue
    end
    return true
end

check_inverse(mt) = sum(mt[:, :, 1]) == size(mt, 1)

function check_associativity(mt::Array{Int, 3})
    r = size(mt, 1)
    for a in 1:r, b in 1:r, c in 1:r, d in 1:r
        lhs = sum(mt[a, b, e] * mt[e, c, d] for e in 1:r)
        rhs = sum(mt[a, f, d] * mt[b, c, f] for f in 1:r)
        lhs == rhs || return false
    end
    true
end

check_labels(mt, names) = length(names) == size(mt, 1)


function fusion_ring(
    mt;
    labels                              = [],
    names                               = [], 
    texnames                            = [],
    barcode                             = missing,
    anyonwiki_code                      = missing,
    characters                          = missing,
    sub_fusion_rings                    = missing,
    projective_SL2Z_reps                = missing,
    frobenius_perron_dimension          = missing,
    frobenius_perron_dimensions         = missing,
    tensor_product_decompositions       = missing,
    numeric_characters                  = missing,
    numeric_frobenius_perron_dimension  = missing,
    numeric_frobenius_perron_dimensions = missing,
    numeric_projective_SL2Z_reps        = missing,
    has_categories_with_props           = missing,
    categorifiable                      = missing,
    categorifications                   = missing,
    references                          = missing,
    software                            = missing,
    comments                            = missing,
    non_cat_reasons                     = missing,
    skip_check                          = false,
    )

    if !skip_check
        check_struct_const(mt)     || error("All structure constants must be non-negative integers")
        check_mt_dims(mt)          || error("multiplication_table must be a 3-tensor with equal side lengths")
        check_unit(mt)             || error("First basis element must act as unit object")
        check_inverse(mt)          || error("Each simple object must have a unique inverse")
        check_associativity(mt)    || error("Structure constants violate associativity")
        (labels == [] || check_labels(mt, labels)) ||
            error("labels length ≠ rank")
    end

    labels == [] && (labels = String[bold_integer(i) for i in 1:size(mt, 1)])

    FusionRing(
        mt
        ,names
        ,texnames
        ,labels
        ,barcode
        ,anyonwiki_code
        ,characters
        ,sub_fusion_rings
        ,projective_SL2Z_reps
        ,frobenius_perron_dimension
        ,frobenius_perron_dimensions
        ,tensor_product_decompositions
        ,numeric_characters
        ,numeric_projective_SL2Z_reps
        ,numeric_frobenius_perron_dimension
        ,numeric_frobenius_perron_dimensions
        ,has_categories_with_props
        ,categorifiable
        ,categorifications
        ,references
        ,software
        ,comments
        ,non_cat_reasons
    )
end
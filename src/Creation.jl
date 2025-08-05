include("GeneralFunctions.jl")

using LinearAlgebra   # for `I`

struct FusionRing
    multiplication_table
    names
    texnames
    element_names::Vector{String}
    barcode
    formal_code
    direct_product_decompositions
    sub_fusion_rings
    frobenius_perron_dimensions
    modular_data
    characters
end

export fusion_ring

check_struct_const(mt) = all(x -> x >= 0, mt)

function check_mt_dims(mt)
    dims = size(mt)
    length(dims) == 3 && is_constant_array(dims)
end

check_unit(mt) = mt[1, :, :] == mt[:, 1, :] == I

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

check_element_names(mt, names) = length(names) == size(mt, 1)


function fusion_ring(mt; names = missing, texnames = missing, element_names = missing,
                     barcode = missing, formal_code = missing,
                     tensor_product_decompositions = missing, sub_fusion_rings = missing,
                     frobenius_perron_dimensions = missing, modular_data = missing,
                     characters = missing)

    check_struct_const(mt)     || error("All structure constants must be non‑negative integers")
    check_mt_dims(mt)          || error("multiplication_table must be a 3‑tensor with equal side lengths")
    check_unit(mt)             || error("First basis element must act as unit object")
    check_inverse(mt)          || error("Each simple object must have a unique inverse")
    check_associativity(mt)    || error("Structure constants violate associativity")
    (element_names === missing || check_element_names(mt, element_names)) ||
        error("element_names length ≠ rank")

    element_names === missing && (element_names = [bold_integer(i) for i in 1:size(mt, 1)])

    FusionRing(
        Int.(mt), names, texnames, element_names, barcode, formal_code,
        tensor_product_decompositions, sub_fusion_rings, frobenius_perron_dimensions,
        modular_data, characters)
end

export psu2k_fusion_ring, su2k_fusion_ring, son2_fusion_ring, metaplectic_fusion_ring,
       fusion_ring_from_group, zn_fusion_ring, group_rep_fusion_ring, hi_fusion_ring,
       ty_fusion_ring

range_psu2k(i, j, k) = abs(i - j):2:min(i + j, 2k - i - j)

# PSU(2)_k
function psu2k_fusion_ring(k::Int)::FusionRing
    rk = div(k, 2) + 1
    mt = fill(0, rk, rk, rk)
    for a in 0:2:k, b in 0:2:k, c in 0:2:k
        c in range_psu2k(a, b, k) && (mt[div(a, 2)+1, div(b, 2)+1, div(c, 2)+1] = 1)
    end
    fusion_ring(mt,
                names = ["PSU(2)" * subscript_integer(k)],
                element_names = [denominator((i-1)//2) == 1 ? string((i-1)//2) :
                                 string(numerator((i-1)//2))*"/"*string(denominator((i-1)//2))
                                 for i in 1:rk])
end

# SU(2)_k
function su2k_fusion_ring(k::Int)::FusionRing
    rk = k + 1
    mt = fill(0, rk, rk, rk)
    for a in 0:k, b in 0:k, c in 0:k
        c in range_psu2k(a, b, k) && (mt[a+1, b+1, c+1] = 1)
    end
    fusion_ring(mt, names = ["SU(2)" * subscript_integer(k)])
end


function zn_fusion_ring(n::Int)::FusionRing
    mt = fill(0, n, n, n)
    for i in 0:n-1, j in 0:n-1
        k = mod(i + j, n)
        mt[i+1, j+1, k+1] = 1
    end
    fusion_ring(mt,
                names = ["Z_" * string(n)],
                element_names = string.(0:n-1))
end

# fusion‑ring creation from a group multiplication table

function is_cayley_table(gmt::Array{Int, 2})
    r = size(gmt, 1)
    size(gmt, 2) == r || return false
    # each row/col is a permutation of 1:r
    all(all(sort(gmt[i, :]) == 1:r for i in 1:r)) || return false
    all(all(sort(gmt[:, j]) == 1:r for j in 1:r)) || return false
    # each element appears exactly once in its own row/col diag -> inverses
    # associativity check via fusion_ring constructor later
    true
end

function fusion_ring_from_group(gmt::Array{Int, 2}; skipcheck::Bool = false)::FusionRing
    !skipcheck && is_cayley_table(gmt) || error("Provided table is not a valid Cayley table")
    r = size(gmt, 1)
    mt = fill(0, r, r, r)
    for i in 1:r, j in 1:r
        mt[i, j, gmt[i, j]] = 1
    end
    fusion_ring(mt)
end

# Overload for actual group objects (needs character data)
function fusion_ring_from_group(grp)
    throw(ErrorException("fusion_ring_from_group(grp) not yet implemented — require group algebra / character data"))
end

# Tambara–Yamagami rings

function ty_fusion_ring(G::AbstractVector)::FusionRing
    n = length(G)
    rank = n + 1
    mt = fill(0, rank, rank, rank)
    # group object fusion
    for i in 1:n, j in 1:n
        k = mod(i + j - 2, n) + 1
        mt[i, j, k] = 1
    end
    m = rank
    for i in 1:n
        mt[i, m, m] = 1; mt[m, i, m] = 1
    end
    for i in 1:n
        mt[m, m, i] = 1
    end
    fusion_ring(mt,
                names = ["TY(" * join(G, ",") * ")"],
                element_names = vcat(string.(G), ["m"]))
end

#to do:
son2_fusion_ring(n::Int) = throw(ErrorException("son2_fusion_ring requires SO(n)_2 fusion rules (TODO)"))
metaplectic_fusion_ring(m::Int) = throw(ErrorException("metaplectic_fusion_ring not yet implemented"))

group_rep_fusion_ring(grp) = throw(ErrorException("group_rep_fusion_ring needs character tables (TODO)"))
hi_fusion_ring(grp)       = throw(ErrorException("hi_fusion_ring (Haagerup–Izumi) pending implementation"))

groupname(grp) = try string(grp) catch; "Unknown Group" end


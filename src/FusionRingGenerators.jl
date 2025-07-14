# FusionRingGenerators.jl — utilities for constructing and validating fusion rings

include("GeneralFunctions.jl")   # utility helpers (bold_integer, subscript_integer, …)
using LinearAlgebra               # for `I` (uniform scaling identity)


struct FusionRing
    multiplication_table               # 3‑tensor of structure constants (Int)
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


is_constant_array(dims) = length(unique(dims)) == 1

function check_struct_const(mt)
    all(x -> x >= 0, mt)
end

function check_mt_dims(mt)
    dims = size(mt)
    length(dims) == 3 && is_constant_array(dims)
end

function check_unit(mt)
    # first simple object acts as the unit ↔ first index slice equals identity
    mt[1, :, :] == mt[:, 1, :] == I
end

function check_inverse(mt)
    r = size(mt, 1)
    sum(mt[:, :, 1]) == r            # each simple object appears exactly once in column/row for the unit
end

function check_associativity(mt::Array{Int,3})
    r = size(mt, 1)
    for a in 1:r, b in 1:r, c in 1:r, d in 1:r
        lhs = sum(mt[a, b, e] * mt[e, c, d] for e in 1:r)
        rhs = sum(mt[a, f, d] * mt[b, c, f] for f in 1:r)
        lhs == rhs || return false
    end
    true
end

function check_element_names(mt, names)
    size(mt, 1) == length(names)
end

# Fusion‑ring constructor 

function fusion_ring(mt;
    names = missing,
    texnames = missing,
    element_names = missing,
    barcode = missing,
    formal_code = missing,
    tensor_product_decompositions = missing,
    sub_fusion_rings = missing,
    frobenius_perron_dimensions = missing,
    modular_data = missing,
    characters = missing)

    !check_struct_const(mt)        && error("All structure constants should be non‑negative integers")
    !check_mt_dims(mt)             && error("multiplication_table must be a 3‑tensor with equal side lengths")
    !check_unit(mt)                && error("Unit object must be first basis element (index 1)")
    !check_inverse(mt)             && error("Some objects lack inverses (or have multiple)")
    !check_associativity(mt)       && error("Multiplication tensor is not associative")
    element_names !== missing && !check_element_names(mt, element_names) &&
        error("Length of element_names must equal rank")

    if element_names === missing
        element_names = [bold_integer(i) for i in 1:size(mt,1)]
    end

    FusionRing(
        Int.(mt),              # ensure plain integers
        names,
        texnames,
        element_names,
        barcode,
        formal_code,
        tensor_product_decompositions,
        sub_fusion_rings,
        frobenius_perron_dimensions,
        modular_data,
        characters)
end


export psu2k_fusion_ring, su2k_fusion_ring, son2_fusion_ring, metaplectic_fusion_ring,
       fusion_ring_from_group, fusion_ring_from_group, zn_fusion_ring, group_rep_fusion_ring,
       hi_fusion_ring, ty_fusion_ring

# PSU(2)_k

range_psu2k(i, j, k) = abs(i - j):2:min(i + j, 2 * k - i - j)

function psu2k_fusion_ring(k::Int)::FusionRing
    rk = div(k, 2) + 1                                 # rank
    mt = fill(0, rk, rk, rk)
    for a in 0:2:k, b in 0:2:k, c in 0:2:k
        c in range_psu2k(a, b, k) &&
            (mt[div(a, 2) + 1, div(b, 2) + 1, div(c, 2) + 1] = 1)
    end

    fusion_ring(mt,
        names = ["PSU(2)" * subscript_integer(k)],
        element_names = [to_psu2_elname((i - 1)//2) for i in 1:rk])
end

to_psu2_elname(a::Rational) = denominator(a) == 1 ? string(numerator(a)) : string(numerator(a)) * "/" * string(denominator(a))

# SU(2)_k

function su2k_fusion_ring(k::Int)::FusionRing
    rk = k + 1
    mt = fill(0, rk, rk, rk)
    for a in 0:k, b in 0:k, c in 0:k
        c in range_psu2k(a, b, k) && (mt[a + 1, b + 1, c + 1] = 1)
    end
    fusion_ring(mt, names = ["SU(2)" * subscript_integer(k)])
end

# Cyclic group

function zn_fusion_ring(n::Int)::FusionRing
    mt = fill(0, n, n, n)
    for i in 0:n-1, j in 0:n-1
        k = mod(i + j, n)
        mt[i + 1, j + 1, k + 1] = 1
    end
    fusion_ring(mt,
        names = ["Z_" * string(n)],
        element_names = [string(i) for i in 0:n-1])
end

# Fusion ring from  group multiplication table

function fusion_ring_from_group(gmt::AbstractMatrix{Int}; skipcheck::Bool = false)::FusionRing
    rank = size(gmt, 1)
    @assert size(gmt,2) == rank "Group multiplication table must be square"
    mt = fill(0, rank, rank, rank)
    for i in 1:rank, j in 1:rank
        k = gmt[i,j]
        mt[i, j, k] = 1
    end
    fusion_ring(mt)   # will validate associativity etc.
end

# Tambara‑Yamagami fusion rings for an abelian group G (labels in Vector)

function ty_fusion_ring(G::AbstractVector)::FusionRing
    n = length(G)
    rank = n + 1
    mt = fill(0, rank, rank, rank)

    # group object fusion
    for i in 1:n, j in 1:n
        k = mod(i + j - 2, n) + 1   # additive group law on indices
        mt[i, j, k] = 1
    end

    m = rank                       # label index for the extra "m" object

    # m ⊗ g = m = g ⊗ m
    for i in 1:n
        mt[i, m, m] = 1
        mt[m, i, m] = 1
    end

    # m ⊗ m = ⊕ g
    for i in 1:n
        mt[m, m, i] = 1
    end

    fusion_ring(mt,
        names = ["TY(" * join(G, ",") * ")"],
        element_names = vcat(string.(G), ["m"]))
end

# Placeholders 

function son2_fusion_ring(n::Int)::FusionRing
    throw(ErrorException("son2_fusion_ring not yet implemented — requires SO(n)_2 fusion rules"))
end

function metaplectic_fusion_ring(m::Int)::FusionRing
    throw(ErrorException("metaplectic_fusion_ring not yet implemented"))
end

function fusion_ring_from_group(grp)  # overloaded for actual group object
    throw(ErrorException("fusion_ring_from_group(grp) requires group algebra or character data"))
end

function group_rep_fusion_ring(grp)
    throw(ErrorException("group_rep_fusion_ring not yet implemented (needs character tables)"))
end

function hi_fusion_ring(grp)
    throw(ErrorException("hi_fusion_ring (Haagerup‑Izumi) needs concrete construction details"))
end


groupname(grp) = try string(grp) catch; "Unknown Group" end

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


export psu2k_fusion_ring, su2k_fusion_ring, son2_fusion_ring, metaplectic_fusion_ring,
       fusion_ring_from_group, zn_fusion_ring, group_rep_fusion_ring, hi_fusion_ring,
       ty_fusion_ring

range_psu2k(i, j, k) = abs(i - j):2:min(i + j, 2k - i - j)

# TODO: add missing information
# TODO: code for labels is a bit too dense
# PSU(2)_k
function psu2k_fusion_ring(k::Int)::FusionRing
    rk = div(k, 2) + 1
    mt = fill(0, rk, rk, rk)
    for a in 0:2:k, b in 0:2:k, c in 0:2:k
        c in range_psu2k(a, b, k) && (mt[div(a, 2)+1, div(b, 2)+1, div(c, 2)+1] = 1)
    end

    elnames = 
        [
            denominator((i-1)//2) == 1 ? 
            string((i-1)//2) :
            string(numerator((i-1)//2))*"/"*string(denominator((i-1)//2))
            for i in 1:rk
        ]
    
    fusion_ring(
        mt,
        names = ["PSU(2)" * subscript_integer(k)],
        labels = elnames
    )
end

# TODO: add missing information
# SU(2)_k
function su2k_fusion_ring(k::Int)::FusionRing
    rk = k + 1
    mt = fill(0, rk, rk, rk)
    for a in 0:k, b in 0:k, c in 0:k
        c in range_psu2k(a, b, k) && (mt[a+1, b+1, c+1] = 1)
    end
    fusion_ring(
        mt, 
        names = ["SU(2)" * subscript_integer(k)],
        labels = string.(0:k)
    )
end


# TODO: add missing information
function zn_fusion_ring(n::Int)::FusionRing
    mt = fill(0, n, n, n)
    for i in 0:n-1, j in 0:n-1
        k = mod(i + j, n)
        mt[i+1, j+1, k+1] = 1
    end
    fusion_ring(
        mt,
        names = ["Z_" * string(n)],
        labels = string.(0:n-1)
    )
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

# TODO: add missing information
function fusion_ring_from_group(gmt::Array{Int, 2}; skipcheck::Bool = false)::FusionRing
    !skipcheck && is_cayley_table(gmt) || error("Provided table is not a valid Cayley table")
    r = size(gmt, 1)
    mt = fill(0, r, r, r)
    for i in 1:r, j in 1:r
        mt[i, j, gmt[i, j]] = 1
    end
    fusion_ring(mt, skip_check = skipcheck)
end

# TODO: Overload for actual group objects (needs character data)
function fusion_ring_from_group(grp)
    throw(ErrorException("fusion_ring_from_group(grp) not yet implemented — require group algebra / character data"))
end

# TODO: add missing information
# TODO: this doesn't look correct. It should work for any group,
# not just cyclic ones
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
    fusion_ring(
        mt,
        names = ["TY(" * join(G, ",") * ")"],
        labels = vcat(string.(G), ["m"])
    )
end

# TODO: implement 
son2_fusion_ring(n::Int) = throw(ErrorException("son2_fusion_ring requires SO(n)_2 fusion rules (TODO)"))
# TODO: implement 
metaplectic_fusion_ring(m::Int) = throw(ErrorException("metaplectic_fusion_ring not yet implemented"))

# TODO: implement 
group_rep_fusion_ring(grp) = throw(ErrorException("group_rep_fusion_ring needs character tables (TODO)"))
# TODO: implement 
hi_fusion_ring(grp)       = throw(ErrorException("hi_fusion_ring (Haagerup–Izumi) pending implementation"))




# Haagerup–Izumi (HI) and Tambara–Yamagami (TY) fusion rings
# - Input `tab` is  n×n group multiplication table on {1,…,n} with identity = 1.
# - Output `mt` is a rank×rank×rank multiplication tensor with structure constants
#     mt[i,j,k] = multiplicity of simple k in i ⊗ j.
# must already have:
#   - struct FusionRing with fields `multiplication_table`, `names`, `labels` (etc.)
#   - `fusion_ring(mt; names=..., labels=...)` constructor
#
# This file provides:
#   FusionRingHI(tab; names=...)
#   FusionRingTY(tab; names=...)
=

export FusionRingHI, FusionRingTY



"""
    _is_group_table(tab) -> Bool

Very explicit check that `tab` is a group multiplication table on {1..n}
with identity element 1.

Checks:
- tab is n×n Int
- entries are in 1..n
- 1 acts as identity: tab[1,i]=i and tab[i,1]=i
- each row and column is a permutation of 1..n
- associativity: tab[ tab[i,j], k ] == tab[ i, tab[j,k] ]
"""
function _is_group_table(tab::AbstractMatrix{<:Integer})::Bool
    n = size(tab, 1)
    size(tab, 2) == n || return false
    n ≥ 1 || return false

    # Entries in 1..n
    @inbounds for i in 1:n, j in 1:n
        x = tab[i, j]
        (1 <= x <= n) || return false
    end

    # Identity is 1
    @inbounds for i in 1:n
        tab[1, i] == i || return false
        tab[i, 1] == i || return false
    end

    # Latin square: each row/col is a permutation of 1..n
    seen = falses(n)
    @inbounds for i in 1:n
        fill!(seen, false)
        for j in 1:n
            seen[tab[i, j]] = true
        end
        all(seen) || return false

        fill!(seen, false)
        for j in 1:n
            seen[tab[j, i]] = true
        end
        all(seen) || return false
    end

    # Associativity
    @inbounds for i in 1:n, j in 1:n, k in 1:n
        tab[tab[i, j], k] == tab[i, tab[j, k]] || return false
    end

    return true
end

"""
    _inverse_vector(tab) -> inv

Return inv[1..n] where inv[a] is the (unique) inverse of `a` in the group-table `tab`,
i.e. tab[a, inv[a]] == 1.
"""
function _inverse_vector(tab::AbstractMatrix{<:Integer})::Vector{Int}
    n = size(tab, 1)
    inv = zeros(Int, n)
    @inbounds for a in 1:n
        found = 0
        for b in 1:n
            if tab[a, b] == 1
                found = b
                break
            end
        end
        found == 0 && error("Group table has no inverse for element $a (no b with tab[a,b]=1).")
        inv[a] = found
    end
    return inv
end



"""
    _mats_to_mt(mats) -> mt

Given mats[a] = N_a (rank×rank), return mt[a,b,c] = (N_a)[b,c].
"""
function _mats_to_mt(mats::Vector{<:AbstractMatrix{<:Integer}})::Array{Int,3}
    r = length(mats)
    r ≥ 1 || error("_mats_to_mt: empty mats")
    mt = zeros(Int, r, r, r)
    @inbounds for a in 1:r
        A = mats[a]
        size(A,1) == r && size(A,2) == r || error("_mats_to_mt: mat $a has wrong size $(size(A)) (expected $r×$r)")
        mt[a, :, :] .= A
    end
    return mt
end


"""
    FusionRingTY(tab; names=String[]) -> FusionRing

Build the Tambara–Yamagami fusion ring for a group with multiplication table `tab`.
Rank is n+1 (group elements + one extra object).
"""
function FusionRingTY(tab::AbstractMatrix{<:Integer}; names::Vector{String}=String[])
    _is_group_table(tab) || throw(ArgumentError("FusionRingTY: tab must be a group multiplication table (identity=1, associative, latin square)."))
    n = size(tab, 1)
    r = n + 1

    mats = Matrix{Int}[]

    # For each simple object i=1..r, build its fusion matrix N_i.
    # This mirrors the Mathematica Which[...] table.
    @inbounds for i in 1:r
        Ni = zeros(Int, r, r)
        for j in 1:r
            if i <= n && j <= n
                k = tab[i, j]
                Ni[j, k] += 1
            elseif i <= n && j > n
                # group element ⊗ m = m
                Ni[j, r] += 1
            elseif i > n && j <= n
                # m ⊗ group element = m
                Ni[j, r] += 1
            else
                # m ⊗ m = sum_{g in G} g
                for k in 1:n
                    Ni[j, k] += 1
                end
            end
        end
        push!(mats, Ni)
    end

    mt = _mats_to_mt(mats)

    # labels: 1..n are group elements, last is "m"
    labels = [string(i) for i in 1:n]
    push!(labels, "m")

    default_names = isempty(names) ? String[] : names
    return fusion_ring(mt; names=default_names, labels=labels)
end

# Haagerup–Izumi

"""
    FusionRingHI(tab; names=String[]) -> FusionRing

Build the Haagerup–Izumi fusion ring from a *symmetric* group multiplication table `tab`.

Rank is 2n. Objects are:
- 1..n   : group elements
- n+1..2n: "rho*g" sector (s X_g), indexed by g=1..n as n+g.


"""
function FusionRingHI(tab::AbstractMatrix{<:Integer}; names::Vector{String}=String[])
    _is_group_table(tab) || throw(ArgumentError("FusionRingHI: tab must be a group multiplication table (identity=1, associative, latin square)."))
    issymmetric(tab) || throw(ArgumentError("FusionRingHI: multiplication table must be symmetric."))

    n = size(tab, 1)
    r = 2n
    inv = _inverse_vector(tab)

    mats = Matrix{Int}[]

    # For i in 1..2n build N_i as in  Mathematica Which cases.
    @inbounds for i in 1:r
        Ni = zeros(Int, r, r)
        for j in 1:r
            if i <= n && j <= n
                # k == tab[[i,j]]
                k = tab[i, j]
                Ni[j, k] += 1

            elseif i <= n && j > n
                # k == n + tab[[i, j-n]]
                k = n + tab[i, j - n]
                Ni[j, k] += 1

            elseif i > n && j <= n
                # k == n + tab[[ inv[[j]], i-n ]]
                k = n + tab[inv[j], i - n]
                Ni[j, k] += 1

            else
                # i>n && j>n:
                # If[ k == tab[[ i-n, inv[[j-n]] ]] || k > n, 1, 0 ]
                # => all "rho-sector" (k>n) appear with multiplicity 1,
                #    plus exactly one group element tab[i-n, inv[j-n]].
                k0 = tab[i - n, inv[j - n]]
                Ni[j, k0] += 1
                for k in (n+1):r
                    Ni[j, k] += 1
                end
            end
        end
        push!(mats, Ni)
    end

    mt = _mats_to_mt(mats)

    # Labels: group elements then rho-sector
    labels = [string(i) for i in 1:n]
    append!(labels, ["ρ_$i" for i in 1:n])  
    default_names = isempty(names) ? String[] : names
    return fusion_ring(mt; names=default_names, labels=labels)
end

groupname(grp) = try string(grp) catch; "Unknown Group" end

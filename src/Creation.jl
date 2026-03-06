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
group_rep_fusion_ring(grp) = throw(ErrorException("group_rep_fusion_ring needs character tables (TODO)"))
# TODO: implement 
hi_fusion_ring(grp)       = throw(ErrorException("hi_fusion_ring (Haagerup–Izumi) pending implementation"))




# TODO 
# 1. Could use unicode to make everything more readable
# 2. Could use dictionaries rather than elseif statements 

# Anyonica rulesodd[m_]  (metaplectic / SO(m)_2, m odd)
# 
# - Anyonica builds matX as a Table of vectors (length=rank),
#   then does Transpose /@ at the end.
# - i follow that exactly: build "row_table" (rank×rank),
#   then transpose before converting to mt.

# basis vector e_i in ℤ^rank
@inline function _e(i::Int, rank::Int)::Vector{Int}
    v = zeros(Int, rank)
    v[i] = 1
    return v
end

# convert a list of fusion matrices mats[a][b,c] into mt[a,b,c]
function _mats_to_mt(mats::Vector{Matrix{Int}})::Array{Int,3}
    r = length(mats)
    mt = zeros(Int, r, r, r)
    @inbounds for a in 1:r
        mt[a, :, :] .= mats[a]
    end
    return mt
end

function _son2_rules_odd(m::Integer)::Array{Int,3}
    isodd(m) || throw(ArgumentError("_son2_rules_odd expects odd m, got m=$m"))
    m ≥ 5    || throw(ArgumentError("_son2_rules_odd expects m≥5 (odd), got m=$m"))

    r    = (m - 1) ÷ 2
    rank = (m + 7) ÷ 2 

    # convenience
    ar(i) = _e(i, rank)

    # mat1 = IdentityMatrix[rank]
    mat1 = Matrix{Int}(I, rank, rank)

    # matZ = Table[ Which[...], {i,rank} ]
    matZ = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == 1
            ar(2)
        elseif i == 2
            ar(1)
        elseif 3 <= i <= 4
            ar(3 + mod(i, 2))   # i=3 -> 4, i=4 -> 3
        else
            ar(i)
        end
        matZ[i, :] .= v
    end

    # matXe1 = Table[ Which[...], {i,rank} ]
    matXe1 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == 1
            ar(3)
        elseif i == 2
            ar(4)
        elseif i == 3
            # ar[1] + Sum[ar[j], {j, 5, rank}]
            v0 = zeros(Int, rank)
            v0[1] = 1
            for j in 5:rank
                v0[j] += 1
            end
            v0
        elseif i == 4
            # ar[2] + Sum[ar[j], {j, 5, rank}]
            v0 = zeros(Int, rank)
            v0[2] = 1
            for j in 5:rank
                v0[j] += 1
            end
            v0
        else
            ar(3) .+ ar(4)
        end
        matXe1[i, :] .= v
    end

    # matXe2 = Table[ Which[...], {i,rank} ]
    matXe2 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == 1
            ar(4)
        elseif i == 2
            ar(3)
        elseif i == 3
            # ar[2] + Sum[ar[j], {j, 5, rank}]
            v0 = zeros(Int, rank)
            v0[2] = 1
            for j in 5:rank
                v0[j] += 1
            end
            v0
        elseif i == 4
            # ar[1] + Sum[ar[j], {j, 5, rank}]
            v0 = zeros(Int, rank)
            v0[1] = 1
            for j in 5:rank
                v0[j] += 1
            end
            v0
        else
            ar(3) .+ ar(4)
        end
        matXe2[i, :] .= v
    end

    # matY[j_] := Table[ Which[...], {i,rank} ]
    function matY(j::Int)::Matrix{Int}
        M = zeros(Int, rank, rank)
        @inbounds for i in 1:rank
            v = if i == 1
                ar(j + 4)
            elseif i == 2
                ar(j + 4)
            elseif i == 3
                ar(3) .+ ar(4)
            elseif i == 4
                ar(3) .+ ar(4)
            else
                ii = i - 4
                if ii == j
                    # ar[1] + ar[2] + ar[ Min[2j, m-2j] + 4 ]
                    t = min(2*j, m - 2*j) + 4
                    ar(1) .+ ar(2) .+ ar(t)
                else
                    # ar[ Abs[ii-j] + 4 ] + ar[ Min[ii+j, m-ii-j+4] + 4 ]
                    t1 = abs(ii - j) + 4
                    t2 = min(ii + j, m - i - j + 4) + 4
                    ar(t1) .+ ar(t2)
                end
            end
            M[i, :] .= v
        end
        return M
    end

    #   Transpose /@ Join[{mat1, matZ, matXe1, matXe2}, matY /@ Range[r]]
    mats = Matrix{Int}[]
    push!(mats, transpose(mat1))
    push!(mats, transpose(matZ))
    push!(mats, transpose(matXe1))
    push!(mats, transpose(matXe2))
    for j in 1:r
        push!(mats, transpose(matY(j)))
    end

    # Convert fusion matrices to multiplication table
    return _mats_to_mt(mats)
end



#  rulesdiv2[p_] and rulesdiv4[p_]  (SO(m)_2 even cases)
#
# Used by son2_fusion_ring:
#   if m % 4 == 0  => rulesdiv4(m/2)
#   elseif m % 2 == 0 => rulesdiv2(m/2)
#
# 
#   build fusion matrices mats[a] as Integer matrices,
#   apply transpose to match Anyonica's Transpose /@,
#   convert mats -> multiplication table mt[a,b,c].



function _mats_to_mt(mats::Vector{Matrix{Int}})::Array{Int,3}
    r = length(mats)
    mt = zeros(Int, r, r, r)
    @inbounds for a in 1:r
        mt[a, :, :] .= mats[a]
    end
    return mt
end

# index convention (matches  Evaluate[...] = IdentityMatrix[rank]):
# 1: Id
# 2: Θ
# 3: Φ1
# 4: Φ2
# 5: σ1
# 6: σ2
# 7: τ1
# 8: τ2
# 9..rank:  Φ[j] with j = 1..(rank-8)
@inline _Id() = 1
@inline _Th() = 2
@inline _Phi1() = 3
@inline _Phi2() = 4
@inline _Sig1() = 5
@inline _Sig2() = 6
@inline _Tau1() = 7
@inline _Tau2() = 8
@inline _Phi(j::Int) = 8 + j

# rulesdiv2[p_]
function _son2_rules_div2(p::Integer)::Array{Int,3}
    p ≥ 1 || throw(ArgumentError("_son2_rules_div2 expects p≥1, got p=$p"))
    rank = p + 7
    maxphi = rank - 8  # = p-1

    ar(i) = _e(i, rank)

    # sums: sumEvenΛs = Σ_{i=2,4,...,p-1} Φ[i], sumOddΛs = Σ_{i=1,3,...,p-1} Φ[i]
    sumEven = zeros(Int, rank)
    sumOdd  = zeros(Int, rank)
    for i in 1:maxphi
        (isodd(i) ? (sumOdd[_Phi(i)] += 1) : (sumEven[_Phi(i)] += 1))
    end

    mats = Matrix{Int}[]

    # matId = IdentityMatrix[rank]
    matId = Matrix{Int}(I, rank, rank)
    push!(mats, transpose(matId))

    # matΘ
    matTh = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Th())
        elseif i == _Th()
            ar(_Id())
        elseif i == _Phi1()
            ar(_Phi2())
        elseif i == _Phi2()
            ar(_Phi1())
        elseif i == _Sig1()
            ar(_Tau1())
        elseif i == _Sig2()
            ar(_Tau2())
        elseif i == _Tau1()
            ar(_Sig1())
        elseif i == _Tau2()
            ar(_Sig2())
        else
            ar(i) # Φ-lambdas fixed
        end
        matTh[i, :] .= v
    end
    push!(mats, transpose(matTh))

    # matΦ1
    matPhi1 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Phi1())
        elseif i == _Th()
            ar(_Phi2())
        elseif i == _Phi1()
            ar(_Th())
        elseif i == _Phi2()
            ar(_Id())
        elseif i == _Sig1()
            ar(_Sig2())
        elseif i == _Sig2()
            ar(_Tau1())
        elseif i == _Tau1()
            ar(_Tau2())
        elseif i == _Tau2()
            ar(_Sig1())
        else
            # Φ[p - (i-8)]
            j = i - 8
            ar(_Phi(p - j))
        end
        matPhi1[i, :] .= v
    end
    push!(mats, transpose(matPhi1))

    # matΦ2
    matPhi2 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Phi2())
        elseif i == _Th()
            ar(_Phi1())
        elseif i == _Phi1()
            ar(_Id())
        elseif i == _Phi2()
            ar(_Th())
        elseif i == _Sig1()
            ar(_Tau2())
        elseif i == _Sig2()
            ar(_Sig1())
        elseif i == _Tau1()
            ar(_Sig2())
        elseif i == _Tau2()
            ar(_Tau1())
        else
            j = i - 8
            ar(_Phi(p - j))
        end
        matPhi2[i, :] .= v
    end
    push!(mats, transpose(matPhi2))

    # matσ1
    matSig1 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Sig1())
        elseif i == _Th()
            ar(_Tau1())
        elseif i == _Phi1()
            ar(_Sig2())
        elseif i == _Phi2()
            ar(_Tau2())
        elseif i == _Sig1()
            ar(_Phi2()) .+ sumOdd
        elseif i == _Sig2()
            ar(_Id()) .+ sumEven
        elseif i == _Tau1()
            ar(_Phi1()) .+ sumOdd
        elseif i == _Tau2()
            ar(_Th()) .+ sumEven
        else
            # If[OddQ[i], σ2+τ2, σ1+τ1]
            isodd(i) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
        end
        matSig1[i, :] .= v
    end
    push!(mats, transpose(matSig1))

    # matσ2
    matSig2 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Sig2())
        elseif i == _Th()
            ar(_Tau2())
        elseif i == _Phi1()
            ar(_Tau1())
        elseif i == _Phi2()
            ar(_Sig1())
        elseif i == _Sig1()
            ar(_Id()) .+ sumEven
        elseif i == _Sig2()
            ar(_Phi1()) .+ sumOdd
        elseif i == _Tau1()
            ar(_Th()) .+ sumEven
        elseif i == _Tau2()
            ar(_Phi2()) .+ sumOdd
        else
            # If[EvenQ[i], σ2+τ2, σ1+τ1]
            iseven(i) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
        end
        matSig2[i, :] .= v
    end
    push!(mats, transpose(matSig2))

    # matτ1
    matTau1 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Tau1())
        elseif i == _Th()
            ar(_Sig1())
        elseif i == _Phi1()
            ar(_Tau2())
        elseif i == _Phi2()
            ar(_Sig2())
        elseif i == _Sig1()
            ar(_Phi1()) .+ sumOdd
        elseif i == _Sig2()
            ar(_Th()) .+ sumEven
        elseif i == _Tau1()
            ar(_Phi2()) .+ sumOdd
        elseif i == _Tau2()
            ar(_Id()) .+ sumEven
        else
            isodd(i) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
        end
        matTau1[i, :] .= v
    end
    push!(mats, transpose(matTau1))

    # matτ2
    matTau2 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Tau2())
        elseif i == _Th()
            ar(_Sig2())
        elseif i == _Phi1()
            ar(_Sig1())
        elseif i == _Phi2()
            ar(_Tau1())
        elseif i == _Sig1()
            ar(_Th()) .+ sumEven
        elseif i == _Sig2()
            ar(_Phi2()) .+ sumOdd
        elseif i == _Tau1()
            ar(_Id()) .+ sumEven
        elseif i == _Tau2()
            ar(_Phi1()) .+ sumOdd
        else
            iseven(i) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
        end
        matTau2[i, :] .= v
    end
    push!(mats, transpose(matTau2))

    # matΦ[j] for j = 1..(rank-8)
    function matPhi(j::Int)::Matrix{Int}
        M = zeros(Int, rank, rank)
        @inbounds for i in 1:rank
            v = if i == _Id()
                ar(_Phi(j))
            elseif i == _Th()
                ar(_Phi(j))
            elseif i == _Phi1()
                ar(_Phi(p - j))
            elseif i == _Phi2()
                ar(_Phi(p - j))
            elseif i == _Sig1()
                isodd(j) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
            elseif i == _Sig2()
                iseven(j) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
            elseif i == _Tau1()
                isodd(j) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
            elseif i == _Tau2()
                iseven(j) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
            else
                ii = i - 8
                if ii == j && (2*j < p)
                    ar(_Id()) .+ ar(_Th()) .+ ar(_Phi(2*j))
                elseif ii == j && (2*j > p)
                    ar(_Id()) .+ ar(_Th()) .+ ar(_Phi(2*(p - j)))
                elseif ii + j < p
                    ar(_Phi(abs(ii - j))) .+ ar(_Phi(ii + j))
                elseif ii + j > p
                    ar(_Phi(abs(ii - j))) .+ ar(_Phi(2*p - ii - j))
                else
                    # ii == p - j
                    ar(_Phi1()) .+ ar(_Phi2()) .+ ar(_Phi(abs(p - 2*ii)))
                end
            end
            M[i, :] .= v
        end
        return M
    end

    for j in 1:maxphi
        push!(mats, transpose(matPhi(j)))
    end

    return _mats_to_mt(mats)
end


# rulesdiv4[p_]
function _son2_rules_div4(p::Integer)::Array{Int,3}
    p ≥ 1 || throw(ArgumentError("_son2_rules_div4 expects p≥1, got p=$p"))
    rank = p + 7
    maxphi = rank - 8  # = p-1

    ar(i) = _e(i, rank)

    sumEven = zeros(Int, rank)
    sumOdd  = zeros(Int, rank)
    for i in 1:maxphi
        (isodd(i) ? (sumOdd[_Phi(i)] += 1) : (sumEven[_Phi(i)] += 1))
    end

    mats = Matrix{Int}[]

    matId = Matrix{Int}(I, rank, rank)
    push!(mats, transpose(matId))

    # matΘ
    matTh = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Th())
        elseif i == _Th()
            ar(_Id())
        elseif i == _Phi1()
            ar(_Phi2())
        elseif i == _Phi2()
            ar(_Phi1())
        elseif i == _Sig1()
            ar(_Tau1())
        elseif i == _Sig2()
            ar(_Tau2())
        elseif i == _Tau1()
            ar(_Sig1())
        elseif i == _Tau2()
            ar(_Sig2())
        else
            ar(i)
        end
        matTh[i, :] .= v
    end
    push!(mats, transpose(matTh))

    # matΦ1
    matPhi1 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Phi1())
        elseif i == _Th()
            ar(_Phi2())
        elseif i == _Phi1()
            ar(_Id())
        elseif i == _Phi2()
            ar(_Th())
        elseif i == _Sig1()
            ar(_Tau1())
        elseif i == _Sig2()
            ar(_Sig2())
        elseif i == _Tau1()
            ar(_Sig1())
        elseif i == _Tau2()
            ar(_Tau2())
        else
            j = i - 8
            ar(_Phi(p - j))
        end
        matPhi1[i, :] .= v
    end
    push!(mats, transpose(matPhi1))

    # matΦ2
    matPhi2 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Phi2())
        elseif i == _Th()
            ar(_Phi1())
        elseif i == _Phi1()
            ar(_Th())
        elseif i == _Phi2()
            ar(_Id())
        elseif i == _Sig1()
            ar(_Sig1())
        elseif i == _Sig2()
            ar(_Tau2())
        elseif i == _Tau1()
            ar(_Tau1())
        elseif i == _Tau2()
            ar(_Sig2())
        else
            j = i - 8
            ar(_Phi(p - j))
        end
        matPhi2[i, :] .= v
    end
    push!(mats, transpose(matPhi2))

    # matσ1
    matSig1 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Sig1())
        elseif i == _Th()
            ar(_Tau1())
        elseif i == _Phi1()
            ar(_Tau1())
        elseif i == _Phi2()
            ar(_Sig1())
        elseif i == _Sig1()
            ar(_Id()) .+ ar(_Phi2()) .+ sumEven
        elseif i == _Sig2()
            sumOdd
        elseif i == _Tau1()
            ar(_Th()) .+ ar(_Phi1()) .+ sumEven
        elseif i == _Tau2()
            sumOdd
        else
            isodd(i) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
        end
        matSig1[i, :] .= v
    end
    push!(mats, transpose(matSig1))

    # matσ2
    matSig2 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Sig2())
        elseif i == _Th()
            ar(_Tau2())
        elseif i == _Phi1()
            ar(_Sig2())
        elseif i == _Phi2()
            ar(_Tau2())
        elseif i == _Sig1()
            sumOdd
        elseif i == _Sig2()
            ar(_Id()) .+ ar(_Phi1()) .+ sumEven
        elseif i == _Tau1()
            sumOdd
        elseif i == _Tau2()
            ar(_Th()) .+ ar(_Phi2()) .+ sumEven
        else
            iseven(i) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
        end
        matSig2[i, :] .= v
    end
    push!(mats, transpose(matSig2))

    # matτ1
    matTau1 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Tau1())
        elseif i == _Th()
            ar(_Sig1())
        elseif i == _Phi1()
            ar(_Sig1())
        elseif i == _Phi2()
            ar(_Tau1())
        elseif i == _Sig1()
            ar(_Th()) .+ ar(_Phi1()) .+ sumEven
        elseif i == _Sig2()
            sumOdd
        elseif i == _Tau1()
            ar(_Id()) .+ ar(_Phi2()) .+ sumEven
        elseif i == _Tau2()
            sumOdd
        else
            isodd(i) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
        end
        matTau1[i, :] .= v
    end
    push!(mats, transpose(matTau1))

    # matτ2
    matTau2 = zeros(Int, rank, rank)
    @inbounds for i in 1:rank
        v = if i == _Id()
            ar(_Tau2())
        elseif i == _Th()
            ar(_Sig2())
        elseif i == _Phi1()
            ar(_Tau2())
        elseif i == _Phi2()
            ar(_Sig2())
        elseif i == _Sig1()
            sumOdd
        elseif i == _Sig2()
            ar(_Th()) .+ ar(_Phi2()) .+ sumEven
        elseif i == _Tau1()
            sumOdd
        elseif i == _Tau2()
            ar(_Id()) .+ ar(_Phi1()) .+ sumEven
        else
            iseven(i) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
        end
        matTau2[i, :] .= v
    end
    push!(mats, transpose(matTau2))

    # matΦ[j]
    function matPhi(j::Int)::Matrix{Int}
        M = zeros(Int, rank, rank)
        @inbounds for i in 1:rank
            v = if i == _Id()
                ar(_Phi(j))
            elseif i == _Th()
                ar(_Phi(j))
            elseif i == _Phi1()
                ar(_Phi(p - j))
            elseif i == _Phi2()
                ar(_Phi(p - j))
            elseif i == _Sig1()
                isodd(j) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
            elseif i == _Sig2()
                iseven(j) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
            elseif i == _Tau1()
                isodd(j) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
            elseif i == _Tau2()
                iseven(j) ? (ar(_Sig2()) .+ ar(_Tau2())) : (ar(_Sig1()) .+ ar(_Tau1()))
            else
                ii = i - 8
                if ii == j && (2*j < p)
                    ar(_Id()) .+ ar(_Th()) .+ ar(_Phi(2*j))
                elseif ii == j && (2*j == p)
                    ar(_Id()) .+ ar(_Th()) .+ ar(_Phi1()) .+ ar(_Phi2())
                elseif ii == j && (2*j > p)
                    ar(_Id()) .+ ar(_Th()) .+ ar(_Phi(2*(p - j)))
                elseif ii + j < p
                    ar(_Phi(abs(ii - j))) .+ ar(_Phi(ii + j))
                elseif ii + j > p
                    ar(_Phi(abs(ii - j))) .+ ar(_Phi(2*p - ii - j))
                else
                    # ii == p - j
                    ar(_Phi1()) .+ ar(_Phi2()) .+ ar(_Phi(abs(p - 2*ii)))
                end
            end
            M[i, :] .= v
        end
        return M
    end

    for j in 1:maxphi
        push!(mats, transpose(matPhi(j)))
    end

    return _mats_to_mt(mats)
end




# SO(m)_2 / Metaplectic(m) t)

# Uses:
#   _son2_rules_odd(m)      # for odd m
#   _son2_rules_div2(p)     # for m ≡ 2 (mod 4), with p = m÷2
#   _son2_rules_div4(p)     # for m ≡ 0 (mod 4), with p = m÷2

export son2_fusion_ring


# odd m: rank = (m+7)/2, elements are [1, Z, X_e1, X_e2, Y_1, ..., Y_r], r=(m-1)/2
function _son2_labels_odd(m::Int)::Vector{String}
    r = (m - 1) ÷ 2
    labels = String["1", "Z", "Xₑ₁", "Xₑ₂"]
    for j in 1:r
        push!(labels, "Y_$j")
    end
    return labels
end

# even m: rank = p+7 with p=m/2, elements are [Id, Θ, Φ1, Φ2, σ1, σ2, τ1, τ2, Φ_1..Φ_{p-1}]
function _son2_labels_even(p::Int)::Vector{String}
    labels = String[
        "1", "Θ", "Φ₁", "Φ₂", "σ₁", "σ₂", "τ₁", "τ₂"
    ]
    for j in 1:(p - 1)
        push!(labels, "Φ_$j")
    end
    return labels
end


"""
    son2_fusion_ring(m::Int) -> FusionRing

Return fusion ring SO(N)_2 (metaplectic) .
"""

#- odd `N`: uses `_son2_rules_odd(m)`
#- even `N ≡ 0 (mod 4)`: uses `_son2_rules_div4(N÷2)`
#- even `N ≡ 2 (mod 4)`: uses `_son2_rules_div2(N÷2)`

function son2_fusion_ring(m::Int)::FusionRing
    m ≥ 4 || throw(ArgumentError("son2_fusion_ring(m): requires integer m ≥ 4, got m=$m"))

    mt::Array{Int,3}
    labels::Vector{String}

    if isodd(m)
        mt = _son2_rules_odd(m)
        labels = _son2_labels_odd(m)
    else
        p = m ÷ 2
        if m % 4 == 0
            mt = _son2_rules_div4(p)
        else
            mt = _son2_rules_div2(p)
        end
        labels = _son2_labels_even(p)
    end

    #  label count must match rank
    size(mt, 1) == length(labels) || error("son2_fusion_ring: label length mismatch with mt rank")

    R = fusion_ring(
        mt;
        names  = ["SO($m)_2", "Metaplectic($m)"],
        labels = labels,
    )

    return replace_by_known(R)
end

export metaplectic_fusion_ring

metaplectic_fusion_ring( n::Int )::FusionRing = son2_fusion_ring(n)


groupname(grp) = try string(grp) catch; "Unknown Group" end


#Haagerup–Izumi (HI) and Tambara–Yamagami (TY) fusion rings
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
export FusionRingHI, FusionRingTY


#Added: from izumi
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


#Added: from izumi
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


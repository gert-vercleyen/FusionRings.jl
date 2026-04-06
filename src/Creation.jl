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
        if c ∈ range_psu2k(a, b, k) 
            mt[div(a, 2)+1, div(b, 2)+1, div(c, 2)+1] = 1
        else
            continue
        end
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
        if c ∈ range_psu2k(a, b, k) 
            mt[a+1, b+1, c+1] = 1
        else 
            continue
        end
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
        names = ["ℤ" * subscript_integer(n)],
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



# TODO: implement 
group_rep_fusion_ring(grp) = throw(ErrorException("group_rep_fusion_ring needs character tables (TODO)"))
# TODO: implement 
hi_fusion_ring(grp)        = throw(ErrorException("hi_fusion_ring (Haagerup–Izumi) pending implementation"))




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
"""
    _mats_to_mt(mats) -> mt

Given mats[a] = N_a (rank×rank), return mt[a,b,c] = (N_a)[b,c].
"""
function _mats_to_mt(mats::Vector{<:AbstractMatrix{<:Integer}})::Array{Int,3}
    r = length(mats)
    r ≥ 1 || error("_mats_to_mt: empty list of matrices")
    mt = zeros(Int, r, r, r)
    @inbounds for a in 1:r
        A = mats[a]
        size(A,1) == r && size(A,2) == r || error("_mats_to_mt: mat $a has wrong size $(size(A)) (expected $r×$r)")
        mt[a, :, :] .= A
    end
    return mt
end

function _son2_rules_odd(m::Integer)::Array{Int,3}
    isodd(m) || throw(ArgumentError("_son2_rules_odd expects odd N, got N=$m"))
    m ≥ 5    || throw(ArgumentError("_son2_rules_odd expects N≥5 (odd), got N=$m"))

    r    = (m - 1) ÷ 2
    rank = (m + 7) ÷ 2 

    # convenience
    ar(i) = _e(i, rank)

    # mat1 = IdentityMatrix[rank]
    mat1 = Matrix{Int}(0, rank, rank)
    for i in 1:r 
        mat1[i,i] = 1
    end

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

function son2_fusion_ring(N::Int)::FusionRing
    N ≥ 4 || throw(ArgumentError("son2_fusion_ring(N): requires integer N ≥ 4, got N=$N"))

    if isodd(N)
        mt     = _son2_rules_odd(N)
        labels = _son2_labels_odd(N)
    else
        p = N ÷ 2
        if N % 4 == 0
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
        names  = ["SO($N)_2", "Metaplectic($N)"],
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


#Added: from Iazumi  
export HI_fusion_ring
"""
    FusionRingHI(tab; names=String[]) -> FusionRing

Build the Haagerup–Izumi fusion ring from a *symmetric* group multiplication table `tab`.

Rank is 2n. Objects are:
- 1..n   : group elements
- n+1..2n: "rho*g" sector (s X_g), indexed by g=1..n as n+g.


"""
function HI_fusion_ring(tab::AbstractMatrix{<:Integer}; names::Vector{String}=String[])
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
#TODO: this is a bit redundant with _is_group_table, but we need the inverse vector for HI_fusion_ring anyway, so might as well compute it here.
#TODO: could also check that inv[inv[a]] == a, and that tab[a, inv[a]] == 1, to be extra sure.
#TODO: could also check that the group is abelian, since HI_fusion_ring requires a symmetric table.
#TODO: For dr. vercleyen: are we referring to TY categories in the standard sense, which are built from a finite abelian group plus extra data?
#TODO:  If so, should we  want to enforce commutativity of the table `tab` in `_is_group_table` or in `TY_fusion_ring`?
export TY_fusion_ring
"""
    FusionRingTY(tab; names=String[]) -> FusionRing

Build the Tambara–Yamagami fusion ring for a group with multiplication table `tab`.
Rank is n+1 (group elements + one extra object).
"""
function TY_fusion_ring(tab::AbstractMatrix{<:Integer}; names::Vector{String}=String[])
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


#Added: from izumi
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





#Added: from izumi

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
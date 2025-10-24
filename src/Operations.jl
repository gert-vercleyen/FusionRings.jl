
module Operations

using ..Types: FusionRing, fusion_tensor, labels, rank
using ..GeneralFunctions: indexmap
using LinearAlgebra: eigvals
using Combinatorics: permutations

export fusion_matrix, fusion_coeff, tensor_product, decompose, decompose_all
export multiplication_table, product_string, print_multiplication_table
export permute_mult_tab, permute, is_equivalent

labelstr(s) = String(s)

function fusion_matrix(fr::FusionRing, a)::Array{Int,2}
    A = fusion_tensor(fr)
    amap = indexmap(fr)
    idx = a isa Integer ? a : amap[String(a)]
    @views A[idx, :, :]
end

function fusion_coeff(fr::FusionRing, a, b, c)::Int
    imap = indexmap(fr)
    normalize(x) = x isa Integer ? x : imap[String(x)]
    ai = normalize(a); bi = normalize(b); ci = normalize(c)
    fusion_tensor(fr)[ai,bi,ci]
end

function tensor_product(fr::FusionRing, a, b)
    imap = indexmap(fr)
    normalize(x) = x isa Integer ? x : imap[String(x)]
    ai = normalize(a); bi = normalize(b)
    N = fusion_tensor(fr)[ai,bi,:]
    out = Dict{String,Int}()
    L = labels(fr)
    for (ci,m) in enumerate(N)
        m==0 && continue
        out[L[ci]] = m
    end
    out
end

decompose(fr::FusionRing, a, b) = [ (k,v) for (k,v) in tensor_product(fr,a,b) ]

function decompose_all(fr::FusionRing, a)
    out = Dict{String, Vector{Tuple{String,Int}}}()
    for b in labels(fr)
        out[b] = decompose(fr, a, b)
    end
    out
end

function multiplication_table(fr::FusionRing; include_zeros::Bool=false)
    L = labels(fr); r = length(L)
    cell(a::String, b::String) = begin
        d = tensor_product(fr, a, b)
        if include_zeros
            parts = String[]
            for c in L
                m = get(d, c, 0)
                if m==0
                    push!(parts, "0 "*labelstr(c))
                elseif m==1
                    push!(parts, labelstr(c))
                else
                    push!(parts, string(m," ",labelstr(c)))
                end
            end
            join(parts, " ⊕ ")
        else
            isempty(d) && return "0"
            join([ m==1 ? labelstr(c) : string(m," ",labelstr(c)) for (c,m) in d ], " ⊕ ")
        end
    end
    T = Array{String}(undef, r, r)
    for i in 1:r, j in 1:r
        T[i,j] = cell(L[i], L[j])
    end
    T
end

function print_multiplication_table(fr::FusionRing; include_zeros::Bool=false)
    L = labels(fr); r = length(L)
    head = "⊗ │ " * join(labelstr.(L), " │ ")
    sep  = "──┼" * "───┼"^(r-1) * "──"
    println(head); println(sep)
    T = multiplication_table(fr; include_zeros)
    for i in 1:r
        println(labelstr(L[i]), " │ ", join(T[i,:], " │ "))
    end
    nothing
end

function product_string(fr::FusionRing, a, b)
    L = labels(fr)
    amap = indexmap(fr)
    norm(x) = x isa Integer ? L[x] : String(x)
    aS = norm(a); bS = norm(b)
    rhs = multiplication_table(fr)[amap[aS], amap[bS]]
    string(labelstr(aS), " ⊗ ", labelstr(bS), " = ", rhs)
end

function permute_mult_tab(N::Array{Int,3}, p::Vector{Int})
    p[1]==1 || error("Permutation must fix the unit at index 1")
    r = size(N,1)
    M = fill(0, r, r, r)
    for a in 1:r, b in 1:r, c in 1:r
        M[p[a], p[b], p[c]] = N[a,b,c]
    end
    M
end

function permute(fr::FusionRing, p::Vector{Int})
    M = permute_mult_tab(fusion_tensor(fr), p)
    L = labels(fr)[invperm(p)]
    return FusionRing(M, L, fr.name*"/perm")
end

function is_equivalent(fr1::FusionRing, fr2::FusionRing)
    r1 = rank(fr1); r2 = rank(fr2)
    r1 == r2 || return false
    r = r1
    N1 = fusion_tensor(fr1); N2 = fusion_tensor(fr2)
    sum(N1) == sum(N2) || return false
    if r ≤ 8
        for p in permutations(2:r)
            perm = vcat(1, collect(p))
            if permute_mult_tab(N1, perm) == N2
                return true
            end
        end
        return false
    else
        S1 = zeros(Int, r, r); S2 = zeros(Int, r, r)
        for a in 1:r
            @views S1 .+= N1[a,:,:]
            @views S2 .+= N2[a,:,:]
        end
        sort(eigvals(Matrix(S1))) == sort(eigvals(Matrix(S2)))
    end
end

end


# Back-compat alias for table of strings
# const tensor_table = multiplication_table

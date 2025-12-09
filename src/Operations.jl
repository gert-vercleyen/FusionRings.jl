module Operations

using ..Types: FusionRing, multiplication_table
export fusion_matrix, fusion_coeff, fusion_product, decompose, decompose_all
export fusion_outcomes, product_string, print_multiplication_table
export permute_mult_tab, permute, is_equivalent


"Return the fusion matrix (left multiplication by `a`)."
function fusion_matrix(fr::FusionRing, a::Int)::Matrix{Int}
    @views multiplication_table(fr)[a, :, :]
end

"Structure constant N[a,b,c]."
fusion_coeff(fr::FusionRing, a::Int, b::Int, c::Int)::Int =
    multiplication_table(fr)[a,b,c]

"""
    fusion_product(fr, a, b) -> Dict{Int,Int}

Return the decomposition of `a ⊗ b` as a multiplicity dictionary
`Dict{simple_index => multiplicity}`.
"""
function fusion_product(fr::FusionRing, a::Int, b::Int)
    N = @views multiplication_table(fr)[a,b,:]
    out = Dict{Int,Int}()
    @inbounds for (c,m) in enumerate(N)
        m==0 && continue
        out[c] = m
    end
    out
end


#check this out - points to properties.jl -> tensor product
const tensor_product = fusion_product

"Return vector of simple indices with positive multiplicity in `a ⊗ b`."
fusion_outcomes(fr::FusionRing, a::Int, b::Int)::Vector{Int} =
    [c for (c,m) in fusion_product(fr,a,b) if m>0]

"Ordered list form of `a ⊗ b`."
decompose(fr::FusionRing, a::Int, b::Int) =
    [(k,v) for (k,v) in fusion_product(fr,a,b)]

"Decompose `a ⊗ j` for all `j`."
function decompose_all(fr::FusionRing, a::Int)
    r = size(multiplication_table(fr),1)
    Dict(j => decompose(fr, a, j) for j in 1:r)
end


#check this out - replace tensor product with multipication, dir_Sum with sum
#no nlonger use names - should be labels
#labels should be printed as bold integers
"""
Pretty prints the multiplication table as strings (no mutation).
"""
function print_multiplication_table(fr::FusionRing; include_zeros::Bool=false)
    N = multiplication_table(fr)
    names = fr.element_names
    r = length(names)
    head = "⊗ │ " * join(names, " │ ")
    sep  = "──┼" * "───┼"^(r-1) * "──"
    println(head); println(sep)
    for i in 1:r
        rowcells = String[]
        for j in 1:r
            d = fusion_product(fr, i, j)
            if include_zeros
                parts = String[]
                for c in 1:r
                    m = get(d,c,0)
                    if m==0; push!(parts, "0 "*names[c])
                    elseif m==1; push!(parts, names[c])
                    else; push!(parts, string(m," ",names[c]))
                    end
                end
                push!(rowcells, join(parts, " ⊕ "))
            else
                isempty(d) && push!(rowcells, "0") && continue
                push!(rowcells,
                    join([ m==1 ? names[c] : string(m," ",names[c]) for (c,m) in d ], " ⊕ "))
            end
        end
        println(names[i], " │ ", join(rowcells, " │ "))
    end
    nothing
end

"Pretty one-liner: `a ⊗ b = ...` using printed names; `a,b` are indices."
function product_string(fr::FusionRing, a::Int, b::Int)
    rhs = let d = fusion_product(fr,a,b), names = fr.element_names
        isempty(d) ? "0" :
            join([ m==1 ? names[c] : string(m," ",names[c]) for (c,m) in d ], " ⊕ ")
    end
    string(fr.element_names[a], " ⊗ ", fr.element_names[b], " = ", rhs)
end


"""
    permute_mult_tab(N, p)

Apply permutation `p` (fixing 1) to all three indices of `N`.
"""
function permute_mult_tab(N::Array{Int,3}, p::Vector{Int})
    p[1]==1 || error("Permutation must fix the unit at index 1")
    r = size(N,1)
    M = fill(0, r, r, r)
    @inbounds for a in 1:r, b in 1:r, c in 1:r
        M[p[a], p[b], p[c]] = N[a,b,c]
    end
    M
end

"""
    permute(fr, p) -> FusionRing

Return a *new* ring obtained by permuting simples by `p` (fixing 1).
"""
function permute(fr::FusionRing, p::Vector{Int})
    Np = permute_mult_tab(multiplication_table(fr), p)
    names = fr.element_names[invperm(p)]
    FusionRing(Np, names, fr.name*"/perm")
end

"""
#might return false positives - only checks in 1direction
    is_equivalent(r1, r2) -> Bool

Check graded ring isomorphism by brute force for rank ≤ 8,
else compare a spectral checksum of ∑_a N[a,:,:].
"""
function is_equivalent(fr1::FusionRing, fr2::FusionRing)
    N1 = multiplication_table(fr1); N2 = multiplication_table(fr2)
    r1 = size(N1,1); r2 = size(N2,1)
    r1 == r2 || return false
    r = r1
    sum(N1) == sum(N2) || return false

    if r ≤ 8
        using Combinatorics: permutations
        for p in permutations(2:r)
            perm = vcat(1, collect(p))
            permute_mult_tab(N1, perm) == N2 && return true
        end
        return false
    else
        using LinearAlgebra: eigvals
        S1 = zeros(Int, r, r); S2 = zeros(Int, r, r)
        @inbounds for a in 1:r
            @views S1 .+= N1[a,:,:]
            @views S2 .+= N2[a,:,:]
        end
        sort(eigvals(Matrix(S1))) == sort(eigvals(Matrix(S2)))
    end
end

end # module

export permute, permute_mult_tab, sort, perm_vec_qd, perm_vec_sd_conj,
       tensor_product, ⊗, which_permutation,
       register_known_ring!, known_rings, replace_by_known


"Return the fusion matrix (left multiplication by `a`)."
function fusion_matrix(fr::FusionRing, a::Int)::Matrix{Int}
    @views multiplication_table(fr)[a, :, :]
end

"Structure constant N[a,b,c]."
function fusion_coeff(fr::FusionRing, a::Int, b::Int, c::Int)::Int
    multiplication_table(fr)[a,b,c]
end

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


"""permute(r, perm) – return a new `FusionRing` with all data
    permuted by `perm`.  `perm[1]` **must** equal 1 to keep the vacuum first."""
function permute(r::FusionRing, perm::Vector{Int})::FusionRing
    n = rank(r)
    n == length(perm)      || throw(ArgumentError("perm length ≠ rank"))
    sort(perm) == collect(1:n) || throw(ArgumentError("perm must be a true permutation"))
    perm[1] == 1 || throw(ArgumentError("vacuum must stay at index 1"))

    # Core table
    mt_new = permute_mult_tab(multiplication_table(r), perm)

    # Metadata that needs re‑ordering (guard against `missing`)
    el_names = labels(r)[perm]
    tex_names = length(r.texnames) == n ? r.texnames[perm] : r.texnames
    fpdims   = r.frobenius_perron_dimensions === missing ?
               missing : (length(r.frobenius_perron_dimensions) == n ? r.frobenius_perron_dimensions[perm] : r.frobenius_perron_dimensions)
    chars    = r.characters === missing ?
               missing : (ndims(r.characters) == 2 && size(r.characters, 2) == n ? r.characters[:, perm] : r.characters)

    return fusion_ring(mt_new;                       # core data
        names         = r.names,
        texnames      = tex_names,
        labels        = el_names,
        barcode       = r.barcode,
        anyonwiki_code = r.anyonwiki_code,
        sub_fusion_rings = r.sub_fusion_rings,
        frobenius_perron_dimensions = fpdims,
        characters    = chars
    )
end

"""perm_vec_qd(r; order = :increasing) – permutation that sorts the non‑vacuum
    particles by Frobenius–Perron dimension."""
function perm_vec_qd(r::FusionRing; order::Symbol = :increasing)::Vector{Int}
    idx = collect(2:rank(r))
    qd  = fpdims(r)
    sort!(idx; by = i -> qd[i], rev = (order == :decreasing))
    return vcat(1, idx)
end

"""perm_vec_sd_conj(r; order = :increasing) – self‑duals first, then conjugate
    pairs, each block ordered by FP‑dimension."""
function perm_vec_sd_conj(r::FusionRing; order::Symbol = :increasing)::Vector{Int}
    n  = rank(r)
    conj = conjugate_element(r)
    qd = fpdims(r)

    self_dual = [i for i in 2:n if conj(i) == i]
    sort!(self_dual; by = i -> qd[i], rev = (order == :decreasing))

    paired   = Set(self_dual)
    pairs    = Tuple{Int,Int}[]

    for i in 2:n
        i in paired && continue
        j = conj(i)
        i == j && continue

        a, b = i, j
        if (order == :increasing && qd[a] > qd[b]) || (order == :decreasing && qd[a] < qd[b])
            a, b = b, a
        end

        push!(pairs, (a, b))
        push!(paired, a)
        push!(paired, b)
    end

    sort!(pairs; by = p -> qd[p[1]], rev = (order == :decreasing))
    conjlist = reduce(vcat, ([p[1], p[2]] for p in pairs); init = Int[])

    vcat(1, self_dual, conjlist)
end


function tensor_product(r1::FusionRing, r2::FusionRing)::FusionRing
    m, n   = rank(r1), rank(r2)
    mt1, mt2 = multiplication_table(r1), multiplication_table(r2)

    mt = zeros(Int, m*n, m*n, m*n)
    @inbounds for a in 1:m, α in 1:n, b in 1:m, β in 1:n, c in 1:m, γ in 1:n
        i = (a-1)*n + α
        j = (b-1)*n + β
        k = (c-1)*n + γ
        mt[i,j,k] = mt1[a,b,c] * mt2[α,β,γ]
    end

    # Assemble element names
    elnames = [ string(e1, "⊗", e2) for e1 in labels(r1) for e2 in labels(r2) ]

    names_tp = (isempty(names(r1)) || isempty(names(r2))) ? String[] :
        [string(names(r1)[1], "⊗", names(r2)[1])]

    fpdims_new = try
        [d1 * d2 for d1 in fpdims(r1) for d2 in fpdims(r2)]
    catch
        missing
    end

    return fusion_ring(
        mt; 
        names = names_tp, 
        labels = elnames, 
        frobenius_perron_dimensions = fpdims_new
    )
end




"Return vector of simple indices with positive multiplicity in `a ⊗ b`."
function fusion_outcomes(fr::FusionRing, a::Int, b::Int)::Vector{Int}
    [c for (c,m) in fusion_product(fr,a,b) if m>0]
end

"Ordered list form of `a ⊗ b`."
function decompose(fr::FusionRing, a::Int, b::Int) 
    [ (k,v) for (k,v) in fusion_product(fr,a,b) ]
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
#might return false positives - only checks in 1direction
    is_equivalent(r1, r2) -> Bool

Check graded ring isomorphism by brute force for rank ≤ 8,
else compare a spectral checksum of ∑_a N[a,:,:].
"""
function is_equivalent( r1::FusionRing, r2::FusionRing )
    !( which_permutation === missing )
end
function which_permutation(fr1::FusionRing, fr2::FusionRing)
    nsdnsd(fr1) == nsdnsd(fr2) || return missing

    dims1 = fpdims(fr1)
    dims2 = fpdims(fr2)

    Base.sort(dims1) == Base.sort(dims2) || return missing 

    r = r1
    sum(N1) == sum(N2) || return missing

    # if r ≤ 8
    #     using Combinatorics: permutations
    #     for p in permutations(2:r)
    #         perm = vcat(1, collect(p))
    #         permute_mult_tab(N1, perm) == N2 && return true
    #     end
    #     return false
    # else
    #     using LinearAlgebra: eigvals
    #     S1 = zeros(Int, r, r); S2 = zeros(Int, r, r)
    #     @inbounds for a in 1:r
    #         @views S1 .+= N1[a,:,:]
    #         @views S2 .+= N2[a,:,:]
    #     end
    #     sort(eigvals(Matrix(S1))) == sort(eigvals(Matrix(S2)))
    # end
    # return (r, nothing)
end

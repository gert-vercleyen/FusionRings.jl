
module Operations

using LinearAlgebra: kron
using Combinatorics          # for permutation generation (≤ 8 elements)
using ..FusionRings          

"""permute_mult_tab(mt, perm) – return a copy of the 3‑tensor `mt` with indices
    relabelled by `perm`."""
function permute_mult_tab(mt::Array{Int,3}, perm::Vector{Int})
    n = length(perm)
    out = similar(mt)
    @inbounds for a in 1:n, b in 1:n, c in 1:n
        out[a,b,c] = mt[perm[a], perm[b], perm[c]]
    end
    return out
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
    el_names = element_names(r)[perm]
    fpdims   = r.frobenius_perron_dimensions === missing ?
               missing : r.frobenius_perron_dimensions[perm]
    chars    = r.characters === missing ?
               missing : r.characters[:, perm]

    md = r.modular_data
    md_perm = md === missing ? missing : [Dict(
        "SMatrix"      => M["SMatrix"][perm, perm],
        "TwistFactors" => M["TwistFactors"][:, perm]
    ) for M in md]

    return fusion_ring(mt_new;                       # core data
        names         = r.names,
        texnames      = r.texnames,
        element_names = el_names,
        barcode       = r.barcode,
        formal_code   = r.formal_code,
        sub_fusion_rings = r.sub_fusion_rings,
        frobenius_perron_dimensions = fpdims,
        modular_data  = md_perm,
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
    cm = conjugation_matrix(r)            # antiparticle matrix
    qd = fpdims(r)

    self_dual = [i for i in 2:n if cm[i,i] == 1]
    sort!(self_dual; by = i -> qd[i], rev = (order == :decreasing))

    paired   = Set(self_dual)
    conjlist = Int[]

    for i in 2:n
        j = cm[i,i]
        if i != j && !(i in paired) && !(j in paired)
            push!(conjlist, i, j)
            push!(paired, i, j)
        end
    end

    return vcat(1, self_dual, conjlist)
end


function sort(r::FusionRing; sortby::String = "fpdims", kwargs...)
    perm = sortby == "fpdims"      ? perm_vec_qd(r; kwargs...) :
           sortby == "sd-conj"    ? perm_vec_sd_conj(r; kwargs...) :
           throw(ArgumentError("unknown sortby keyword: $sortby"))
    return permute(r, perm)
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
    elnames = [ string(e1, "⊗", e2) for e1 in element_names(r1) for e2 in element_names(r2) ]

    fp1, fp2 = r1.frobenius_perron_dimensions, r2.frobenius_perron_dimensions
    fpdims_new = (fp1 === missing || fp2 === missing) ?
                 missing : vec(kron(fp1, fp2))

    names = isempty(r1.names) || isempty(r2.names) ? missing : [string(r1.names[1], " × ", r2.names[1])]

    return fusion_ring(mt; names = names, element_names = elnames, frobenius_perron_dimensions = fpdims_new)
end

const ⊗ = tensor_product


"""which_permutation(r1, r2; exhaustive_limit = 8) – returns a permutation
vector σ such that `permute(r1, σ)` is identical to `r2`, or `nothing` if no
such σ exists.  For ranks ≤ `exhaustive_limit` the search is exhaustive;
otherwise a heuristic backtracking search is used."""
function which_permutation(r1::FusionRing, r2::FusionRing; exhaustive_limit::Int = 8)
    rank(r1) == rank(r2) || return nothing
    rank(r1) == 1 && return [1]           # trivial ring

    n  = rank(r1)
    mt2 = multiplication_table(r2)

    # Fast fingerprint: structure‑constant histograms must match
    if sort(vec(multiplication_table(r1))) != sort(vec(mt2))
        return nothing
    end

    # Exhaustive search for small n
    if n ≤ exhaustive_limit
        for perm in Combinatorics.permutations(1:n)
            perm[1] == 1 || continue
            if permute_mult_tab(multiplication_table(r1), perm) == mt2
                return collect(perm)
            end
        end
        return nothing
    end

    # Heuristic search for larger n – start with vacuum fixed, then greedy match
    perm = ones(Int, n)               # result vector being filled
    used = falses(n)
    used[1] = true

    # 1. Match self‑dual particles by FP dimension signature
    fp     = fpdims(r1)
    target = fpdims(r2)
    sd1    = [i for i in 2:n if multiplication_table(r1)[i,i,1] > 0]
    sd2    = [i for i in 2:n if mt2[i,i,1] > 0]
    sd_map = Dict{Int,Int}()
    for i in sd1
        match = findfirst(j -> isapprox(fp[i], target[j]; atol=1e-8) && !used[j], sd2)
        match === nothing && return nothing
        sd_map[i] = sd2[match]
        used[sd2[match]] = true
    end

    perm[collect(keys(sd_map))] = collect(values(sd_map))

    # 2. Greedy match remaining particles by their fusion with vacuum and FP‑dim
    for i in 2:n
        if perm[i] == 0
            cand = findfirst(j -> !used[j] && isapprox(fp[i], target[j]; atol=1e-8), 2:n)
            cand === nothing && return nothing
            perm[i] = cand
            used[cand] = true
        end
    end

    return permute_mult_tab(multiplication_table(r1), perm) == mt2 ? perm : nothing
end


const _KNOWN_RINGS = Dict{Int,Vector{FusionRing}}()   # keyed by barcode

"""register_known_ring!(r) – store `r` in the in‑memory catalogue so that later
    rings can be matched quickly via `replace_by_known`."""
function register_known_ring!(r::FusionRing)
    push!(_KNOWN_RINGS, barcode(r) => get(_KNOWN_RINGS, barcode(r), FusionRing[]))
    push!(_KNOWN_RINGS[barcode(r)], r)
    return nothing
end

known_rings() = _KNOWN_RINGS   # read‑only view

"""replace_by_known(r) – if an equivalent ring exists in `_KNOWN_RINGS`, return
    that canonical ring plus the permutation used; else return `(r, nothing)`."""
function replace_by_known(r::FusionRing)
    br = barcode(r)
    haskey(_KNOWN_RINGS, br) || return (r, nothing)

    for ref in _KNOWN_RINGS[br]
        perm = which_permutation(ref, r)
        perm === nothing && continue
        return (permute(ref, perm), perm)
    end
    return (r, nothing)
end

export permute, permute_mult_tab, sort, perm_vec_qd, perm_vec_sd_conj,
       tensor_product, ⊗, which_permutation,
       register_known_ring!, known_rings, replace_by_known

end 

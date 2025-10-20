
##############################
# GroupTheory.jl
# GAP/Oscar-backed functionality for fusion rings.
# This file integrates "Category 2" features that rely on group theory.
# Assumes Oscar.jl is available in the environment.
##############################

module GroupTheory

export fusion_ring_automorphisms,
       universal_grading,
       group_rep_fusion_ring

using ..Types
using ..Creation
using ..Operations
using ..Properties

# We only import Oscar at runtime inside functions so the package still loads
# on platforms where Oscar/GAP is not available.
const _OSCAR_HINT = "This function requires Oscar.jl (which bundles GAP). Please `using Oscar` in a Linux/WSL environment."

"""
    _require_oscar()

Internal helper. Throws a helpful error if Oscar cannot be loaded.
"""
function _require_oscar()
    try
        @eval using Oscar
        return nothing
    catch e
        error(_OSCAR_HINT * "\nOriginal error: $(sprint(showerror, e))")
    end
end

"""
    group_rep_fusion_ring(G; labelstyle = :index)

Construct the fusion ring of **irreducible representations** of a finite group `G`
using Oscar/GAP character tables. The simples are the irreps, fusion coefficients
are multiplicities in tensor products of irreducibles computed via the standard
inner product on characters.

- `labelstyle = :index`   → labels are `"χ1"`, `"χ2"`, …
- `labelstyle = :degree`  → labels are `"χ(1)"`, `"χ(2)"`, … where the number in parentheses is the degree

This mirrors the character ring construction in Anyonica.
"""
function group_rep_fusion_ring(G; labelstyle=:index)
    _require_oscar()
    # Access Oscar/GAP character table
    T = Oscar.CharacterTable(G)
    irreps = T.irreducibles  # GAP characters
    n = length(irreps)

    # Build labels
    deg = [ Int(Oscar.GAP.Globals.DegreeOfCharacter(c)) for c in irreps ]
    labels = if labelstyle === :degree
        ["χ($(deg[i]))" for i in 1:n]
    else
        ["χ$(i)" for i in 1:n]
    end

    # The multiplication tensor N[i,j,k] = < χ_i * χ_j, χ_k >
    N = fill(0, n, n, n)
    for i in 1:n, j in 1:n
        prod = Oscar.GAP.Globals.CharacterProduct(irreps[i], irreps[j])
        for k in 1:n
            # multiplicity is the inner product of prod with χ_k
            mult = Int(Oscar.GAP.Globals.ScalarProduct(T, prod, irreps[k]))
            N[i,j,k] = mult
        end
    end

    return fusion_ring(N; labels, name = "Rep(" * string(G) * ")")
end

"""
    fusion_ring_automorphisms(R; method=:auto, maxrank_bruteforce=8)

Return a description of the **label automorphism group** of `R`,
i.e. permutations of simples preserving the multiplication tensor.

- `method=:auto` → use GAP when available; fall back to brute force for small ranks.
- `method=:bruteforce` → try all permutations up to `maxrank_bruteforce`.
- `method=:gap` → compute via GAP/Oscar (recommended for rank ≥ 8).

Returns a struct-like `Dict` containing:
- `:generators` → permutations as vectors of integers (1-based)
- `:size`       → group size (if GAP path used)
- `:gap`        → GAP group object (if available)
"""
function fusion_ring_automorphisms(R::FusionRing; method=:auto, maxrank_bruteforce=8)
    r = rank(R)
    N = fusion_tensor(R)

    if method === :bruteforce || (method === :auto && r <= maxrank_bruteforce)
        # Brute force search: keep permutations p that satisfy N[p(i), p(j), p(k)] == N[i,j,k]
        using Combinatorics
        gens = Int[][]
        # quick generators: find all label swaps that preserve unit row/column
        unit = 1
        for p in permutations(1:r)
            ok = true
            if p[unit] != unit
                continue
            end
            @inbounds for i in 1:r, j in 1:r, k in 1:r
                if N[p[i], p[j], p[k]] != N[i,j,k]
                    ok = false
                    break
                end
            end
            if ok
                push!(gens, collect(p))
            end
        end
        return Dict(:generators => gens, :size => length(gens))
    else
        _require_oscar()
        # Provide GAP with a "colored" 3-tensor via list of structure constants
        # We encode the set S = { (i,j,k) | N[i,j,k] != 0 } and compute its automorphisms.
        triples = Oscar.GAP.GapObj([Oscar.GAP.GapObj([i,j,k]) for i in 1:r for j in 1:r for k in 1:r if N[i,j,k] != 0])
        # Symmetric group on r letters acts by permuting coordinates of each triple
        S_r = Oscar.GAP.Globals.SymmetricGroup(r)
        act = Oscar.GAP.Globals.Action(S_r, triples, Oscar.GAP.Globals.OnTuples)
        G = Oscar.GAP.Globals.Stabilizer(act, triples)
        # Extract generators as Julia permutations
        gens = Oscar.GAP.Globals.GeneratorsOfGroup(G)
        toperm = p -> begin
            v = [Int(Oscar.GAP.Globals.Image(p, i)) for i in 1:r]
            v
        end
        return Dict(:generators => [ toperm(g) for g in gens ],
                    :size => Int(Oscar.GAP.Globals.Size(G)),
                    :gap => G)
    end
end

"""
    universal_grading(R::FusionRing)

Compute a **coarse universal grading** of `R` by the adjoint subring:
group simples by their cosets modulo the adjoint closure (i.e. generated by `x⊗x̄`).

Returns a dictionary:
- `:grading` → Dict{String,Int} mapping labels to grade indices
- `:numgrades` → number of grading components
- `:adjoint_irreps` → labels contained in the adjoint
This matches the idea of universal grading in Anyonica; a fully functorial construction
can be added via GAP later.
"""
function universal_grading(R::FusionRing)
    A = Set(adjoint_irreps(R))
    labs = labels(R)
    # naive partition: grade 0 = adjoint; other grades by multiplication cosets
    grade = Dict{String,Int}()
    for a in labs
        grade[a] = (a in A) ? 0 : 1
    end
    # refine: if a ⊗ ā contains unit, keep in grade 0
    for a in labs
        b = conjugate_label(R, a)
        for (c,m) in decompose(R, a, b)
            if c == labs[1] && m>0
                grade[a] = 0
            end
        end
    end
    ngrades = maximum(values(grade); init=0) + 1
    return Dict(:grading => grade, :numgrades => ngrades, :adjoint_irreps => collect(A))
end

end # module

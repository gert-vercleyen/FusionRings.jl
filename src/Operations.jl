module Operations
using LinearAlgebra
using ..Types, ..Creation

export fusion_matrix, fusion_coeff, tensor_product

export decompose, decompose_all, tensor_table, print_tensor_table

function fusion_matrix(fr::Types.FusionRing, a)
    ia = Creation._idx(fr,a)
    r = size(fr.N,1)
    M = Matrix{Int}(undef, r, r)
    @inbounds for b in 1:r, c in 1:r
        M[b,c] = fr.N[ia,b,c]
    end
    return M
end

fusion_coeff(fr::Types.FusionRing,a,b,c) = fr.N[Creation._idx(fr,a),Creation._idx(fr,b),Creation._idx(fr,c)]

function tensor_product(fr::Types.FusionRing,a,b)
    ia, ib = Creation._idx.(Ref(fr), (a,b))
    r = size(fr.N,1)
    out = Dict{Symbol,Int}()
    @inbounds for c in 1:r
        n = fr.N[ia,ib,c]
        n==0 && continue
        out[fr.labels[c]] = n
    end
    return out
end

end # module


# Return the decomposition a⊗b as a list of (label, mult).
function decompose(fr::Types.FusionRing, a, b; include_zeros::Bool=false)
    ia, ib = Creation._idx.(Ref(fr), (a,b))
    r = size(fr.N,1)
    out = Vector{Tuple{Symbol,Int}}()
    @inbounds for c in 1:r
        n = fr.N[ia,ib,c]
        (include_zeros || n!=0) && push!(out, (fr.labels[c], n))
    end
    return out
end

# Return a dictionary mapping each right factor to its decomposition with `a`.
function decompose_all(fr::Types.FusionRing, a; include_zeros::Bool=false)
    r = size(fr.N,1)
    res = Dict{Symbol, Vector{Tuple{Symbol,Int}}}()
    for b in 1:r
        res[fr.labels[b]] = decompose(fr, a, b; include_zeros)
    end
    return res
end

# Pretty table (matrix of strings) for all a,b
function tensor_table(fr::Types.FusionRing; include_zeros::Bool=false)
    r = size(fr.N,1)
    tbl = Array{String}(undef, r, r)
    for a in 1:r, b in 1:r
        parts = String[]
        for (lab, n) in decompose(fr, a, b; include_zeros)
            if n==0 && !include_zeros; continue; end
            push!(parts, n==1 ? string(lab) : string(n)*" "*string(lab))
        end
        tbl[a,b] = isempty(parts) ? "0" : join(parts, " ⊕ ")
    end
    return tbl
end


export product_string

# Return a pretty string like "a ⊗ b = x ⊕ 2 y ⊕ ..."
function product_string(fr::Types.FusionRing, a, b)
    parts = String[]
    for (lab, n) in decompose(fr, a, b)
        n == 0 && continue
        push!(parts, n==1 ? string(lab) : string(n) * " " * string(lab))
    end
    return string(a) * " ⊗ " * string(b) * " = " * (isempty(parts) ? "0" : join(parts, " ⊕ "))
end


# Pretty-print the full tensor table with row/column headers.
function print_tensor_table(fr::Types.FusionRing; include_zeros::Bool=false)
    labs = string.(Types.labels(fr))
    T = tensor_table(fr; include_zeros=include_zeros)

    # Build rows with headers
    rows = Vector{Vector{String}}(undef, length(labs) + 1)
    rows[1] = vcat(["⊗"], labs)                  # header row
    for i in 1:length(labs)
        rows[i+1] = vcat([labs[i]], T[i, :])     # row label + data
    end

    # Column widths
    ncol = length(rows[1])
    widths = [maximum(length.(getindex.(rows, j))) for j in 1:ncol]

    # printing helpers
    function _pad(s, w)
        n = w - length(s)
        n <= 0 && return s
        return s * " "^n
    end

    # Print header + separator + rows
    header = join([_pad(rows[1][j], widths[j]) for j in 1:ncol], " │ ")
    println(header)
    sep = join(["─"^widths[j] for j in 1:ncol], "─┼─")
    println(sep)
    for i in 2:length(rows)
        line = join([_pad(rows[i][j], widths[j]) for j in 1:ncol], " │ ")
        println(line)
    end
    nothing
end

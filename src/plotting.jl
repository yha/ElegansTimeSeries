using Plots
using BlockArrays
using NamedArrays
using StatsBase: midpoints

# heatmap for block-matrix of binned matrix blocks with named blocks
labeled_binmat_heatmap(m; block_separator_color = "white", kwargs...) = _labeled_binmat_heatmap(m; block_separator_color, kwargs...)

#function _labeled_binmat_heatmap(m::AbstractBlockArray; block_separator_color, kwargs...)
function _labeled_binmat_heatmap(m; block_separator_color, kwargs...)
    res = heatmap(m; kwargs...)
    ends1 = blocklasts(axes(m,1))
    ends2 = blocklasts(axes(m,2))
    hline!(ends1 .+ 1/2, c=block_separator_color, label="")
    vline!(ends2 .+ 1/2, c=block_separator_color, label="")
    if (n = _names(blocks(m))) !== nothing
        yticks!((midpoints([0; ends1]), n[1]))
        xticks!((midpoints([0; ends2]), n[2]))
    end
    res
end
#_labeled_binmat_heatmap(m; kwargs...) = heatmap(m; kwargs)
_names(a::NamedArray) = names(a)
_names(a) = nothing

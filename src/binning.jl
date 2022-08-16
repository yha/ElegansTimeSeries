using Statistics
#using OnlineStats

function bins(ax, nbins)
    edges = round.(Int, range(first(ax), last(ax)+1, length=nbins+1))
    [edges[i]:edges[i+1]-1 for i in 1:nbins]
end
fin(x) = !ismissing(x) && isfinite(x)
bin_means(v::AbstractVector, nbins) = [mean(view(v,i)) for i in bins(eachindex(v), nbins)]
bin_mean_fin(v::AbstractVector, nbins) = [mean(filter(fin,view(v,i))) for i in bins(eachindex(v), nbins)]

# # naive implementation
# function bin_means_weights(v, nbins)
#     b = bins(eachindex(v), nbins)
#     f = [filter(fin, view(v,i)) for i in b]
#     means = mean.(f)
#     weights = length.(f) ./ length.(b)
#     (; means, weights)
# end

function bin_means_weights(v, nbins)
    means = similar(v, (nbins,))
    weights = similar(v, (nbins,))
    bin_means_weights!(means, weights, v, nbins)
    (; means, weights)
end

# fast non-allocating computation of mean after filtering
# and weight ( = proportion filtered ) in each bin
function bin_means_weights!(m, w, v, nbins; filter_f = fin)
    b = bins(eachindex(v), nbins)
    @assert length(b) == nbins
    for i in eachindex(m,w,b)
        indexes = b[i]
        bin = view(v, indexes)
        bin_sum, bin_count = zero(eltype(v)), 0
        for x in bin
            filter_f(x) || continue
            bin_count += 1
            bin_sum += x
        end
        m[i] = bin_sum / bin_count
        w[i] = bin_count / length(indexes)
    end
    m, w
end

function bin_means_weights_per_col(a, nbins)
    means = similar(a, (nbins, size(a,2)))
    weights = similar(a, (nbins, size(a,2)))
    for i in axes(a,2)
        bin_means_weights!(view(means,:,i), view(weights,:,i), view(a,:,i), nbins)
    end
    (; means, weights)
end

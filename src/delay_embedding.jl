using JuliennedArrays

"""
	y, indices = delay_embed(a, k; dim=ndims(a), cond)
Delay-embedding (DE) matrix from an array `a` with time along dimension `dim`.
Each column of `y` is a vectorized slice of `a` of length `k` along dimension `dim`.
If `cond` is specified, only slices where `cond(slice)` holds are included.
Returns the DE matrix `y` and corresponding indices into `axes(a,dim)` (indices 
where `cond` holds), where an index of 1 corresponds to the first possible slice.
"""
function delay_embed(a, k; dim=ndims(a), cond=Returns(true))
	# slice(i) = selectdim(a,dim,i+1:i+k)
	# min_i, max_i = firstindex(a,dim), lastindex(a,dim)
	# slices = [vec(slice(i)) for i in (min_i - 1):(max_i - k)]
	# i = findall(cond, slices)
	# out = Align(slices[i], 1; slice_axes = (1:(length(a)÷size(a,dim))*k,))
	# out, i
    vecs, i = delay_embed_vecs(a, k; dim, cond)
    aligned = Align(vecs, 1; slice_axes = (1:(length(a)÷size(a,dim))*k,))
    aligned, i
end

function delay_embed_vecs(a, k; dim=ndims(a), cond=Returns(true))
    slice(i) = selectdim(a,dim,i+1:i+k)
	min_i, max_i = firstindex(a,dim), lastindex(a,dim)
	slices = [vec(slice(i)) for i in (min_i - 1):(max_i - k)]
	i = findall(cond, slices)
	slices[i], i
end

allfinite(v) = all(isfinite,v)

function well_bin_de(experiment, wellname, statname, nbins, bin_i, winlen; 
                     stagedict = loadstages(),
                     mid_stats,
                     midpoints_cache, 
                     n_trim, 
                     cond = ElegansTimeSeries.allfinite)
    mids = midpoints_cache(experiment, wellname)

    stage_i, bin_in_stage = fldmod1(bin_i, nbins)
    bin_range = let
        stage_ends = stagedict[experiment][wellname]
        stage_range = stage_ends[stage_i]+1:stage_ends[stage_i+1]
        stage_bin_ranges = bins(stage_range, nbins)
        stage_bin_ranges[bin_in_stage]		
    end
    bin_mids = mids[stage_i][Base.IdentityUnitRange(bin_range), :]

    mat = mid_stats(bin_mids)[statname][:,begin+n_trim:end-n_trim]

    de, indices = delay_embed(mat', winlen; cond)
    #delay_embedding_pca(mat', winlen)
    # window_μ = mean(posture_pca)
    # window_P = window_pca_rescale ? loadings(posture_pca) : projection(posture_pca)
end

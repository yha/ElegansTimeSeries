using MultivariateStats
using StatsBase: sample
using ProgressLogging

function split_mat(mat, p=0.1)
    i = rand(size(mat,1)) .< √p # rand(Bool, size(mat,1))
    j = rand(size(mat,2)) .< √p # rand(Bool, size(mat,2))
    A, B, C, D = mat[i,j], mat[i,.!j], mat[.!i,j], mat[.!i, .!j]
end

function pca_cv_decomposition(B,C,D; maxoutdim)
    pca = fit(PCA, D; maxoutdim, pratio = 1)
    D_scores = MultivariateStats.transform(pca, D)
    D_axes = projection(pca)
    C_scores = D_axes \ C
    B_axes = B / D_scores
    (; B_axes, C_scores)
end

function cv_pca_errs(mat, nshuffles, ndims)
    function pca_cv_err(A, B_axes, C_scores, dim)
        mean((A .- B_axes[:,1:dim] * C_scores[1:dim,:]) .^ 2)
    end

    errs = zeros(only(axes(ndims)), nshuffles)
    maxoutdim = maximum(ndims)
    # var_A = fill(NaN, ndims, nshuffles)
    @progress for i in 1:nshuffles
        A,B,C,D = split_mat(mat)
        err_denom = var(A)
        #var_A[:,i] = cumsum(principalvars(fit(PCA,A))[dims]) ./ (size(A,1) * err_denom)
        #pca = fit(PCA,A)
        #var_A[:,i] = 1 .- cumsum(principalvars(pca)[1:ndims]) ./ tvar(pca)
        B_axes, C_scores = pca_cv_decomposition(B,C,D; maxoutdim)
        out_d = size(B_axes, 2)
        @assert out_d == size(C_scores, 1)
        for (k, dim) in pairs(ndims)
            errs[k,i] = if dim > out_d
                NaN
            else
                pca_cv_err(A, B_axes, C_scores, dim) ./ err_denom
            end
        end
    end
    # (; errs)#, var_A)
    errs
end	

function logdecrease_thr_d_estimates(e, n_resamples, threshold)
	_, asamples = axes(e.err)
	
	map(1:n_resamples) do _
		i = sample(asamples, length(asamples))
		μ = mean(e.err[:,i]; dims = 2)[:,1]
		decrease = diff(.-log10.(μ)) ./ diff(e.ndims)
		e.ndims[something(findfirst(<(threshold), decrease), length(decrease)) + 1]
	end
end
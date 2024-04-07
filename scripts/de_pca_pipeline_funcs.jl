function depca_vars(de_mat, pca)
	if isempty(de_mat)
		vars = fill(NaN, outdim(pca))
		tvar = NaN
	else
		tvar = sum(var, eachrow(de_mat))
		mat_t = MultivariateStats.transform(pca, de_mat)
		vars = var.(eachrow(mat_t))
	end
	(; vars, tvar, nwindows = size(de_mat,1))
end

sampledim(a, nsamples, d; replace=true, ordered=false) = 
            selectdim(a, d, sample(axes(a,d), nsamples; replace, ordered))
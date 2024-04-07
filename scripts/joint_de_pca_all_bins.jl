# Compute population DE-PCA for all bins jointly by sampling some windows from each bin,
# then compute individuals' score variances and cross-validation errors as in the main
# pipeline.
# Used as a one-off for N2 population. If needed more often or for more conditions, 
# probably should be merged into main pipeline at `de_pca_pipeline_dist.jl` to distribute
# conditions across worker processes.

using Elegans
using ElegansTimeSeries
using MultivariateStats
using StatsBase
using BlockArrays
using JLD2
using ProgressLogging
using Unzip

include("de_pca_pipeline_funcs.jl")

#conditions = [(strain = "N2", ds = 0)]
condition = (strain = "N2", ds = 0)
condname = "$(condition.strain) $(condition.ds)"

nbins = 10
nbins_total = nbins * length(Elegans.stage_names)

pca_winlen = 60
ntrim = 1
confidence_threshold = 0.05

# windows sampled per bin for population PCA
nsamples_per_bin = 10000

# PCA cross-validation parameters
cv_pca_nsamples = 10
#cv_pca_nsamples = 1
cv_pca_nwindows = 1000

de_path = "U:/cached-data/DEs"
outdir = "U:/cached-data/DE-PCA"

##


# TODO move filename logic inside package? (repeated in `de_pca_pipeline_dist.jl`)
condname = "$(condition[1]) $(condition[2])DS"
path = "$de_path/DEs with conf $condname (nbins=$nbins, winlen=$pca_winlen, ntrim=$ntrim).jld2"
@assert isfile(path)
@info "Loading DEs from file `$path`"
@time d = load(path)

# TODO list of wells taken from DE file prevents skipping loading of DE file when PCs 
# are already computed. Share well list logic with `de_pca_pipeline_dist`?
des = d["pop_des"]
wells = d["wells"]

##

nwells = length(wells)
@assert nwells == size(des,2)

depca_params = (; nwells, ntrim, pca_winlen, confidence_threshold)

# Output paths are modified versions of those defined by functions in `src/paths.jl`
outpath_depca = "$outdir/DE-PCA all bins $condname, $nsamples_per_bin samples per bin, $(ElegansTimeSeries.meta_string(depca_params)).jld2"
outpath_vars = "$outdir/DE-PCA vars all bins $condname, $nsamples_per_bin samples per bin, $(ElegansTimeSeries.meta_string(depca_params)).jld2"
outpath_pop_vars = "$outdir/DE-PCA population vars all bins $condname, $nsamples_per_bin samples per bin, $(ElegansTimeSeries.meta_string(depca_params)).jld2"
outpath_cv = "$outdir/CV PCA errors all bins, $nsamples_per_bin samples per bin, $cv_pca_nwindows windows, $cv_pca_nsamples samples, $condname, $(ElegansTimeSeries.meta_string(depca_params)).jld2" 

##

@info "Sampling $nsamples_per_bin DE windows from each time bin"
@time begin
    all_des_each_bin = [mortar(reshape(r,1,:)) for r in eachrow(des)]

    # Seperate `map` and then `reduce` seems to be faster
    # than `mapreduce` here
    sampled_des_per_bin = map(all_des_each_bin) do m
        i = sample(axes(m, 2), nsamples_per_bin; replace = false)
        m[:,i]
    end

    sampled_des = reduce(hcat, sampled_des_per_bin)
end

##

@info "Fitting PCA..."

@time joint_pca = fit(PCA, sampled_des)

jldsave(outpath_depca; joint_pca)

##

@info "Computing population score variances..."

pop_score_vars = let
    # iterate over wells in the outer loop to use the midpoints cache 
    # more efficiently, then permute dims 
    @progress "DEPCA variances" res = [ depca_vars( reduce(hcat, des[bin,:]), joint_pca ) 
        for bin in 1:nbins_total
    ]
    vars, tvars, nwindows = unzip(permutedims(res))
    (; vars, tvars, nwindows)
end

jldsave(outpath_pop_vars; 
        pop_score_vars.vars, pop_score_vars.tvars, pop_score_vars.nwindows,
        wells
)

##

@info "Computing individual score variances..."

score_vars = let
    # iterate over wells in the outer loop to use the midpoints cache 
    # more efficiently, then permute dims 
    @progress "DEPCA variances" res = [ depca_vars( convert(Matrix, des[bin,well_i]), joint_pca ) 
        for bin in 1:nbins_total, well_i in 1:nwells
    ]
    vars, tvars, nwindows = unzip(permutedims(res))
    (; vars, tvars, nwindows)
end

jldsave(outpath_vars; 
        score_vars.vars, score_vars.tvars, score_vars.nwindows,
        wells
)

##

@info "Computing PCA cross-validation errors..."

ndims = [1:50; 60:10:400]

# Sample a subset of `cv_pca_nwindows` windows from the DEs sampled for the full PCA,
# and collect to a matrix
mat = Matrix(sampledim(sampled_des, cv_pca_nwindows, 2; replace = false))
err = cv_pca_errs(mat, cv_pca_nsamples, ndims)
#(; ndims, err)
    

@info "Saving CV PCA errors to `$outpath_cv`"
jldsave(outpath_cv; ndims, err)
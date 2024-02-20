using JLD2, FileIO
using TOML

function data_paths(; expected_params = ["video", "contours", "midpoints", "area", "depca", "depca_individual", "depca_bs", "stages"])
	paths_dict = TOML.parsefile("$(@__DIR__)/../paths.toml")["paths"]
	@assert keys(paths_dict) == Set(expected_params)

	NamedTuple((Symbol(k), normpath(joinpath(@__DIR__, v))) for (k,v) in paths_dict)
end

paths = data_paths()

function load_depcas(source)
	d = load(source)
	pcas = d["pcas"]
	indices = d["de_indices"]
	if haskey(d, "conditions")
		conds = d["conditions"]
		wells = d["wells"]
	else
		# Sanitized data format: 
		# single condition specified as strain& DS, and well count in each experiment rather than well list
		strain, ds = d["strain"], d["DS"]
		n_individuals = d["n_individuals"] # vector of pairs of (experiment, n)
		conds = [(strain, ds)]
		wells = [(; experiment, well = "well $i") for (experiment,n) in n_individuals for i = 1:n]
	end
	@info "Loaded $(length(pcas)) PCAs for $(join(join.(conds, "-"), ", ")) (n=$(length(wells)))"
	(; pcas, indices, conds, wells)
end

function make_depca_cache(; cache = Dict{String, Any}(), path = data_paths().depca)
    return get_depcas!(source) = get!(cache, source) do
        @info "Loading PCAs from `$source`"
        load_depcas(joinpath(path, source))
    end
end

meta_string(params) = meta_string(; params...)
meta_string(; nwells, ntrim, pca_winlen, confidence_threshold) = 
	"n=$nwells, n_trim=$ntrim, windowlength=$pca_winlen, conf_th=$confidence_threshold"

depca_filename(condname, depca_params) = "DE-PCA $condname, $(meta_string(depca_params)).jld2"
depca_vars_filename(condname, depca_params)  = "DE-PCA vars $condname, $(meta_string(depca_params)).jld2"
depca_cv_filename(condname, nwindows, nsamples, depca_params) = 
	"CV PCA errors $condname, $nwindows windows, $nsamples sample, $(meta_string(depca_params)).jld2" 
depca_individual_filename(well, maxoutdim, depca_params) = 
	"DE-PCA $(well.experiment), $(well.well), $(meta_string(depca_params)), maxoutdim=$maxoutdim.jld2"
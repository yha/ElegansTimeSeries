using JLD2, FileIO
"""
A somewhat hacky way to load correct paths for various machines
"""
paths = let hostname = gethostname()
    if hostname == "bi-sternlab21"
        coord_stats = "H:/coord-stats-binned"
        video = "U:/experiments/reemy"
        contours = "U:/cached-data/contours"
        midpoints = "F:/cached-data/midpoints"
        pop_des = "E:/cached-data/DEs"
        depca = "H:/cached-data/DE-PCA"
        
    elseif hostname == "bi-sternlab33"
        coord_stats = "U:/cached-data/coord-stats-binned"
        video = "U:/experiments/reemy"
        contours = "U:/cached-data/contours"
        midpoints = "U:/cached-data/midpoints"
        pop_des = "U:/cached-data"
        depca = "U:/cached-data/DE-PCA"

    elseif hostname == "BI-Sternlab3"
        coord_stats = "U:/cached-data/coord-stats-binned"
        video = "U:/experiments/reemy"
        contours = "U:/cached-data/contours"
        midpoints = "F:/midpoints"
        pop_des = "F:/cached-data/DEs"
        depca = "F:/cached-data/DE-PCA"

    else
        coord_stats = "../../coord-stats-binned"
        video =  "$(homedir())/experiments/reemy"
        contours = "$(homedir())/cached-data/contours"
        midpoints = "$(homedir())/cached-data/midpoints"
        pop_des = "$(homedir())/cached-data/DEs"
        depca = "$(homedir())/cached-data/DE-PCA"        
    end

    depca_individual = "$depca/individuals"
    depca_bs = "$depca/bootstraps"
    (; #shape_base, shape, 
        depca, depca_individual, depca_bs,
        coord_stats, video, contours, midpoints, pop_des
        )
end

function load_depcas(source)
	d = load(source)
	pcas = d["pcas"]
	conds = d["conditions"]
	indices = d["de_indices"]
	wells = d["wells"]
	@info "Loaded $(length(pcas)) PCAs for $(join(join.(conds, "-"), ", ")) (n=$(length(wells)))"
	(; pcas, indices, conds, wells)
end
function make_depca_cache()
    #cache = LRU(; maxsize = 20)
    cache = Dict{String,Any}()
    return get_depcas!(source) = get!(cache, source) do
        @info "Loading PCAs from `$source`"
        res = load_depcas(joinpath(paths.depca, source))
        #@assert length(res.pcas) == nbins_total
        #res
    end
end

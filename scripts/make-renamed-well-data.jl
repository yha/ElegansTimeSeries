"""
Create a "name-sanitized" copy of coordinates, midlines and individuals' DE-PCA data for publishing,
where experiments and wells are referred to by consecutive numerals.
Also saves the mapping from lab-local name to sanitized names/numerals in CSV files.

Input midpoints directory should have a subdirectory per condition.
Other (non-individual) DE-PCA data is not copied by this script.
"""

using CSV, DataFrames
using Elegans
using OrderedCollections
using FileIO, JLD2
using ProgressLogging
using MultivariateStats # For loading PCAs from JLD2 files

well_list_from_notebook_path = "G:/cached-data/all-depca-wells.csv"

midpoints_path = "G:/cached-data/midpoints/"
depca_path     = "G:/cached-data/DE-PCA/"
individual_depca_path = "$depca_path/individuals"

out_path         = "G:/sanitized/data"
out_path_mapping = "G:/sanitized"

coordinate_paths = ["U:/experiments/reemy", "U:/experiments/manal/Coords"]
ex_list_files = "$(@__DIR__)/../" .* ["exp-list-filtered-RA.csv", "exp-list-filtered-NG.csv"]

##


# Set of contour methods that should have midpoint files
contour_methods = [Thresholding(1.0,0.35), Thresholding(1.0,0.34)]


##

ex_dfs = [DataFrame(CSV.File(f)) for f in ex_list_files]
exs_from_lists = DataFrame(
    name   = reduce(vcat, df.name for df in ex_dfs),
    strain = reduce(vcat, df.strain for df in ex_dfs),
    DS     = reduce(vcat, df.DS for df in ex_dfs)
)
@assert allunique(exs_from_lists.name)

##

stagedict = loadstages()

##

conditions = unique((; r.strain, r.DS) for r in eachrow(exs_from_lists))

let
    exs_without_stages = filter(ex -> !haskey(stagedict,ex), exs_from_lists.name)
    if !isempty(exs_without_stages)
        @info "Some experiments have no stage data" exs_without_stages
    end
end


## Get list of wells for each experiment from midpoint file names

# condition name as used in midpoint directory names
condname(cond) = cond.DS > 0 ? "$(cond.strain) $(cond.DS)DS" : cond.strain
# full name is used in individual DE-PCA file names
condname_full(cond) = "$(cond.strain) $(cond.DS)DS"

function midpoint_filename_parse(name)
    m = match(r"midpoints-(.+?)-(cam.+?) (.+?)-(.+?) \(t = (.+?), s = (.+?)\) .*.jld"i, name)
    m === nothing && return nothing
    experiment, well, threshold, σ, t, s = m.captures
    threshold, σ, t, s = parse.(Float64, (threshold, σ, t, s))
    (; experiment, well, threshold, σ, t, s)
end

cond2wells = map(conditions) do cond
    files = readdir("$midpoints_path/$(condname(cond))")
    midfile_meta = midpoint_filename_parse.(files)
    # all midpoint files should use the same parameters outside of contour method
    @assert length(unique((m.s, m.t) for m in midfile_meta)) == 1
    wells = unique((; m.experiment, m.well) for m in midfile_meta)
    # each well should have one midpoint file for each method in `contour_methods`
    contour_params = Set((c.level, c.σ) for c::Thresholding in contour_methods)
    for (; experiment, well) in wells
        well_mid_meta = filter( m -> (m.experiment, m.well) == (experiment, well), midfile_meta )
        if Set((m.threshold, m.σ) for m in well_mid_meta) != contour_params
            error("Contour params mismatch: $(condname(cond)): $experiment-$well. Found midpoints for $([(m.threshold, m.σ) for m in well_mid_meta])")
        end
    end
    cond => wells
end |> OrderedDict
##
cond2exs = OrderedDict(cond => unique(w.experiment for w in wells) for (cond, wells) in cond2wells)

exs_with_midfiles = DataFrame(
    (; cond.strain, cond.DS, experiment) for cond in conditions for experiment in cond2exs[cond]
)


##
let 
    no_midfile_exs   = Vector{String}(setdiff(exs_from_lists.name, exs_with_midfiles.experiment))
    i = exs_from_lists.name .∈ Ref(no_midfile_exs)
    no_midfile_conds = [(; r.strain, r.DS) for r in eachrow(exs_from_lists[i,:])]
    isempty(no_midfile_exs) || @info """
        Experiments without midpoint files: $(join(("$ex ($(condname(cond)))" for (ex, cond) in zip(no_midfile_exs, no_midfile_conds)), ", "))
    """
    only_midfile_exs   = Vector{String}(setdiff(exs_with_midfiles.experiment, exs_from_lists.name))
    i = exs_with_midfiles.experiment .∈ Ref(only_midfile_exs)
    only_midfile_conds = [(; r.strain, r.DS) for r in eachrow(exs_with_midfiles[i,:])]
    isempty(only_midfile_exs) || @info """
        Experiments with midpoint files not in experiment list: $(join(("$ex ($(condname(cond)))" for (ex, cond) in zip(only_midfile_exs, only_midfile_conds)), ", "))
    """
end


##
# This should match list generated from notebook
wells_with_midfiles_df = DataFrame(
    (; cond.strain, cond.DS, well.experiment, well.well) for cond in conditions for well in cond2wells[cond]
)
mkpath(out_path_mapping)
CSV.write("$out_path_mapping/wells_with_midpoint_files.csv", wells_with_midfiles_df)


let depca_wells_from_notebook = CSV.read(well_list_from_notebook_path, DataFrame)
    @assert depca_wells_from_notebook == sort!(DataFrame(
        (; w.strain, ds = w.DS, w.experiment, w.well) for w in eachrow(wells_with_midfiles_df)))
end
##

ex_name_map = mapreduce(vcat, cond2exs) do (cond, exs)
    map(enumerate(exs)) do (i, ex)
        #ex => "$(condname(cond)), experiment $i"
        ex => "$(condname(cond)), $i"
    end
end |> OrderedDict

CSV.write("$out_path_mapping/ex-name-mapping.csv", [(; experiment=ex, sanitized=s) for (ex,s) in ex_name_map])
##

# map each well to a unique index within its experiment
well_index_map = mapreduce(vcat, cond2exs) do (cond, exs)
    cond_wells = cond2wells[cond]
    mapreduce(vcat, exs) do experiment
        wells = filter( well -> well.experiment == experiment, cond_wells )
        [well => i for (i, well) in enumerate(wells)]
    end
end |> OrderedDict

CSV.write("$out_path_mapping/well-indexes.csv", [(; w.experiment, w.well, index) for (w, index) in well_index_map])

##

well_mids_filename(cond, ex, wellname; contour_method) = joinpath(midpoints_path, condname(cond), 
    Elegans.midpoints_filename(ex, wellname; midpoints_path="", contour_method)
)


# Comment out desired version:
# For debugging
_cp(source, dest) = (@assert isfile(source); println(dest))
#_cp(source, dest) = println("CP: \n\t  $source \n\t⇒ $dest")
# # For actual run
# _cp(source, dest) = isfile(dest) || cp(source, dest)



individual_depca_filenames = readdir(individual_depca_path)

@progress "conditions" for cond in conditions
    #cond_path = "$midpoints_path/$(condname(cond))"
    wells = cond2wells[cond]

    cname = condname(cond)
    cname_full = condname_full(cond)
    println()
    println(cname)
    println()

    cond_coord_out_path = joinpath(out_path, "coords-and-size", cname)
    cond_mids_out_path = joinpath(out_path, "midpoints", cname)
    cond_depcas_out_path = joinpath(out_path, "individual-DE-PCAs", cname)
    mkpath(cond_coord_out_path)
    mkpath(cond_mids_out_path)
    mkpath(cond_depcas_out_path)

    # let well = first(wells)
    @progress "wells" for well in wells
        well_name = "$(well.experiment)-$(well.well)"

        ex_out = ex_name_map[well.experiment]
        well_i = well_index_map[well]

        println("$well_name => $ex_out-$well_i")

        possible_coord_paths = [joinpath(c, well.experiment, well.well, "coords_and_size.jld2") for c in coordinate_paths]
        coord_file_path = possible_coord_paths[only(findall(isfile, possible_coord_paths))]
        coord_file_out_path = joinpath(cond_coord_out_path, "$ex_out-$(well_i) coords-and-size.jld2")
        
       
        _cp(coord_file_path, coord_file_out_path)

        # `Elegans.midpoints_filename` assumes midpoint files are directly in `midpoints_path`.
        # Here, we instead expect them do be organized by in directories by condition.
        for contour_method in contour_methods
            mids_filename = Elegans.midpoints_filename(well.experiment, well.well; midpoints_path="", contour_method)
            mids_path = joinpath(midpoints_path, condname(cond), mids_filename)
            @assert isfile(mids_path)

            mids_out_filename = Elegans.midpoints_filename(ex_out, string(well_i); midpoints_path="", contour_method)
            mids_out_path = joinpath(cond_mids_out_path, mids_out_filename)
            
            _cp(mids_path, mids_out_path)
        end

        depca_filename = only(filter(startswith("DE-PCA $(well.experiment) $(well.well), $cname"), individual_depca_filenames))
        depca_path = joinpath(individual_depca_path, depca_filename)
        #depca_out_path = joinpath(cond_depcas_out_path, "DE-PCA $ex_out-$well_i.jld2")
        r = r"DE-PCA " * "$(well.experiment) $(well.well), $cname_full, " * r"(.*)"
        m = match(r, depca_filename)
        m === nothing && error("filename `$depca_filename`` does not match $r")
        depca_out_path = joinpath(cond_depcas_out_path, "DE-PCA $ex_out, well $well_i, $(only(m.captures))")

        _cp(depca_path, depca_out_path)
    end
end    


## Open population DE-PCA files and save versions with sanitized names

[cond => length(wells) for (cond, wells) in cond2wells]

##

# Suffix expected on all DE-PCA files, based on params that were used in the pipepline 
depca_params_suffix = "n_trim=1, windowlength=30, conf_th=0.05"
@assert Set(conditions) == keys(cond2wells)
depca_filename(cond) = "DE-PCA $(condname_full(cond)), n=$(length(cond2wells[cond])), $depca_params_suffix.jld2"
depca_file_path(cond) = joinpath(depca_path, depca_filename(cond))

for cond in conditions
    @assert isfile(depca_file_path(cond))
end


##

depca_out_path = "$out_path/DE-PCAs"
mkpath(depca_out_path)

@progress for cond in conditions
    @info condname(cond)

    path = depca_file_path(cond)
    dest = "$depca_out_path/$(depca_filename(cond))"
    isfile(dest) && continue

    @info "Loading `$path`"
    d = load(path)

    (; strain, DS) = cond
    @assert (strain, DS) == only(d["conditions"])

    de_indices, wells, pcas = d["de_indices"], d["wells"], d["pcas"]
    @assert issorted(wells)
    @assert wells == sort(cond2wells[(;strain, DS)])
    # In the sanitized data, store sanitized experiment names and a count of wells in each experiment,
    # instead of well list
    experiments = unique(w.experiment for w in wells)
    n_individuals = [ex_name_map[ex] => count(w.experiment == ex for w in wells) for ex in experiments]

    @assert sum( n for (_,n) in n_individuals ) == length(wells) == size( d["de_indices"], 2 )

    @info "Saving to `$dest`"
    jldsave(dest; strain, DS, n_individuals, de_indices, pcas)
end

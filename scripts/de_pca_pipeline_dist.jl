using Pkg

local_project_root = normpath("$(@__DIR__)/..")

Pkg.activate(local_project_root)
@info "Installing packages..."
Pkg.instantiate()

##
using TOML

args = TOML.parsefile("$(@__DIR__)/args.toml")
paths, dist = args["paths"], args["dist"]

ex_list_file = normpath(joinpath(local_project_root, paths["ex_list_file"]))
stage_path   = normpath(joinpath(local_project_root, paths["stage_path"]))

local_root = paths["local_root"]
logdir = paths["logdir"]

remote_root = dist["remote_root"]
remote_project_root = "$remote_root/$(dist["remote_project_root"])"

midpoints_relpath = paths["midpoints_relpath"]
des_relpath       = paths["des_relpath"]
out_relpath       = paths["out_relpath"]

local_midpoints_path  = "$local_root/$midpoints_relpath"
remote_midpoints_path = "$remote_root/$midpoints_relpath"

assert_isdir(path)  = isdir(path)  || error("Not found or not a directory: $path")
assert_isfile(path) = isfile(path) || error("Not found or not a regular file: $path")

assert_isdir(local_root)
assert_isdir(local_project_root)
assert_isdir(local_midpoints_path)
assert_isfile(ex_list_file)
assert_isfile(stage_path)

##
using Distributed

remote = dist["remote"]
n_local, n_remote = dist["n_local"], dist["n_remote"]

@info "Starting workers..." remote remote_root remote_project_root local_root local_project_root n_local n_remote logdir

rmprocs(setdiff(workers(), [1])...) # remove any existing workers, to allow running this script twice in the same Julia session
local_procs = cd(local_project_root) do
    addprocs(n_local; exeflags = "--project")
end
remote_procs = n_remote == 0 ?
                Int[] : # `addprocs` has a bug where where it adds 1 proc when asked for 0
                addprocs([(remote, n_remote)];
                      dir = "$remote_project_root", 
                      exename = "$remote_root/julia-$VERSION/bin/julia",
                      exeflags = "--project"
                      )

##

@info "Loading packages..."

# @everywhere using Pkg

# println(join(
#     map(procs()) do i
#         s = remotecall_fetch(()->sprint(io->Pkg.status(;io)), i)
#         "Worker $i:\n$s"
#     end, 
#     "\n"))

# ##

using Elegans
##


@everywhere begin
    using Elegans
    using ElegansTimeSeries
    using ElegansTimeSeries: depca_filename, depca_vars_filename, depca_cv_filename, depca_individual_filename
    using CSV, DataFrames
    using JLD2
    using OffsetArrays
    using Unzip
    using ProgressLogging
    using BlockArrays
    using MultivariateStats
    using OffsetArrays: no_offset_view
end

##

contour_methods = Dict( 1 => Thresholding(1.0,0.35), 
        ( i => Thresholding(1.0,0.34) for i=2:5 )... )

stages = sort!(collect(keys(contour_methods)))

##

# `stringtype=String` avoids having CSV's custom string types in the resulting 
# .jld2 files (which would require loading CSV.jl to open).
ex_df = DataFrame(CSV.File(ex_list_file; stringtype = String))

@info "$(nrow(ex_df)) experiments listed in $ex_list_file..."

ex_dict = Dict(strain => Dict(ds => 
                        ex_df.name[(ex_df.DS .== ds) .& (ex_df.strain .== strain)]
                    for ds in unique(ex_df.DS))
                for strain in unique(ex_df.strain))

stagedict = loadstages(stage_path)


conditions = sort!(collect(Iterators.flatten(((str, cond) for cond in keys(ex_dict[str])) for str in keys(ex_dict))))

exs = [ex_dict[str][cond] for (str,cond) in conditions]

## Keep experiments with stage data

for cond_exs in exs
    filter!( ex -> haskey(stagedict,ex), cond_exs )
end
exs

##

exs_all = reduce(vcat, exs)

@info "Processing $(length(exs_all)) experiments that have stage data..."


##

well_ids = [ 
    [(; experiment, well) for experiment in cond_exs for well in sort!(collect(keys(stagedict[experiment])))]
    for cond_exs in exs]

@info "... $(sum(length, well_ids)) wells"

## Use only wells with midpoints files

headtail_method = Elegans.SpeedHTCM(5,0)
end_assignment_params = Elegans.EndAssignmentParams()

foreach(well_ids) do v
    filter!(v) do w
        all(isfile(Elegans.midpoints_filename(w.experiment, w.well; midpoints_path = local_midpoints_path, 
                    contour_method, headtail_method, end_assignment_params))
                for contour_method in unique(values(contour_methods)))
    end
end
well_ids

@info "Processing $(sum(length, well_ids)) wells with midpoints..."

##

@everywhere begin
    function well_bin_mid_angles( experiment, wellname, bin, nbins; 
                                midpoints_cache = midpoints_cache, stagedict )
        # assumes the midpoint_cache was created with `full = true`
        (; mids, conf, iters) = midpoints_cache(experiment, wellname)
        bin_range, stage_i, _ = well_bin_range(experiment, wellname, bin, nbins; stagedict )

        angles = Elegans.mids2turns(mids[stage_i][bin_range,:])
        conf  = conf[stage_i][bin_range]
        iters = iters[stage_i][bin_range]
        (; angles, conf, iters)
    end

    node_trim(mat; ntrim) = mat[:, begin+ntrim:end-ntrim]
    allfinite(v) = all(isfinite,v)

    function well_bin_de_mat( experiment, wellname, bin, nbins, winlen; 
            midpoints_cache = midpoints_cache,
            stagedict, ntrim = ntrim, cond=allfinite )
            
        (; angles, conf, iters) = well_bin_mid_angles( experiment, wellname, bin, nbins; 
                midpoints_cache, stagedict )

        vecs, i = ElegansTimeSeries.delay_embed(node_trim(angles; ntrim)', winlen; cond)
        # `i` indexes into possible slices, starting with index 1. 
        # `conf` and `iters` are offset, indexed by frame number.
        # de-offseting them and indexing by `i` gives the the conf/iters at the first frame of each DE window in `vecs`.
        conf, iters = no_offset_view(conf)[i], no_offset_view(iters)[i]
        (; vecs, i, conf, iters)
    end

    function well_bin_range( experiment, wellname, bin, nbins; stagedict )
        stage_i, bin_in_stage = fldmod1(bin, nbins)

        stage_range = stage_frames(experiment, wellname, stage_i; stagedict)
        
        stage_bin_ranges = ElegansTimeSeries.bins(stage_range, nbins)
        bin_range = stage_bin_ranges[bin_in_stage]
        bin_range, stage_i, bin_in_stage
        #@info "DEPCA for $bin_range"
    end
end

##

params = args["params"]
nbins, pca_winlen, ntrim = params["nbins"], params["pca_winlen"], params["ntrim"]

confidence_threshold = params["confidence_threshold"]

maxoutdim = params["maxoutdim"] # Applied only to individuals. Condition PCAs are computed to full dimensionality

cv_pca_nsamples, cv_pca_nwindows = params["cv_pca_nsamples"], params["cv_pca_nwindows"]

@show nbins pca_winlen ntrim confidence_threshold maxoutdim cv_pca_nsamples cv_pca_nwindows

##

@everywhere include("de_pca_pipeline_funcs.jl")

@everywhere function load_or_compute_des(condition, wells, outpath; nbins, pca_winlen, ntrim, midpoints_cache, stagedict)
    condname = "$(condition[1]) $(condition[2])DS"
    fname = "$outpath/DEs with conf $condname (nbins=$nbins, winlen=$pca_winlen, ntrim=$ntrim).jld2"
    return if isfile(fname)
        @info "... loading from file `$fname`"
        d = load(fname)
        d["wells"] == wells || error("""
                Well list from file does not match input list.
                Expected: $wells
                In file: $(d["wells"])
                """)
        d["pop_des"], d["pop_de_indices"], d["conf"], d["iters"]
    else
        @info "... computing delay embeddings... `$fname`"
        nbins_total = nbins * length(Elegans.stage_names)
        # wells are in the outer loop and bins are in the inner loop, 
        # so that each well's midlines are loaded from the cache if necessary 
        # in an iteration of the outer loop, ensuring it's still in the cache
        # for later bins, and giving a better progress bar percentage
        @progress "Delay embeddings" res = [well_bin_de_mat(w.experiment, w.well, 
                                    bin, nbins, pca_winlen; 
                                    ntrim,
                                    midpoints_cache,
                                    stagedict
        ) for 
            bin in 1:nbins_total, 
            w in wells
        ]
        des, indices, conf, iters = unzip(res)
        @info "... saving to `$fname`"
        mktemp(dirname(fname)) do path, io
            @info "...      saving to `$path`"
            jldsave(path; pop_des = des, pop_de_indices = indices, conf, iters, wells)
            @info "...      rename `$path` => `$fname`"
            mv(path, fname)
        end
        #jldsave(fname; pop_des = des, pop_de_indices = indices, conf, iters, wells)
        des, indices, conf, iters
    end
end


@everywhere using StatsBase: sample

##

# TODO copied from Parallelism.jl

function retry_check(delay_state, err)
    # Below each condition to retry is listed along with an explanation about why
    # retrying should/might work.
    should_retry = (
        # Worker death is normally stocastic, if not then doesn't matter how many
        # retries as it will rapidly kill all workers
        err isa ProcessExitedException ||
        # If we are in the middle of fetching data and the process is killed we could
        # get an ArgumentError saying that the stream was closed or unusable.
        # So same as above.
        err isa ArgumentError && occursin("stream is closed or unusable", err.msg) ||
        # In general IO errors can be transient and related to network blips
        err isa Base.IOError
    )
    if should_retry
        @info "Retrying computation that failed due to a $(typeof(err)): $err"
    else
        @warn "Non-retryable $(typeof(err)) occurred: $err"
    end
    return should_retry
end



##


using WebIO, CSSUtil, Mux
@everywhere using ObservablePmap
@everywhere using LoggingExtras: TeeLogger
@everywhere using TerminalLoggers
@everywhere using Dates

##

@everywhere using ElegansTimeSeries: progress_throttle_logger, logger_with_timestamps

# Do `f` with a log file per worker
function with_logfiles(f, base_path)
    
    logfile_streams = Dict{Int,IOStream}()
    
    # `base_path` must be relative, as it is used as-is on both local and remote workers
    @assert !isabspath(base_path)
    run_id = Dates.format(now(), "yyyy-mm-dd HH.MM.SS")
    path = "$base_path/$run_id"

    function logger_f(io)
        # A terminal logger for the io that `ologpmap` will read from
        terminal_logger = TerminalLogger(io)
        # Also log to a file with timestamps
        id = myid()
        mkpath(path)
        logfile_streams[id] = open("$path/worker-$id.log"; append=true)
        file_logger = progress_throttle_logger(Second(1), logger_with_timestamps(logfile_streams[id]))
        TeeLogger(terminal_logger, file_logger)
    end

    try
        f(logger_f)
    finally
        for id in keys(logfile_streams)
            (s = get(logfile_streams, id, nothing)) === nothing || close(s)
        end
    end
end

##        

num_retries = 3

summ, task = with_logfiles(logdir) do logger_f 
    ologpmap(eachindex(conditions); logger_f, on_error = identity, 
            retry_check, retry_delays = ExponentialBackOff(n=num_retries)) do cond_i

        condition = conditions[cond_i]
        wells = well_ids[cond_i]

        isremote = myid() in remote_procs
        root = isremote ? remote_root : local_root
        midpoints_path = isremote ? remote_midpoints_path : local_midpoints_path

        # We can create a separate midpoint cache in each condition, since they don't share anything
        midpoints_cache = Elegans.well2midpoints_cache(contour_methods; midpoints_path, full = true)

        strain, ds = condition
        condname = "$strain $(ds)DS"
        condname_short = ds > 0 ? condname : strain
        nwells = length(wells)
        nbins_total = nbins * length(Elegans.stage_names)

        pop_des_path = "$root/$des_relpath"
        out_path = "$root/$out_relpath"
        individual_out_path = "$out_path/individuals/$condname_short"
        mkpath(pop_des_path)
        mkpath(out_path)
        mkpath(individual_out_path)

        
        @info "$condname: $nwells wells"
        
        depca_params = (; nwells, ntrim, pca_winlen, confidence_threshold)
        outfile_pca  = joinpath(out_path, depca_filename(condname, depca_params))
        outfile_vars = joinpath(out_path, depca_vars_filename(condname, depca_params))
        outfile_cv   = joinpath(out_path, depca_cv_filename(condname, cv_pca_nwindows, cv_pca_nsamples, depca_params))
        outfiles_individual_pca = [joinpath(individual_out_path, depca_individual_filename(well, maxoutdim, depca_params))
                                for well in wells]
        outfiles = [[outfile_pca, outfile_vars, outfile_cv]; outfiles_individual_pca]

        # @show outfiles
        # @show filter(!isfile, outfiles)
        if all(isfile, outfiles)
            # Bail out early if all output files exist. No need to compute/load DEs.
            @info "Skipping $condname. All output files exists."
            #continue
            return
        end

        @info "Delay embeddings..."
        des, indices, conf, iters = nothing, nothing, nothing, nothing # For GC
        des, indices, conf, iters = load_or_compute_des(condition, wells, pop_des_path; nbins, pca_winlen, ntrim, midpoints_cache, stagedict)
        @assert (nbins_total, nwells) == size(des)

        # apply confidence threshold (for each bin and well)
        @info "... filtering by confidence"
        @progress "filtering" des = [
            let
                #@show i size(des[i])
                res = view(des[i], :, conf[i] .> confidence_threshold)
                #@show size(res)
                res
            end
            for i in eachindex(IndexCartesian(), des, conf)
        ]
        # TODO filter `indices` as well. Currently indices corresponds to `des` before confidence filtering
        
        if isfile(outfile_pca)
            @info "Loading PCAs for $condname from file `$outfile_pca`"
            depcas_multi = load(outfile_pca, "pcas")
        else
            @info "Computing PCAs..."
 
            # DEs for each bin, each individual in a block
            des_multi_for_pca = mortar.(reshape.(eachrow(des),1,:))
            #@time @progress depcas_multi = [Elegans.try_return(()->fit(PCA, Matrix(m))) for m in des_multi_for_pca[1:1]]
            depcas_multi = nothing # For GC
            @info "... fitting"
            @progress "DE PCAs" depcas_multi = [let
                    @info "DE PCAs bin $bin" size(m)
                    fit(PCA, Matrix(m)) 
                end 
                for (bin,m) in enumerate(des_multi_for_pca)]
            
            @info "Saving to `$outfile_pca`..."
            jldsave(outfile_pca; 
                    pcas = depcas_multi, 
                    conditions = [condition],
                    wells,
                    de_indices = indices
            )
        end

        @progress "Individuals" for i in eachindex(wells, outfiles_individual_pca)
            well, dest = wells[i], outfiles_individual_pca[i]
            if isfile(dest)
                @info "Skipping individual DE-PCA for $well. Output file already exists: $dest"
                continue
            end
            
            depcas = [
                let 
                    well_bin_des = des[bin,i]
                    if isempty(well_bin_des)
                        @warn "No valid DE windows" well bin
                        nothing
                    else
                        fit(PCA, Matrix(well_bin_des); maxoutdim)
                    end    
                end
                for bin in 1:nbins_total]
            jldsave(dest; depcas)
            GC.gc()
        end

        if isfile(outfile_vars)
            @info "Skipping score variances for $condname. Output file already exists: $outfile_vars"
        else
            @info "Computing PCA score variances..."
            depca_vars_multibin = let
                # iterate over wells in the outer loop to use the midpoints cache 
                # more efficiently, then permute dims 
                @progress "DEPCA variances" res = [ depca_vars( convert(Matrix, des[bin,well_i]), depcas_multi[bin] ) 
                    for bin in 1:nbins_total, well_i in 1:nwells
                ]
                vars, tvars, nwindows = unzip(permutedims(res))
                (; vars, tvars, nwindows)
            end

            jldsave(outfile_vars; 
                depca_vars_multibin.vars, depca_vars_multibin.tvars, depca_vars_multibin.nwindows,
                wells
            )
        end

        if isfile(outfile_cv)
            @info "Skipping CV errors for $condname. Output file already exists: $outfile_cv"
        else
            @info "Computing PCA cross-validation errors..."
            ndims = [1:50; 60:10:400]
			@progress "CV PCA" cvpca_errs = [
				let
                    # Sample `cv_pca_nwindows` windows in this bin from the whole population,
                    # and collect to a matrix
                    m = mortar(des[bin:bin,:])
                    size(m, 2) < cv_pca_nwindows && error("Not enough DE windows ($(size(m, 2))) in bin $bin for condition $condname")
					mat = Matrix(sampledim(m, cv_pca_nwindows, 2; replace = false))
					err = cv_pca_errs(mat, cv_pca_nsamples, ndims)
					(; ndims, err)
				end
				for bin in 1:nbins_total]
			
			@info "Saving CV PCA errors to `$outfile_cv`"
			jldsave(outfile_cv; cvpca_errs)
        end
    end
end


h = map(x -> style(HTML("<pre>$x</pre>")), summ)

port = rand(8100:8200)
WebIO.webio_serve(page("/", vbox(h)), port)
address = "localhost:$port"
println()
println("Progress available at $address (copied to clipboard)")
clipboard(address)

Threads.@spawn begin
    wait(task)
    @info "Done"
end

nothing
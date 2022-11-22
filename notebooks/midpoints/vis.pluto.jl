### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 70c6f550-633f-11ed-1910-436db7422453
begin
	using Pkg
	Pkg.activate(".")
	Pkg.status
	using Revise
	#using PlutoPlotly
	using PlutoUI
end

# ╔═╡ bd0ac901-02c1-4411-afa8-a13679750c59
using Elegans

# ╔═╡ abe38ea8-9b92-4376-baba-698aba1ae6e8
using ElegansTimeSeries

# ╔═╡ 71368a38-0e33-4e1e-b91a-c8f5bde4afe1
 begin
	 using Markdown
	 #using MarkdownLiteral: @markdown
	 using IterTools
	 using DataFrames
	 using MosaicViews
	 using ProgressLogging
	 using FileIO
	 using Images
	 using MultivariateStats
	 using Statistics
	 using LinearAlgebra
	 using LazyArrays
	 using ImageFiltering
	 using ImageFiltering.KernelFactors: gaussian
 end

# ╔═╡ 8ba09f88-2ad1-47f8-aa2a-4729a09e8214
begin
	using OffsetArrays
	using OffsetArrays: no_offset_view
end

# ╔═╡ 57b4b915-cb50-4206-8de6-994cb9a4239e
using Plots

# ╔═╡ 9c5efb27-cd27-4720-bba2-cbf109102f25
begin
	using GeometryBasics
	function postures_plot!(f!, plt, i, padsize = nothing; lw=2, la=0.2, c=:plasma, kwargs...)
	    plt_kwargs = Dict{Symbol,Any}(kwargs)
	    push!(plt_kwargs, :lw=>lw, :la=>la, :c=>c)
	    for k in i
	        f!( k, plt_kwargs )
	    end
	    plot!(plt, legend=false, aspect_ratio=1, bg="black", 
	            axis=false, ticks=false, yflip=true, cbar=false; plt_kwargs...)
	    plt
	end
	
	plot_contours!(plt, contours, traj, i; kwargs...) = postures_plot!(plt, i; kwargs...) do k, plt_kwargs
	        p = Point2(traj.x[k], traj.y[k])
	        c = contours(k)
	        if c isa AbstractVector && !isempty(c)
	            plot!(plt, Elegans.vertices_closed(first(c)) .+ p, linez=k; plt_kwargs...)
	        end
	end
	
	function plot_midlines!(plt, midpts, i, padsize; kwargs...)
		postures_plot!(plt, i; kwargs...) do k, plt_kwargs
		    all(ismissing, midpts[k,:]) && return
		    m = replace(midpts[k,:], missing=>Elegans.missingpoint)
		    plot!(plt, m, linez=k; plt_kwargs...)
		end
		zoom_around_mids_center!(plt, midpts[i,:], padsize)
	end
	
	plot_contours(contours, traj, i; kwargs...) = plot_contours!(plot(), contours, traj, i; kwargs...)
	plot_midlines(midpts, i, padsize = nothing; kwargs...) = plot_midlines!(plot(), midpts, i, padsize; kwargs...)

	function zoom_around_mids_center!(plt, midpts, padsize)
		center = mean(filter(p -> all(isfinite,p), midpts))
		if all(isfinite, center)
			xlims!(center[1] .+ (-padsize, padsize))
			ylims!(center[2] .+ (-padsize, padsize))
		end
		plt
	end
end

# ╔═╡ bf16d177-8291-4c7b-8137-b6b17956c82d
Pkg.status("Plots")

# ╔═╡ 7a59acf5-3d22-4190-9857-c7d9c435109b
# Avoid the very slow default `show` for `PCA` (MultivariateStats.jl issue #186)
function Base.show(io::IO, ::MIME"text/plain", M::PCA)
    idim, odim = size(M)
    print(io, "PCA(indim = $idim, outdim = $odim, principalratio = $(MultivariateStats.r2(M)))")
end

# ╔═╡ 92c912e5-ff20-4e5c-a546-c9debd64598d
video_path = Sys.iswindows() ? "U:/experiments/reemy" : "$(homedir())/experiments/reemy"

# ╔═╡ 64799c8d-c74d-4f9b-af9b-353fb0c11be8
begin
	contour_method = Thresholding(1.0, 0.34)
	contour_methods = Dict( 1 => Thresholding(1.0,0.35), 
        ( i => Thresholding(1.0,0.34) for i=2:5 )... )
	contour_path = Sys.iswindows() ? "U:/cached-data/contours" : "$(homedir())/cached-data/contours"
	midpoints_path = Sys.iswindows() ? "H:/cached-data/midpoints" : "$(homedir())/cached-data/midpoints"
	depca_path = Sys.iswindows() ? "U:/ElegansTimeSeries/notebooks/midpoints" : "$(homedir())/ElegansTimeSeries/notebooks/midpoints"
end

# ╔═╡ f14dd6e9-2f23-43a4-9b10-784176d15f20
const STAGE_NAMES = ["L1", "L2", "L3", "L4", "A"]

# ╔═╡ 78e19a2e-af78-4988-9f54-74ff71da1a12
stagedict = loadstages();

# ╔═╡ 7efed74f-d0e5-46fd-b81d-1b08b9fc6188
#exps = [name for name in readdir(video_path) if startswith(name, "RA")]
exps = sort!(intersect!(readdir(video_path), keys(stagedict)))

# ╔═╡ 6ef0573f-ae6b-46f9-a50a-666d166d6292
md"# Frame mosaic"

# ╔═╡ 76109068-6b32-44cf-991e-5df223402665
@bind vv_ex Select(exps)

# ╔═╡ c52254a8-dba8-48cc-97f9-4bc40c164e73
let vv_wells = readdir("$video_path/$vv_ex")
	@bind vv_wellname Select(vv_wells)
end

# ╔═╡ 868921ff-9c62-4e2d-9131-72313b845894
begin
	vv_well = Well(video_path, vv_ex, vv_wellname)
	vv_traj = load_coords_and_size(vv_well)
	vv_vcache = VideoCache(vv_well; traj = vv_traj, maxsize = 100);
end

# ╔═╡ e48c5b54-2872-45db-ae63-2b03d5480f55
begin
	vv_stage_ends = stagedict[vv_ex][vv_wellname]
	vv_stage_ranges = [from+1:to for (from, to) in IterTools.partition(vv_stage_ends,2,1)]
	# # extra space fixes weird @md_str behavior with initial interpolated string
	# md"""
	#  $(STAGE_NAMES[vv_stage_i]) frames: $(vv_stage_ranges[vv_stage_i])
	
	# $(@bind vv_frame_i NumberField(vv_stage_ranges[vv_stage_i]))
	# """
end

# ╔═╡ 350d771b-0151-40cb-939b-edfdc390e0d0
md"per stage: $(@bind per_stage NumberField(2:1000))"

# ╔═╡ 413f1d36-f19c-4f96-8a10-2d80abce6eef
i = reduce(hcat, map(vv_stage_ranges) do r
	round.(Int, range(first(r), last(r), length=per_stage+1)[1:end-1])
end)

# ╔═╡ 5ffc03f8-9123-484f-92de-1133dcd37aef
@progress frames = [get_frame(vv_vcache, i) for i in i'];

# ╔═╡ 6fda44a1-a58c-4d72-81e0-b70725f896a1
begin
	crange_sl = @bind crange RangeSlider(0:0.01:1; default=0.08:0.01:0.30)
	md"color range $crange_sl"
end

# ╔═╡ 82952df6-dcfe-4900-802f-74f61ec27964
#enhance(frame; a=a, b=b) = Gray.(a .* (gray.(Gray.(frame)) .+ b .- 1/2) .+ 1/2)
enhance(frame; a=first(crange), b=last(crange)) = Gray.((gray.(Gray.(frame) .- a)) ./ b)

# ╔═╡ 418730b0-d84c-4856-bc34-253ecffb9838
m = mosaicview(vec(enhance.(frames)); ncol = per_stage, npad=4)

# ╔═╡ 8b40128d-5bbc-4feb-9335-af5186ca8c92
# b = let io = IOBuffer()
# 	save(Stream{format"PNG"}(io), m)
# 	take!(io)
# end

# ╔═╡ ced42a59-2e05-44fc-8a0b-8689d72e13f9
DownloadButton(sprint(io->save(Stream{format"PNG"}(io), m)), "mosaic $vv_ex $vv_wellname ($per_stage per stage).png")

# ╔═╡ 282365b7-0274-48f9-9439-726b3190a12a
md"## midlines"

# ╔═╡ 249e29f6-c504-4520-a595-83221481012f
midpoints_cache = Elegans.well2midpoints_cache(contour_methods; midpoints_path)

# ╔═╡ 75bfd7c3-7092-444b-925e-c4a1eea49ecc
vv_mids = Elegans.load_well_midpoints(vv_well, contour_methods, vv_stage_ranges;
						    midpoints_path);

# ╔═╡ dc694c29-e14b-435d-8e1e-d1d1666c23c3
function traj_plot(i, j, pointvecs, labels, contours, shift, σ; kwd...)
    k = gaussian(σ)
    plt = plot(aspect_ratio=1, yflip=true, legend=false)
    for (ci,(v,l)) in enumerate(zip(pointvecs, labels))
        plot!(imfilter(v[j],k), label=label, c=ci; kwd...)
        scatter!([v[i]], c=ci)
    end
    cs = contours(i)
    cs isa AbstractVector && for c in cs
        plot!(vertices_closed(c) .+ [shift[i]], lw=3, la=0.25, c="black")
    end
    plt
end


# ╔═╡ 67a93bd7-a13f-4aac-95b6-e7a6972d8b54
begin
	# plot a frame at the correct x/y coordinates 
	function plot_frame!(plt::Plots.Plot, i, vcache, traj; 
						 frame_span = (-75:75, -75:75))
		fr = get_frame(vcache, i)
		plot!(plt, traj.x[i] .+ frame_span[1], traj.y[i] .+ frame_span[2], fr)
		plt
	end
	plot_frame!(plt::Plots.Plot, args...; kwargs...) = throw(MethodError(plot_frame!, (plt, args...)))
	plot_frame!(args...; kwargs...) = plot_frame!(current(), args...; kwargs...)
	plot_frame(args...; kwargs...) = plot_frame!(plot(), args...; kwargs...)

	function plot_midline!(plt::Plots.Plot, i, traj, mids, contours, stage_ends)
		stage = searchsortedfirst(stage_ends, i)-1
		midline = mids[stage][i, begin:end]
		plot!(plt, midline, label="", lw=2)
		scatter!(plt, midline[[begin]], ms=5, c = "white", label="head", alpha=0.5)
		scatter!(plt, midline[[end]], ms=5, c="black", label="tail", alpha=0.5)
		if contours !== nothing
			for c in contours[i]
				plot!(plt, (Point2(traj.x[i], traj.y[i]),) .+ c, label="")
			end
		end
		plt
	end
	plot_midline!(plt::Plots.Plot, args...; kwargs...) = throw(MethodError(plot_midline!, (plt, args...)))
	plot_midline!(args...; kwargs...) = plot_midline!(current(), args...; kwargs...)
end

# ╔═╡ 048717ad-2e77-44da-8b65-cf770cd87fb1
@bind winlen NumberField(1:100_000; default=180)

# ╔═╡ 8ffa6fbd-dc60-4eb2-aafb-c95b22774312
ranges = [k:k+winlen-1 for k in i]

# ╔═╡ 72e61581-4e4f-4479-b07a-82b1c5ff46a2
md"zoom $(@bind padsize Slider(1:200; default=60, show_value=true))"

# ╔═╡ c00120f5-8f97-4cec-9e4f-2740865daf9d
plts = reduce(hcat, [plot_midlines(vv_mids[j], r, padsize; foreground_color_border="white", axis=true) for r in ranges[:,j]] for j in 1:5)

# ╔═╡ 180d7688-a344-4059-8476-33e6ca491f2a
plot(plts..., layout = reverse(size(i)), size = 400 .* size(i),
	bg="black")

# ╔═╡ 70a21682-f6a8-4c09-9a4c-fd44d19420c5
# let i = 500_000
# 	plot_frame(i, vv_vcache, vv_traj)
# 	plot_midline!(i, vv_traj, vv_mids, nothing, vv_stage_ends)
# end

# ╔═╡ 25d3f210-4185-4660-9bcb-c47d561903f0
md"### video"

# ╔═╡ 30cf199f-c5cc-4ab9-8a9f-a546f528519e
md"color range $crange_sl"

# ╔═╡ f0ef0f7d-59a1-4fff-9568-26bd8fa6eeab
begin
	function timeline_plot(stage_ranges, t)
		plt = plot(legend=false, ticks=false, axis=false, title = Elegans.frametime_str(t))
		for stage in 1:5
			r = stage_ranges[stage]
			vspan!([first(r), last(r)+1], label=STAGE_NAMES[stage])
			vline!([first(r), last(r)+1], label=STAGE_NAMES[stage], c=Gray(0.2))
			annotate!([middle(r)], [0.5], text(STAGE_NAMES[stage], "white"))
		end
		vline!([t], c="black", lw=2)
		plt
	end
	function make_vid_frame(i, stage_ranges; vcache = vv_vcache)
		plot(plot(enhance(get_frame(vcache, i)),
					
					axis=false, ticks=false),
			timeline_plot(stage_ranges, i),
			layout = grid(2,1, heights=(0.9,0.1))
		)
	end
end

# ╔═╡ 9f4539be-0af9-488c-9ed7-efde0c2d5892
# let
# 	anim = @animate for i=1:5, j=0:0.1:i
# 		bar([5,i,j], legend=false)
# 	end
# 	gif(anim)
# end
#gif(anim, fps=3)

# ╔═╡ 76eba5a3-94b9-4e23-aa8c-9c1d5259bc03
md"""
frame skip: $(@bind dframes NumberField(1:100_000; default=1000))
fps: $(@bind fps NumberField(1:100; default=5))
"""

# ╔═╡ e402928b-71db-484c-a9b9-e3184a78b374
vid_frames = map(1:5) do stage
	r = vv_stage_ranges[stage]
	range(first(r), last(r), step=dframes)
end

# ╔═╡ 0fec7cb0-bbca-499a-acca-fda83bd1e06e
md"preview frame $(@bind preview_frame Slider(first(first(vid_frames)):last(last(vid_frames))))"

# ╔═╡ 119e192b-a453-4433-8dbe-96c9247ba389
make_vid_frame(preview_frame, vv_stage_ranges)

# ╔═╡ 488dd852-1951-4ada-872c-fec993a92a4c
begin
	last_vid_button_count = Ref(0)
	@bind vid_button_count CounterButton("Make video")
end

# ╔═╡ 01f55a1b-ad80-4c16-85c0-a0e3bd34edb4
vid = if vid_button_count > last_vid_button_count[]
	last_vid_button_count[] = vid_button_count
	mktempdir() do dir
		path = "$dir/$(vv_ex)-$(vv_wellname)-shapes.mp4"
		@info "Writing video to $path"
		anim = @withprogress name="animation" begin
			@animate for stage in 1:5, (i,j) in enumerate(vid_frames[stage])
				make_vid_frame(j, vv_stage_ranges)
				# plot(plot(enhance(get_frame(vv_vcache, j)), axis=false, ticks=false),
				# 	timeline_plot(vv_stage_ranges, j),
				# 	layout = grid(2,1, heights=(0.8,0.2))
				# )
				@logprogress (stage-1) / 5 + i / 5length(vid_frames[stage])
			end
		end
		mp4(anim, path; fps)
		read(path)
	end
else
	nothing
end;

# ╔═╡ f958f616-79c4-40fb-b09f-ef6633a46a4d
if vid === nothing
	md"(no video)"
else
	DownloadButton(vid, "$(vv_ex)-$(vv_wellname)-shapes.mp4")
end

# ╔═╡ a9c0edf4-d038-46ac-89b1-d0a4d7e18ecb
md"## angles"

# ╔═╡ 4f64c797-550f-4054-b122-6465a625612b
vv_all_mids = let stage_mids = values(sort(vv_mids))
	for (a,b) in IterTools.partition(stage_mids,2,1)
		@assert lastindex(a,1) + 1 == firstindex(b,1)
	end
	# ApplyArray(vcat, ...) directly of offset arrays is broken
	OffsetArray(ApplyArray(vcat, no_offset_view.(stage_mids)...),
		firstindex(first(stage_mids), 1) - 1, 0
	)
end;

# ╔═╡ 018548ea-b72d-469b-a9cb-d48875ccdebe
vv_angles = Elegans.mids2turns(vv_all_mids)

# ╔═╡ 2e618924-ba89-415e-8808-c9b6aa94b98f
md"""
 $(@bind angles_animate CheckBox()) animate angles 
 - start frame $(@bind angles_anim_start NumberField(1:1_000_000; default=300_001)) length $(@bind angles_anim_len NumberField(1:100_000; default=100))
 - heatmap window half length $(@bind angles_anim_hm_halflen NumberField(1:10_000; default=60))
 - fps $(@bind anim_angles_fps NumberField(1:1000; default=20))
"""

# ╔═╡ 8e59f48d-cdc5-441c-b8b3-0dd59ddb1c7f
md"## DE PCA"

# ╔═╡ e680a293-6af7-47b5-a16c-c383e107f37c
@bind refresh_file_list_count CounterButton("Refresh file list")

# ╔═╡ 3b967101-c45c-4970-ba53-7b711ceadf56
begin
	refresh_file_list_count
	md"""Load population PCAs from
	$(@bind depca_file Select(filter(endswith(".jld2"), readdir(depca_path))))
	"""
	# @bind pc_multi_source confirm(PlutoUI.combine() do Child
	# 	md"""
	# 	Load population PCAs from file $(Child("from_file", CheckBox(default=true)))
	# 	$(Child("filename", Select(filter(endswith(".jld2"), readdir(depca_path)))))
	# 	"""
	# 	end)
end

# ╔═╡ 9117977d-1c23-4b15-af02-562d8b785cda
md"""**TODO**
window length and trim size should be stored in the jld2 file of PCAs
"""

# ╔═╡ c3a4fa1d-cb48-4333-ae49-4d543295bc74
md"""
window length $(@bind depca_winlen NumberField(1:1000; default=30))
trim size $(@bind depca_n_trim NumberField(1:100; default=1))
"""

# ╔═╡ 9a4f2795-11a5-4cec-8a97-83c8dbbf0cf2
depcas_multi, de_indices_multi = let
	d = load(joinpath(depca_path, depca_file))
	pcas = d["pcas"]
	@info "Loaded $(length(pcas)) PCAs for $(join(join.(d["conditions"], "-"), ", ")) (n=$(length(d["wells"])))"
	pcas, d["de_indices"]
end

# ╔═╡ 63715dc6-dcf8-4d61-9a2e-e7d6224883da
md" $(@bind window_reconstruct CheckBox()) reconstruct from $(@bind window_n_pcs NumberField(0:1000; default=10)) PCs"

# ╔═╡ c50c6858-a046-4771-9382-ffa62c672d5b
md"""
 $(@bind window_animate CheckBox()) animate window
-  $(@bind animate_anchor_tail CheckBox(default=true)) anchor tail 
-  $(@bind animate_node_markers CheckBox()) node markers
-  $(@bind animate_gif CheckBox(default=true)) GIF
-  $(@bind animate_fps NumberField(1:100; default=20)) fps
"""

# ╔═╡ 4a788cec-6c81-4617-830d-233889ac41fb
md"""**TODO** 
 - ✓ fixed axes
 - original along reconstructions
 - add original frames/midlines to animation
"""

# ╔═╡ f428bdd9-4e7b-4597-9370-a17e9d3a573a
# TODO animation code duplicated with midpoints notebook.
# Changes here should be migrated to `midpoints`
begin
	angles2midline(ϕ) = [Point(real(z),imag(z)) for z in cumsum([0; cis.(cumsum([0;ϕ]))])]

	function animate_midline_from_angles(ϕ, anchor_tail; node_markers=false)
		maxϕ = maximum(abs, ϕ)
		anchor_tail && (ϕ = reverse!(.-ϕ; dims = 1))
		mids = angles2midline.(eachcol(ϕ))
		x_ext = extrema(p[1] for mid in mids for p in mid) 
		y_ext = extrema(p[2] for mid in mids for p in mid) 
		
		anim = @animate for fr in axes(ϕ,2)
			angles = ϕ[:,fr]
			#anchor_tail && (angles = reverse!(.-angles))
			#mid = angles2midline(angles)
			mid = mids[fr]
			#mid = angles2midline(pcmat[:,fr])
			plt1 = plot(mid, 
				 lw = 3,
				 m = node_markers,
				 ms = 3,
				 mz = [0; angles; 0],
				 lz = [0; angles; 0],
				 c = :vik, 
				 xlims = x_ext, ylims = y_ext,
				 #xlims = (-max_d, max_d),
				 #ylims = (-max_d, max_d),
				 #clims = expand_symmetric,
				 cbar = false,
				 clims = (-maxϕ, maxϕ),
				 ticks = false,
				 aspect_ratio=1, label="")
			scatter!(mid[[begin,end]], 
				ms=12, mc="white", ma=0.5,label="")
			annotate!([mid[begin][1],mid[end][1]],
					  [mid[begin][2],mid[end][2]], 
					  anchor_tail ? ["t", "h"] : ["h", "t"])
		
			#plot(plt1, plot(pcmat[:,fr]), layout=(2,1))
		end
	end
end

# ╔═╡ 36ed3c44-993e-4cd3-87aa-418fa03ffbb1
# reinterpret a column of a DE PCA axis as a stat matrix
pcvec2mat(vec, winlen) = reshape(vec,:,winlen)

# ╔═╡ ccd33b32-e0f3-424d-b4ed-b4b7dfb2b6f6
# TODO
# These are copied from the `midpoints` notebook. Need to find a home
begin
	function expand_symmetric(a)
		a_fin = filter(isfinite, a)
		m = isempty(a_fin) ? -Inf : maximum(abs, extrema(a_fin))
		(-m,m) 
	end
	
	function bins(ax, nbins)
		edges = round.(Int, range(first(ax), last(ax)+1, length=nbins+1))
		[edges[i]:edges[i+1]-1 for i in 1:nbins]
	end
	
	binname(bin, nbins) = let (stage, bin_in_stage) = fldmod1(bin,nbins)
		"$(STAGE_NAMES[stage]): $bin_in_stage"
	end
	binname_long(bin, nbins=nbins) = binname(bin,nbins) * "/$nbins"
	allfinite(v) = all(isfinite, v)
end

# ╔═╡ b85c41fb-1665-445f-9042-016fa8b560fe
if angles_animate
	anim = @withprogress @animate for i in angles_anim_start .+ (0:angles_anim_len-1)
		pad = angles_anim_hm_halflen
		k = i-pad:i+pad
		nnodes = size(vv_angles,2)
		nodes = 2:nnodes-1
		hm = heatmap(k, nodes, vv_angles[k, nodes]', 
			#title = Elegans.frametime_str(i),
			c = :vik, clims=expand_symmetric, cbar = false, xticks=false)
		yticks!([2, nnodes-1], ["head", "tail"])
		vline!([i], label="", c="black")
		mid = vv_all_mids[i,:]
		mid_plt = plot(aspect_ratio=1, legend=:outerright, ticks=false, axis=false)
		plot_midline!(mid_plt, i, vv_traj, vv_mids, nothing, vv_stage_ends)
		# mid_plt = plot(mid, aspect_ratio=1, legend=false, lw=3, ticks=false, axis=false)
		p0 = mean(mid)
		r = 1.1 * maximum(norm.(mid .- p0))
		xlims!(p0[1] .+ (-r,r))
		ylims!(p0[2] .+ (-r,r))
	
		θ = vec(Elegans.mids2turns(reshape(mid,1,:)))
		θf = imfilter(θ[nodes] ./ π, gaussian(1))
		θ_plt = plot(nodes, θf,
				ylims=(-0.25,0.25), ylabel = "curvature", legend=false, yticks=false)
		scatter!([nodes[1]], [θf[1]], c="white", alpha=0.5)
		scatter!([nodes[end]], [θf[end]], c="black", alpha=0.5)
		hline!([0], c="black")
		xticks!([1, nnodes], ["head", "tail"])

		plot(hm, 
			plot(mid_plt, θ_plt), 
			timeline_plot(vv_stage_ranges, i),
			layout = grid(3,1, heights=(0.4,0.4,0.2)))
		
		@logprogress (i - angles_anim_start) / angles_anim_len
	end
	# mktempdir() do dir
	# 	path = "$dir/$(vv_ex)-$(vv_wellname)-angles.mp4"
	# 	@info "Writing video to $path"
	# 	mp4(anim, path; fps=anim_angles_fps)
	# 	#read(path)
	# end

	mp4(anim)
end

# ╔═╡ 7e1d4845-1513-446a-9cd1-270cf4544ce9
begin
	depca_nbins_total = length(depcas_multi)
	depca_nbins = Int(depca_nbins_total / length(STAGE_NAMES))
	bin_ranges = reduce(vcat, bins.(vv_stage_ranges, depca_nbins))
end;

# ╔═╡ 98d714f4-d765-4ec9-80d5-9fb56824df73
let r = 1:depca_nbins_total
	@bind bin Select(r .=> binname.(r, depca_nbins))
end

# ╔═╡ 09e1b99e-2f1e-406d-b796-85062e0926b1
begin
	bin_pca = depcas_multi[bin]
	bin_μ, bin_P = mean(bin_pca), projection(bin_pca)
	bin_range = bin_ranges[bin]
	#@bind win_center Slider(first(bin_range) + winlen÷2 : last(bin_range) - winlen÷2; show_value=true)
	md"bin frames: $bin_range"
end

# ╔═╡ 4b927f98-e456-4ad3-8d30-fa74683468cd
# Delay-embed the transpose to be consistent with how the PCs were computed
# (node varies faster than time within each DE vector)
bin_de, bin_de_indices = ElegansTimeSeries.delay_embed(vv_angles[bin_range, begin+depca_n_trim : end-depca_n_trim]', depca_winlen; cond=allfinite);

# ╔═╡ da709950-1641-421e-9152-cb521b81c621
md"valid window #: $(@bind window_i Slider(eachindex(bin_de_indices); show_value=true))"

# ╔═╡ caf054cb-ea9c-4389-a229-6f88be103456
begin
	window_start = bin_range[bin_de_indices[window_i]]
	md"window = $(window_start:window_start+winlen-1)"
end

# ╔═╡ 94a847d9-89c2-4ea2-8928-aaf85f5b28e7
window_de = bin_de[:,window_i];

# ╔═╡ c3f4aa1d-ed44-405b-ac0e-51b754448b57
angles_to_show = let
	angle_vec = if window_reconstruct
		y = MultivariateStats.transform(bin_pca, window_de)
		y[window_n_pcs+1:end] .= 0
		reconstruct(bin_pca, y)
	else
		window_de
	end
	pcvec2mat(angle_vec, depca_winlen)
end;

# ╔═╡ 9aa4cf37-f9ee-4d65-8594-ed2c43d1b5f7
begin
	bin_de_is_valid = falses(length(bin_range))
	bin_de_is_valid[bin_de_indices] .= true
	
	valid_wins_hm = heatmap(bin_range, 1:1, reshape(bin_de_is_valid,1,:), cbar=false, 		legend=false, xformatter=:plain, yticks=false, c = cgrad(["salmon", "lightgreen"]), axis=false
	)
	vline!([window_start], ls=:dash, c = "black", lw = 2)
	plot!(xformatter=:plain) # Plots.jl issue #4538
	
	win_hm = heatmap(angles_to_show; c = :vik, clims=expand_symmetric)

	window_de_t = MultivariateStats.transform(bin_pca, window_de)
	t_plot = bar(window_de_t[1:window_n_pcs], legend=false, ticks=false, axis=false)
	hline!([0], c="black")
	
	plot(valid_wins_hm, win_hm, t_plot,
		layout=grid(3, 1, heights=(0.1,0.8, 0.1)))
end

# ╔═╡ 68b891be-6cf0-45f1-a766-ccd262168437
if window_animate
	let fps = animate_fps, node_markers = animate_node_markers
		anim = animate_midline_from_angles(angles_to_show, animate_anchor_tail; 
			node_markers)
		animate_gif ? gif(anim; fps) : mp4(anim; fps)
	end
else
	nothing
end

# ╔═╡ Cell order:
# ╟─70c6f550-633f-11ed-1910-436db7422453
# ╠═bf16d177-8291-4c7b-8137-b6b17956c82d
# ╠═bd0ac901-02c1-4411-afa8-a13679750c59
# ╠═abe38ea8-9b92-4376-baba-698aba1ae6e8
# ╠═71368a38-0e33-4e1e-b91a-c8f5bde4afe1
# ╠═8ba09f88-2ad1-47f8-aa2a-4729a09e8214
# ╠═7a59acf5-3d22-4190-9857-c7d9c435109b
# ╟─57b4b915-cb50-4206-8de6-994cb9a4239e
# ╟─92c912e5-ff20-4e5c-a546-c9debd64598d
# ╟─64799c8d-c74d-4f9b-af9b-353fb0c11be8
# ╟─7efed74f-d0e5-46fd-b81d-1b08b9fc6188
# ╟─f14dd6e9-2f23-43a4-9b10-784176d15f20
# ╟─78e19a2e-af78-4988-9f54-74ff71da1a12
# ╟─6ef0573f-ae6b-46f9-a50a-666d166d6292
# ╟─76109068-6b32-44cf-991e-5df223402665
# ╟─c52254a8-dba8-48cc-97f9-4bc40c164e73
# ╟─868921ff-9c62-4e2d-9131-72313b845894
# ╟─e48c5b54-2872-45db-ae63-2b03d5480f55
# ╟─350d771b-0151-40cb-939b-edfdc390e0d0
# ╟─413f1d36-f19c-4f96-8a10-2d80abce6eef
# ╟─5ffc03f8-9123-484f-92de-1133dcd37aef
# ╟─6fda44a1-a58c-4d72-81e0-b70725f896a1
# ╟─82952df6-dcfe-4900-802f-74f61ec27964
# ╟─418730b0-d84c-4856-bc34-253ecffb9838
# ╟─8b40128d-5bbc-4feb-9335-af5186ca8c92
# ╟─ced42a59-2e05-44fc-8a0b-8689d72e13f9
# ╟─282365b7-0274-48f9-9439-726b3190a12a
# ╟─249e29f6-c504-4520-a595-83221481012f
# ╟─75bfd7c3-7092-444b-925e-c4a1eea49ecc
# ╟─dc694c29-e14b-435d-8e1e-d1d1666c23c3
# ╟─67a93bd7-a13f-4aac-95b6-e7a6972d8b54
# ╟─9c5efb27-cd27-4720-bba2-cbf109102f25
# ╟─048717ad-2e77-44da-8b65-cf770cd87fb1
# ╟─8ffa6fbd-dc60-4eb2-aafb-c95b22774312
# ╟─72e61581-4e4f-4479-b07a-82b1c5ff46a2
# ╟─c00120f5-8f97-4cec-9e4f-2740865daf9d
# ╠═180d7688-a344-4059-8476-33e6ca491f2a
# ╟─70a21682-f6a8-4c09-9a4c-fd44d19420c5
# ╟─25d3f210-4185-4660-9bcb-c47d561903f0
# ╟─30cf199f-c5cc-4ab9-8a9f-a546f528519e
# ╟─0fec7cb0-bbca-499a-acca-fda83bd1e06e
# ╟─119e192b-a453-4433-8dbe-96c9247ba389
# ╟─f0ef0f7d-59a1-4fff-9568-26bd8fa6eeab
# ╟─e402928b-71db-484c-a9b9-e3184a78b374
# ╟─9f4539be-0af9-488c-9ed7-efde0c2d5892
# ╟─76eba5a3-94b9-4e23-aa8c-9c1d5259bc03
# ╟─488dd852-1951-4ada-872c-fec993a92a4c
# ╟─f958f616-79c4-40fb-b09f-ef6633a46a4d
# ╟─01f55a1b-ad80-4c16-85c0-a0e3bd34edb4
# ╟─a9c0edf4-d038-46ac-89b1-d0a4d7e18ecb
# ╟─4f64c797-550f-4054-b122-6465a625612b
# ╟─018548ea-b72d-469b-a9cb-d48875ccdebe
# ╟─2e618924-ba89-415e-8808-c9b6aa94b98f
# ╟─b85c41fb-1665-445f-9042-016fa8b560fe
# ╟─8e59f48d-cdc5-441c-b8b3-0dd59ddb1c7f
# ╟─e680a293-6af7-47b5-a16c-c383e107f37c
# ╟─3b967101-c45c-4970-ba53-7b711ceadf56
# ╟─9117977d-1c23-4b15-af02-562d8b785cda
# ╟─c3a4fa1d-cb48-4333-ae49-4d543295bc74
# ╟─9a4f2795-11a5-4cec-8a97-83c8dbbf0cf2
# ╟─7e1d4845-1513-446a-9cd1-270cf4544ce9
# ╟─98d714f4-d765-4ec9-80d5-9fb56824df73
# ╟─09e1b99e-2f1e-406d-b796-85062e0926b1
# ╟─da709950-1641-421e-9152-cb521b81c621
# ╟─caf054cb-ea9c-4389-a229-6f88be103456
# ╟─9aa4cf37-f9ee-4d65-8594-ed2c43d1b5f7
# ╟─63715dc6-dcf8-4d61-9a2e-e7d6224883da
# ╟─c50c6858-a046-4771-9382-ffa62c672d5b
# ╟─68b891be-6cf0-45f1-a766-ccd262168437
# ╟─c3f4aa1d-ed44-405b-ac0e-51b754448b57
# ╟─4a788cec-6c81-4617-830d-233889ac41fb
# ╟─f428bdd9-4e7b-4597-9370-a17e9d3a573a
# ╟─94a847d9-89c2-4ea2-8928-aaf85f5b28e7
# ╟─4b927f98-e456-4ad3-8d30-fa74683468cd
# ╟─36ed3c44-993e-4cd3-87aa-418fa03ffbb1
# ╟─ccd33b32-e0f3-424d-b4ed-b4b7dfb2b6f6

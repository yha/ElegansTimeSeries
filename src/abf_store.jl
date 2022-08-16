using AlignedBinaryFormat
using OffsetArrays: no_offset_view
using GeometryBasics
using Unzip
using Elegans

write_arrays_abf(arrays, filename) = abfopen(filename, "w") do abf
    for (i,a) in pairs(arrays)
        write(abf, string(i), a)
    end
end
write_arrays_abf(filename; arrays...) = write_arrays_abf(arrays, filename)

function load_arrays_abf(filename)
    abf = abfopen(filename, "r")
    n = 0
    while haskey(abf, string(n+1))
        n += 1
    end
    [abf[string(i)] for i in 1:n], abf
end
function load_arrays_abf(f, filename)
    arrays, abf = load_arrays_abf(filename)
    try
        f(arrays)
    finally
        close(abf)
    end
end

# array of `Point2` to arrays of coordinates and back.
# also replace `missing` with `Point2(NaN,NaN)`
_m2nan_pt(pts::AbstractArray{<:Union{Missing,P}}) where {P<:Point} = replace(pts, missing=>P(NaN,NaN))
_m2nan_pt(pts::AbstractArray{<:Point}) = pts
function points2coords(pts)
    pts = _m2nan_pt(pts)
    reinterpret(reshape, eltype(eltype(pts)), pts)
end
coords2points(coords) = reinterpret(reshape, Point2{eltype(coords)}, coords)


# write midpoints to mem-mappable AlignedBinaryFormat
function write_mids_abf(well, path="mids"; 
            contour_methods, midpoints_path, headtail_method, end_assignment_params)
    fname = "$path/$(Elegans.wellname(well)).abf"
    if isfile(fname)
        @info "Skipping well: output file already exists: $fname"
        return fname
    end
    stages = 1:5
    #iranges = [stage_frames(well, i) for i ∈ stages]
    mids = Elegans.load_well_midpoints(well, contour_methods;#, iranges; 
                        midpoints_path, headtail_method, end_assignment_params)
    mids_vec = collect(values(sort(mids)))
    @assert sort!(collect(keys(mids))) == eachindex(mids_vec) == stages
    @info "Writing midpoints to ABF"
    @time abfopen(fname, "w") do abf
        for stage in stages
            #m = replace(mids_vec[stage], missing=>Elegans.GeometryBasics.Point2(NaN,NaN))
            #a = Array(no_offset_view(reinterpret(reshape, Float64, m)))
            a = points2coords(mids_vec[stage])
            write(abf, "$stage", Array(no_offset_view(a)))
            #frames = iranges[stage]
            frames = stage_frames(well, stage)
            #@show frames axes(a)
            @assert frames == axes(a,2)
            write(abf, "$stage-frames", [first(frames), last(frames)])
        end
    end
    fname
end

function load_mids_abf(filename)
    abf = abfopen(filename, "r")
    stages = 1:5
    pts, frames = unzip(map(stages) do stage
        a = abf[string(stage)]
        #pts = reinterpret(reshape, Point2{eltype(a)}, a)
        pts = coords2points(a)
        frames = abf["$stage-frames"]
        pts, frames
    end)
    pts, frames, abf
end
function load_mids_abf(f, filename)
    pts, frames, abf = load_mids_abf(filename)
    try
        f(pts, frames)
    finally
        close(abf)
    end
end

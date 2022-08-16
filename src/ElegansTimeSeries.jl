module ElegansTimeSeries

export points2coords, coords2points, write_arrays_abf, load_arrays_abf, write_mids_abf, load_mids_abf,
        bin_means_weights

include("abf_store.jl")
include("binning.jl")

import Elegans
import ImageFiltering.KernelFactors: gaussian

const ten_sec_filter = gaussian(15,121)
const ten_sec_d_filter = Elegans.velocitives_row_filter(ten_sec_filter)

end # module

module ImageMetrics

export
    abs2dif,
    absdif,
    autocrop,
    bounding_box,
    crop, crop!,
    interpolate, interpolate!,
    map_with_offsets, map_with_offsets!,
    resample,
    resample_slices,
    slice,
    slice_eltype,
    slice_ndims,
    slice_range,
    zerofill!

using TypeUtils, ArrayTools, EasyRanges, OffsetArrays, InterpolationKernels

include("utils.jl")
include("metrics.jl")
include("resample.jl")
include("interpolation.jl")
import .Interpolation: interpolate, interpolate! # NOTE: do not import `resample`

include("iic2024.jl")

end

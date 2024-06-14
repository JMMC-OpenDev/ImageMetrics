module ImageMetrics

export
    autocrop,
    bounding_box,
    crop, crop!,
    interpolate, interpolate!,
    map_with_offsets, map_with_offsets!,
    resample,
    zerofill!

using TypeUtils, ArrayTools, EasyRanges, OffsetArrays

include("utils.jl")
include("metrics.jl")
include("interpolation.jl")
import .Interpolation: resample, interpolate, interpolate!

end

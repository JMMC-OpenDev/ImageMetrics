module ImageMetrics

export
    autocrop,
    bounding_box,
    crop, crop!,
    interpolate, interpolate!,
    resample,
    zerofill!

using TypeUtils, ArrayTools, EasyRanges

include("utils.jl")
include("metrics.jl")
include("interpolation.jl")
import .Interpolation: resample, interpolate, interpolate!

end

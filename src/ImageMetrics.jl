module ImageMetrics

export
    interpolate, interpolate!,
    resample

using Unitless, TypeUtils, ArrayTools

include("metrics.jl")
include("interpolation.jl")
import .Interpolation: resample, interpolate, interpolate!

end

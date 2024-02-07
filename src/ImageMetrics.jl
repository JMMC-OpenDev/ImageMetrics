module ImageMetrics

export
    interpolate, interpolate!,
    resample

using Unitless, TypeUtils

include("interpolation.jl")
import .Interpolation: resample, interpolate, interpolate!

end

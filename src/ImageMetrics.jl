module ImageMetrics

export
    interpolate, interpolate!,
    resample

using TypeUtils, ArrayTools

include("utils.jl")
include("metrics.jl")
include("interpolation.jl")
import .Interpolation: resample, interpolate, interpolate!

end

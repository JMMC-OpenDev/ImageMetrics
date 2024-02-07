module ImageMetrics

export
    interpolate, interpolate!

# Write your package code here.
include("interpolation.jl")
import .Interpolation: interpolate, interpolate!

end

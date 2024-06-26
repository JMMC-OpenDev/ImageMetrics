"""

Module `ImageMetrics.InterferometricImagingContest2024` implements specific metrics methods  for
the 2024 edition of the *Interferometric Imaging Contest*.

"""
module InterferometricImagingContest2024

using ..ImageMetrics

using EasyFITS, TypeUtils, Unitful

struct ContestData{T,N,A<:AbstractArray{T,N},
                   P<:Number,
                   W<:AbstractVector{<:Number}}
    data::A
    pixelsize::P
    wave::W
end

function ContestData(filename::AbstractString)
    openfits(filename) do file
        hdu = file[1]
        @assert  hdu["NAXIS"].integer == 3
        axis1 = get_axis(hdu, 1)
        axis2 = get_axis(hdu, 2)
        @assert abs(step(axis1)) ≈ abs(step(axis2))
        pixelsize = abs(step(axis1))
        data = read(Array{Float64}, hdu)
        if step(axis1) > zero(step(axis1))
            data = data[end:-1:1, :, :]
            axis1 = reverse(axis1)
        end
        if step(axis2) < zero(step(axis2))
            data = data[:, end:-1:1, :, :]
            axis2 = reverse(axis2)
        end
        # Read wavelengths.
        ext = findfirst("WAVELENGTHS", file)
        if ext === nothing
            wave = get_axis(hdu, 3)
        else
            wave = read(Array{Float64}, file[ext], "WAVELENGTH") .* u"μm"
        end
        T = eltype(data)
        N = ndims(data)
        A = typeof(data)
        P = typeof(pixelsize)
        W = typeof(wave)
        return ContestData{T,N,A,P,W}(data, pixelsize, wave)
    end
end

"""
    ImageMetrics.InterferometricImagingContest2024.resample_cube(inp, pixscl, )

* `input_wave` and `output_wave` give the vectors of wavelengths in the spectral
  channels of the input and output cubes of images;

* `input_pixelsize` and `output_pixelsize` give the pixel size in the input and output
  cubes of images;

* `output_resolution` specifies the effective resolution in each spectral channel of the
  resulting image and expressed in the same units as `output_pixelsize`;

"""
function resample_cube(input::AbstractArray;
                       input_pixelsize::Number,
                       input_wave::AbstractVector{<:Number},
                       output_pixelsize::Number,
                       output_wave::AbstractVector{<:Number},
                       output_resolution::Union{Real,AbstractVector{<:Number}})

    if output_resolution isa AbstractVector && axes(output_resolution, 1) != axes(output_wave, 1)
        throw(DimensionMismatch("`output_resolution` and `output_wave` have different indices"))
    end

    # First, interpolate to the same spectral channels.
    output = resample_slices(input, input_wave => output_wave)

    # Second, resample each image to the ouput pixel size and resolution.
    F = floating_point_type(slice_eltype(output))
    rho = as(F, input_pixelsize/output_pixelsize)
    for k in eachindex(output)
        omega = as(F, get_value(output_resolution, k)/output_pixelsize)
        output[k] = resample(output[k]; ratio=rho, blur=omega)
    end
    return output
end

get_value(A::Number, i::Integer) = A
get_value(A::AbstractVector{<:Number}, i::Integer) = getindex(A, i)

get_integer(hdu::FitsImageHDU, key::AbstractString, def) =
    haskey(hdu, key) ? hdu[key].integer : def

get_float(hdu::FitsImageHDU, key::AbstractString, def) =
    haskey(hdu, key) ? hdu[key].float : def

get_string(hdu::FitsImageHDU, key::AbstractString, def) =
    haskey(hdu, key) ? hdu[key].string : def

function to_unit(str::AbstractString)
    str = lowercase(str)
    return str == "deg" ? u"°" :
        str == "rad" ? u"rad" :
        str == "micron" ? u"μm" :
        str == "MICRON" ? u"μm" :
        str == "meter" ? u"m" :
        str == "METER" ? u"m" :
        str == "arcsec" ? u"°"/3_600 :
        str == "mas" ? u"°"/3_600_000 :
        error("unknown units \"$str\"")
end

function get_axis(hdu::FitsImageHDU, i::Integer)
    naxis = hdu["NAXIS$i"].integer
    cdelt = get_float(hdu, "CDELT$i", 1.0)
    crpix = get_float(hdu, "CRPIX$i", 1.0)
    crval = get_float(hdu, "CRVAL$i", 0.0)
    cunit = get_string(hdu, "CUNIT$i", "")
    rng = ((1:naxis) .- crpix).*cdelt .+ crval
    if cunit == ""
        return rng
    else
        return rng.*to_unit(cunit)
    end
end

to_vector(::Type{T}, val::Number, len::Int) where {T} =
    fill!(Vector{T}(undef, len), as(T, val))

function to_vector(::Type{T}, A::AbstractVector, len::Int = length(A)) where {T}
    length(A) == len || throw(DimensionMismatch(
        "incompatible vector length, got $(length(A)), should be $(len)"))
    return convert(Vector{T}, A)
end

end # module

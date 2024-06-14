using Base: tail
using InterpolationKernels: support, infimum, supremum

const AIRY_FWHM = 1.02899397
const CUBICBSPLINE_FWHM = (4 + 8*cos((π + atan(3*sqrt(7)))/3))/3
const SINC_FWHM = 1.2067091288032283
const GAUSSIAN_FWHM = sqrt(log(256))

"""
    B = resample(A; blur, ratio)

resamples array `A` along its dimensions with a given amount of blur and a given sampling
ratio. The geometrical center of the result `B` coincides with the geometric center of
`A`. Keywords are:

* `blur` is the full-width at half-maximum (FWHM) of the blurring kernel (a cubic
  B-spline) expressed in units of the sampling step of the result `B`. However, if the
  absolute value of `blur` is less than 1, no blur is applied and `A` is interpolated
  using a Catmull-Rom cardinal cubic spline. This is not recommended for downsampling,
  that is when `abs(ratio) < 1`.

* `ratio` is the ratio of the sampling rate in `A` to that in `B`.

"""
function resample(A::AbstractArray; blur::Real=1, ratio::Real=1)
    # B[j] = Σ_j ϕ(((j - j0) - (i - i0)/ρ)*(fwhm(ϕ)/ω))*A[i]
    #      = Σ_j ϕ(αj*(j - j0) - αi*(i - i0))*A[i]
    T = floating_point_type(eltype(A))
    if abs(blur) ≥ one(blur)
        # Resample with a certain amount of blur.
        dst_stp = CUBICBSPLINE_FWHM/blur
        src_stp = dst_stp*ratio
        return resample(dst_stp, A, src_stp, BSpline{4,T}(), :)
    else
        return resample(inv(ratio), A, one(T), CatmullRomSpline{T}(), :)
    end
end

function resample(dst_stp::Real, src::AbstractArray,
                  src_stp::Real, ker::Kernel, ::Colon = Colon())
    return resample(dst_stp, src, src_stp, ker, ntuple(identity, Val(ndims(src))))
end

function resample(dst_stp::Real, src::AbstractArray,
                  src_stp::Real, ker::Kernel, dims::Tuple{Vararg{Integer}})
    return resample(dst_stp, src, src_stp, ker, map(Int, dims))
end

function resample(dst_stp::Real, src::AbstractArray,
                  src_stp::Real, ker::Kernel, dims::Tuple{Vararg{Int}})
    dst = resample(dst_stp, src, src_stp, ker, first(dims))
    return resample(dst_stp, dst, src_stp, ker, tail(dims))
end

function resample(dst_stp::Real, src::AbstractArray,
                  src_stp::Real, ker::Kernel, dims::Tuple{})
    return src # end the recursion
end

function resample(dst_stp::Real, src::AbstractArray,
                  src_stp::Real, ker::Kernel, d::Integer)
    1 ≤ d ≤ ndims(src) || error("invalid dimension to resample")
    return resample(dst_stp, src, src_stp, ker, Val(as(Int, d)))
end

function resample(dst_stp::Real, src::AbstractArray{<:Any,N},
                  src_stp::Real, ker::Kernel{T}, ::Val{d}) where {T,N,d}
    src_dims = size(src)
    dst_dim = resampled_length(dst_stp, axes(src, d), src_stp, ker)
    dst_dims = ntuple(i -> i == d ? dst_dim : src_dims[i], Val(N))
    dst = Array{float(eltype(src)),N}(undef, dst_dims)
    return resample!(dst, reference_index(dst, d), dst_stp,
                     src, reference_index(src, d), src_stp, ker, Val(d))
end

function resampled_length(dst_stp::Real, src_axis::AbstractUnitRange,
                          src_stp::Real, ker::Kernel)
    # The resampling formula writes:
    #
    #     dst[i] = sum_j ker((i - dst_ref)*dst_stp - (j - src_ref)*src_stp)*src[j]
    #
    # hence the destination index `i` is:
    #
    #     i = (k + (j - src_ref)*src_stp)/dst_stp + dst_ref
    #
    # for all indices `k` of the resampling kernel and indices `j` of the source.
    #
    Ω = support(ker)
    ja, jb = minmax(minimum(src_axis)*(src_stp/dst_stp),
                    maximum(src_axis)*(src_stp/dst_stp))
    ka, kb = minmax(infimum(Ω)/dst_stp,
                    supremum(Ω)/dst_stp)
    return floor(Int, jb + kb) - ceil(Int, ja + ka) + 1
end

function resampled_range(dst_ref::Real, dst_stp::Real, src_axis::AbstractUnitRange,
                         stp_ref::Real, src_stp::Real, ker::Kernel)
    # The resampling formula writes:
    #
    #     dst[i] = sum_j ker((i - dst_ref)*dst_stp - (j - src_ref)*src_stp)*src[j]
    #
    # hence the destination index `i` is:
    #
    #     i = (k + (j - src_ref)*src_stp)/dst_stp + dst_ref
    #
    # for all indices `k` of the resampling kernel and indices `j` of the source.
    #
    ja, jb = minmax((minimum(src_axis) - src_ref)*(src_stp/dst_stp),
                    (maximum(src_axis) - src_ref)*(src_stp/dst_stp))
    Ω = support(ker)
    ka, kb = minmax(infimum(Ω)/dst_stp,
                    supremum(Ω)/dst_stp)
    return ceil(Int, ja + ka + dst_ref):floor(Int, jb + kb + dst_ref)
end

"""
    ImageMetrics.resample!(dst, [dst_ref,] dst_stp,
                           src, [src_ref,] src_stp, ker, Val(d)) -> dst

overwrites destination array `dst` with source array `src` resampled along `d`-th
dimension with kernel `ker`. All axes of `src` and `dst` but the `d`-th one must be the
same.

Along the dimension of interpolation, the result is given by:

    dst[i] = sum_j ker((i - dst_ref)*dst_stp - (j - src_ref)*src_stp)*src[j]
           = sum_j ker((off - j)*stp)*src[j]

where `dst_ref` and `src_ref` are the fractional indices in the destination and source
arrays of any reference point coinciding in the two arrays (the center of their respective
index range by default) and with:

    off = src_ref + (i - dst_ref)*(dst_stp/src_stp)
    stp = src_stp

The range for `j` in the sum is the intersection of the support of the kernel (scaled and
shifted) and of the source index range:

    (off - j)*stp ∈ support(ker) <=> j ∈ [a + off, b + off] ∩ ℤ

with:

    a, b = minmax(infimum(support(ker))/stp,supremum(support(ker))/stp)

"""
function resample!(dst::AbstractArray{<:Any,N}, dst_stp::Real,
                   src::AbstractArray{<:Any,N}, src_stp::Real,
                   ker::Kernel, ::Val{d}) where {N,d}
    return resample!(dst, reference_index(dst, d), dst_stp,
                     src, reference_index(src, d), src_stp, ker, Val(d))
end

function resample!(dst::AbstractArray{<:Any,N}, dst_ref::Real, dst_stp::Real,
                   src::AbstractArray{<:Any,N}, src_ref::Real, src_stp::Real,
                   ker::Kernel{T}, ::Val{d}) where {T,N,d}
    dst_axes = axes(dst)
    src_axes = axes(src)
    for i in 1:N
        i == d || dst_axes[i] == src_axes[i] || throw(DimensionMismatch(
            "source and destination arrays have different indices along dimension $i"))
    end
    Ipre = CartesianIndices(dst_axes[1:d-1])
    Ipost = CartesianIndices(dst_axes[d+1:end])
    I = dst_axes[d]
    J = src_axes[d]
    # The resampling formula can rewritten as:
    #
    #     dst[i] = sum_j ker((i - dst_ref)*dst_stp - (j - src_ref)*src_stp)*src[j]
    #            = sum_j ker((off - j)*stp)*src[j]
    #
    # with:
    #
    #     off = src_ref + (i - dst_ref)*(dst_stp/src_stp)
    #     stp = src_stp
    #
    # The range for `j` in the sum is the intersection of the support of the kernel (scaled and
    # shifted) and of the source index range:
    #
    #     (off - j)*stp ∈ support(ker) <=> j ∈ [a + off, b + off] ∩ ℤ
    #
    # with:
    #
    #     a, b = minmax(infimum(support(ker))/stp, supremum(support(ker))/stp)
    #
    dst_ref = as(T, dst_ref)
    src_ref = as(T, src_ref)
    stp = as(T, src_stp)
    rho = as(T, dst_stp/src_stp) # magnification
    Ω = support(ker)
    a, b = minmax(as(T, infimum(Ω)/stp), as(T, supremum(Ω)/stp))
    zerofill!(dst)
    @inbounds for ipost in Ipost # assume column-major order
        for i in I
            # Offset and local interval for the discrete convolution.
            off = src_ref + rho*(i - dst_ref)
            Js = UnitRange(max(first(J), ceil(Int, a + off)),
                           min(last(J), floor(Int, b + off)))
            for ipre in Ipre
                s = zero(T)
                for j in Js
                    s += ker(stp*(off - j))*src[ipre,j,ipost]
                end
                dst[ipre,i,ipost] = s
            end
        end
    end
    return dst
end

reference_index(A::AbstractArray, d::Integer) = reference_index(axes(A, d))
reference_index(axis::AbstractUnitRange) = (first(axis) + last(axis))/2

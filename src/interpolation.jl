module Interpolation

export interpolate, interpolate!, resample

using TypeUtils, Unitless
import Base: axes1

# A colon means no interpolation, otherwise a vector of coordinates.
const InterpolationCoordinates = Union{Colon,AbstractVector{<:Real}}

"""
    interpolate(A, x...) -> B

interpolates array `A` along its dimensions at coordinates `x...`. The 1st of
`x...` is for the 1st dimension of `A`, the 2nd of `x...` is for the 2nd
dimension of `A`, and so on. Each of `x...` may be a vector of fractional
coordinates along the corresponding dimension of `A` or a colon `:` to not
interpolate this dimension. The trailing colons may be omitted.

"""
function interpolate(A::AbstractArray, xs::InterpolationCoordinates...)
    # Interpolation is separable, start with the 1st dimension.
    return interpolate(A, Val(1), xs...)
end

function interpolate(A::AbstractArray, ::Val{d}, x::InterpolationCoordinates,
                     xs::InterpolationCoordinates...) where {d}
    # Interpolate along d-th dimension, then proceed with the next dimension.
    return interpolate(interpolate(A, Val(d), x), Val(d+1), xs...)
end

function interpolate(A::AbstractArray, ::Val{d}) where {d}
    # All dimensions have been interpolated, return the result.
    return A
end

@inline function interp_coefs(x::T, I::AbstractUnitRange{Int}) where {T<:AbstractFloat}
    if x < first(I)
        i1 = first(I)
        i2 = i1
        w1 = zero(T)
        w2 = one(T)
    elseif x ≥ last(I)
        i1 = last(I)
        i2 = i1
        w1 = one(T)
        w2 = zero(T)
    else
        i1 = floor(Int, x)
        i2 = i1 + 1
        w2 = x - floor(x)
        w1 = one(T) - w2
    end
    return i1, w1, i2, w2
end

const out_of_range_dimension = BoundsError("out of bounds dimension of interpolation")

"""
    interpolate(A, Val(d), x) -> B

interpolates array `A` at coordinates `x` along dimension `d`. Coordinates `x`
are fractional indices along the `d`-th dimension of `A`. Continuous boundary
conditions are assumed to extrapolate `A`.

Coordinates `x` may also be a colon `:` to skip interpolating `d`-th dimension:

    interpolate(A, Val(d), :) -> A

See also [`interpolate!`](@ref).

"""
function interpolate(A::AbstractArray{<:Any,N}, ::Val{d},
                     x::InterpolationCoordinates) where {N,d}
    1 ≤ d ≤ N || throw(out_of_range_dimension)
    if x isa Colon
        # No interpolation along this dimension.
        return A
    else
        T = float(eltype(A))
        A_inds = axes(A)
        B_inds = ntuple(i -> (i == d ? axes1(x) : A_inds[i]), Val(N))
        B = similar(A, T, B_inds)
        return interpolate!(B, A, Val(d), x)
    end
end

"""
    interpolate!(B, A, Val(d), x) -> B

overwrites `B` with the interpolation of array `A` at coordinates `x` along
dimension `d`. Coordinates `x` are fractional indices along the `d`-th
dimension of `A`. Continuous boundary conditions are assumed to extrapolate
`A`. The leading (before `d`-th dimension) and trailing (after `d`-th
dimension) indices of `A` and `B` must be the same. The indices of `B` along
its `d`-th dimension must be identical to the indices in `x`.

See also [`interpolate`](@ref).

"""
function interpolate!(B::AbstractArray{<:Any,N},
                      A::AbstractArray{<:Any,N},
                      ::Val{d},
                      x::Colon) where {N,d}
    1 ≤ d ≤ N || throw(out_of_range_dimension)
    if B !== A
        axes(B) == axes(A) || throw(DimensionMismatch("axes are not the same"))
        copyto!(B, A)
    end
    return B
end

function interpolate!(B::AbstractArray{<:Any,N},
                      A::AbstractArray{<:Any,N},
                      ::Val{d},
                      x::AbstractVector{<:Real}) where {N,d}
    1 ≤ d ≤ N || throw(out_of_range_dimension)
    A_inds = axes(A)
    B_inds = axes(B)
    A_inds[1:d-1] == B_inds[1:d-1] || throw(DimensionMismatch(
        "leading axes are not the same"))
    A_inds[d+1:N] == B_inds[d+1:N] || throw(DimensionMismatch(
        "trailing axes are not the same"))
    axes1(x) == B_inds[d] || throw(DimensionMismatch(
        "vector of positions has incompatible indices"))

    I = CartesianIndices(B_inds[1:d-1])
    Jb = B_inds[d] # indices along interpolated dimension in B
    Ja = A_inds[d] # indices along interpolated dimension in A
    K = CartesianIndices(B_inds[d+1:N])

    # Assume A and B are in column-major order.
    @inbounds for k ∈ K
        for j ∈ Jb
            j1, w1, j2, w2 = interp_coefs(x[j], Ja)
            for i ∈ I
                B[i,j,k] = w1*A[i,j1,k] + w2*A[i,j2,k]
            end
        end
    end
    #= FIXME This optimization has no impact...
    if d == 1
        @inbounds for k ∈ K
            @simd for j ∈ Jb
                j1, w1, j2, w2 = interp_coefs(x[j], Ja)
                B[j,k] = w1*A[j1,k] + w2*A[j2,k]
            end
        end
    else
        @inbounds for k ∈ K
            for j ∈ Jb
                j1, w1, j2, w2 = interp_coefs(x[j], Ja)
                @simd for i ∈ I
                    B[i,j,k] = w1*A[i,j1,k] + w2*A[i,j2,k]
                end
            end
        end
    end
    =#
    return B
end

"""
    resample(A, dims...; kwds...) -> B

resamples array `A` to size `dims...` as if:

    B[i] ≈ A[(i - δ)/ρ]

with `δ` a translation and `ρ` a magnification factor.

Keywords are:

* `magnification = ρ` is the ratio of the pixel size of `B` divided by the
  pixel size of `A` (default 1).

* `translation = δ` is a tuple of translation Cartesian coordinates (no
  translation by default) for the input array `A`.

The geometrical center of `A`, shifted by `δ`, is mapped to the geometrical
center of `B`

"""
resample(A::AbstractArray, dims::Integer...; kwds...) =
    resample(A, dims; kwds...)
resample(A::AbstractArray, dims::Tuple{Vararg{Integer}}; kwds...) =
    resample(A, as(Dims, dims); kwds...)
function resample(A::AbstractArray{<:Any,N}, dims::Dims{N};
                  magnification::Real = 1,
                  translation::NTuple{N,Real} = ntuple(Returns(0), Val(N))) where {N}
    # FIXME: Base.has_offset_axes(A) && throw(error("array must have 1-based indices"))
    T = floating_point_type(eltype(A))
    ρ = as(T, magnification)
    J = axes(A)
    x = ntuple(d -> resampling_coordinates(
        T, Base.OneTo(dims[d]), J[d], as(T, translation[d]), ρ), Val(N))
    return interpolate(A, x...)
end

"""
    resampling_coordinates(T=Float64, I, J, δ, ρ) -> J′

yields a range `J′` of fractional indices such that (1-D) interpolation of a
vector `A` writes:

    B ≈ A[J′]

with `I` the index range of `B`, `J` the index range of `A`, `δ` a translation,
and `ρ` a magnification factor. To explain the roles of `δ` and `ρ`, the
interpolation is equivalent to:

    B[i] ≈ A[(i - δ)/ρ]       (∀ i ∈ I)

Optional argument `T` is to specify the floating-point type of the fractioanl
indices.

"""
function resampling_coordinates(::Type{T},
                                I::AbstractUnitRange{Int},
                                J::AbstractUnitRange{Int},
                                δ::Real,
                                ρ::Real) where {T<:AbstractFloat}
    # Let `i0` be the central index in `I`:
    #
    i0 = as(T, first(I) + last(I))/2
    #
    # and `j0` be the central index in `J`:
    #
    j0 = as(T, first(J) + last(J))/2
    #
    # Then, the fractional index is:
    #
    #     j = j0 + (i - i0 - δ)/ρ
    #       = i/ρ + (j0 - (i0 + δ)/ρ)
    #
    ρ = as(T, ρ)
    δ = as(T, δ)
    return (I ./ ρ) .+ (j0 - (i0 + δ)/ρ)
end

end # module

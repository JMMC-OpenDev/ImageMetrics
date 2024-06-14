"""
    zerofill!(A) -> A

fills array `A` with zeros and returns `A`.

"""
zerofill!(A::AbstractArray) = fill!(A, zero(eltype(A)))

"""
    zeropad(A, args...) -> B

zero-pads array `A` to dimensions or axes `args...`. The geometric centers (in the same
sense as for `fftshift`) of `A` and `B` are coincident.

"""
zeropad(A::AbstractArray, args::Integer...) = zeropad(A, args)
zeropad(A::AbstractArray, args::AbstractUnitRange{<:Integer}...) = zeropad(A, args)
function zeropad(A::AbstractArray{T,N},
                 args::Union{NTuple{N,Integer},
                             NTuple{N,AbstractUnitRange{<:Integer}}}) where {T,N}
    return zeropad!(similar(A, args), A)
end

"""
    zeropad!(dst, src) -> dst

overwrites destination `dst` with source `src` with zero-padding. The geometric centers
(in the same sense as for `fftshift`) of `dst` and `src` are coincident.

"""
function zeropad!(dst::AbstractArray{<:Any,N}, src::AbstractArray{<:Any,N}) where {N}
    I = axes(src)
    J = axes(dst)
    for d in 1:N
        length(I[d]) ≤ length(J[d]) || error(
            "destination dimensions must not be smaller than source dimensions for zero-padding")
    end
    offs = ntuple(d -> central_axis_index(J[d]) - central_axis_index(I[d]), Val(N))
    pad = zero(eltype(dst))
    return paste!(dst, src; offs, pad)
end

"""
    crop(A, args...) -> B

crops array `A` to dimensions or axes `args...`. The geometric centers (in the same sense
as for `fftshift`) of `A` and `B` are coincident.

"""
crop(A::AbstractArray, args::Integer...) = crop(A, args)
crop(A::AbstractArray, args::AbstractUnitRange{<:Integer}...) = crop(A, args)
function crop(A::AbstractArray{T,N},
              args::Union{NTuple{N,Integer},
                          NTuple{N,AbstractUnitRange{<:Integer}}}) where {T,N}
    return crop!(similar(A, args), A)
end

"""
    crop!(dst, src) -> dst

overwrites destination `dst` with source `src` with cropping. The geometric centers (in
the same sense as for `fftshift`) of `dst` and `src` are coincident.

"""
function crop!(dst::AbstractArray{<:Any,N}, src::AbstractArray{<:Any,N}) where {N}
    I = axes(src)
    J = axes(dst)
    for d in 1:N
        length(I[d]) ≥ length(J[d]) || error(
            "destination dimensions must not be larger than source dimensions for cropping")
    end
    offs = ntuple(d -> central_axis_index(I[d]) - central_axis_index(J[d]), Val(N))
    return paste!(dst, src; offs)
end

"""
    paste!(dst, src; offs, pad) -> dst

paste values of source array `src` into destination array `dst` with offsets `offs` and
padding value `pad`. By default, offsets are all equal to zero and the padding value is
`zero(eltype(dst))`.

"""
function paste!(dst::AbstractArray{<:Any,N}, src::AbstractArray{<:Any,N};
                offs::NTuple{N,Int} = ntuple(Returns(0), Val(N)),
                pad = zero(eltype(dst))) where {N}

    # Determine indices in source and destination of common region:
    #   dst[i + off] = src[i]   ∀ i ∈ Ip
    #   dst[j] = src[j - off]   ∀ j ∈ Jp
    I = axes(src)
    J = axes(dst)
    Ip = clamp_axes(I, J, -, offs) # common indices in source
    Jp = clamp_axes(J, I, +, offs) # common indices in destination

    # Fill destination with padding value if needed.
    @inbounds for d in 1:N
        if first(Jp[d]) > first(J) || last(Jp[d]) < last(J[d])
            fill!(dst, pad)
            break
        end
    end

    # Copy region if not empty.
    if !iszero(prod(map(length, Ip)))
        dst[Jp...] = view(src, Ip...)
    end
    return dst
end

@inline function clamp_axis(I::AbstractUnitRange{Int},
                            J::AbstractUnitRange{Int},
                            pm::Union{typeof(+),typeof(-)},
                            off::Int)
    return max(first(I), pm(first(J), off)) : min(last(I), pm(last(J), off))
end

@inline function clamp_axes(I::NTuple{N,AbstractUnitRange{Int}},
                            J::NTuple{N,AbstractUnitRange{Int}},
                            pm::Union{typeof(+),typeof(-)},
                            offs::NTuple{N,Int}) where {N}
    return ntuple(d -> clamp_axis(I[d], J[d], pm, offs[d]), Val(N))
end

"""
    ImageMetrics.central_axis_index(dim) -> i
    ImageMetrics.central_axis_index(rng) -> i

yield the central index of an array axis of length `dim` or with index range `rng`
following the same conventions as `fftshift`.

"""
central_axis_index(dim::Integer) = div(as(Int, dim), 2) + 1
central_axis_index(rng::AbstractUnitRange{<:Integer}) =
    div(as(Int, first(rng)) + as(Int, last(rng)) + 1, 2)

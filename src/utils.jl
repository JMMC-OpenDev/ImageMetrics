"""
    map_with_offsets(f, A, B, inds...; A_pad=zero(eltype(A)), B_pad=zero(eltype(B))) -> C

yields an array `C` with axes `inds...` such that:

    C[k] =   sum_{j ∈ (I + k) ∩ J} f(A[j-k], B[j])
           + sum_{i ∈ I \\ (I ∩ (J - k))} f(A[i], B_pad)
           + sum_{j ∈ J \\ ((I + k) ∩ J)} f(A_pad, B[j])

where `I`, `J`, and `K` are the respective Cartesian indices of arrays `A`, `B`, and `C`
while `A_pad` and `B_pad` are the padding values assumed for `A[i]` and `B[j]` for
out-of-bound indices `i` and/or `j`.

The function `f` must be such that `iszero(f(A_pad, B_pad))` holds.

"""
function map_with_offsets(f,
                          A::AbstractArray{<:Any,N},
                          B::AbstractArray{<:Any,N},
                          inds::AbstractUnitRange{<:Integer}...; kwds...) where {N}
    return map_with_offsets(f, A, B, inds; kwds...)
end

function map_with_offsets(f,
                          A::AbstractArray{<:Any,N},
                          B::AbstractArray{<:Any,N},
                          inds::NTuple{N,AbstractUnitRange{<:Integer}}; kwds...) where {N}
    return map_with_offsets(f, A, B, map(to_axis, inds); kwds...)
end

function map_with_offsets(f,
                          A::AbstractArray{<:Any,N},
                          B::AbstractArray{<:Any,N},
                          inds::NTuple{N,AbstractUnitRange{Int}}; kwds...) where {N}
    T = float(promote_type(eltype(A), eltype(B)))
    C = OffsetArray(Array{T,N}(undef, map(length, inds)), map(first, inds))
    return map_with_offsets!(C, f, A, B; kwds...)
end

"""
    map_with_offsets!(C, f, A, B; A_pad=zero(eltype(A)), B_pad=zero(eltype(B))) -> C

overwrites array `C` with:

    C[k] =   sum_{j ∈ (I + k) ∩ J} f(A[j-k], B[j])
           + sum_{i ∈ I \\ (I ∩ (J - k))} f(A[i], B_pad)
           + sum_{j ∈ J \\ ((I + k) ∩ J)} f(A_pad, B[j])

where `I`, `J`, and `K` are the respective Cartesian indices of arrays `A`, `B`, and `C`
while `A_pad` and `B_pad` are the padding values assumed for `A[i]` and `B[j]` for
out-of-bound indices `i` and/or `j`.

The function `f` must be such that `iszero(f(A_pad, B_pad))` holds.

"""
function map_with_offsets!(C::AbstractArray{<:Any,N},
                           f,
                           A::AbstractArray{<:Any,N},
                           B::AbstractArray{<:Any,N};
                           A_pad = zero(eltype(A)),
                           B_pad = zero(eltype(B))) where {N}
    A_pad = as(eltype(A), A_pad)
    B_pad = as(eltype(B), B_pad)
    f_zero = f(A_pad, B_pad)
    iszero(f_zero) || throw(AssertionError(
        "`iszero(f(A_pad,B_pad))` does not hold with `A_pad = $(A_pad)` and `B_pad = $(B_pad)`"))
    I = bounding_box(x -> x != A_pad, A)
    J = bounding_box(x -> x != B_pad, B)
    K = CartesianIndices(C)
    A_out = f_zero
    @inbounds for i in I
        A_out += f(A[i], B_pad)
    end
    B_out = f_zero
    @inbounds for j in J
        B_out += f(A_pad, B[j])
    end
    f_out = A_out + B_out
    @inbounds for k in K
        f_in = f_zero
        for j in @range (I + k) ∩ J
            i = j - k
            A_i = A[i]
            B_j = B[j]
            f_in += f(A_i, B_j) - (f(A_i, B_pad) + f(A_pad, B_j))
        end
        C[k] = f_in + f_out
    end
    return C
end

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
    B = autocrop([f = !iszero,] A)

crops array `A` to keep the smallest Cartesian region containing all entries of `A` such
that `f(A[i])` is true.

"""
autocrop(A::AbstractArray) = A[bounding_box(A)]
autocrop(f, A::AbstractArray) = A[bounding_box(f, A)]

"""
    B = bounding_box([f = !iszero,] A)

yields the smallest Cartesian region containing all entries of `A` such that `f(A[i])` is
true.

"""
bounding_box(A::AbstractArray{Bool}) = bounding_box(identity, A)
bounding_box(A::AbstractArray) = bounding_box(!iszero, A)
function bounding_box(f, A::AbstractArray{<:Any,N}) where {N}
    Imin = CartesianIndex(ntuple(d -> typemax(Int), Val(N)))
    Imax = CartesianIndex(ntuple(d -> typemin(Int), Val(N)))
    @inbounds for I in CartesianIndices(A)
        if f(A[I])
            Imin = min(Imin, I)
            Imax = max(Imax, I)
        end
    end
    return Imin:Imax
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

"""
    ImageMetrics.to_axis(dim::Integer)
    ImageMetrics.to_axis(rng::AbstractUnitRange{<:Integer})

yield array axis (as an instance of `AbstractUnitRange{Int}`) corresponding to array
dimension `dim` or to array index range `rng`. In the former case, 1-based indices are
assumed.

"""
to_axis(dim::Integer) = Base.OneTo{Int}(dim)
to_axis(rng::AbstractUnitRange{Int}) = rng
to_axis(rng::AbstractUnitRange{<:Integer}) = convert_eltype(Int, rng)

"""
    ImageMetrics.to_dim(dim::Integer)
    ImageMetrics.to_dim(rng::AbstractUnitRange{<:Integer})

yield array dimension length (an `Int`) corresponding to array dimension `dim` or to array
index range `rng`.

"""
to_dim(dim::Integer) = as(Int, dim)
to_dim(rng::AbstractUnitRange{Int}) = as(Int, length(rng))
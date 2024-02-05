
"""
    interpolate(A, i...) -> B

interpolates array `A` along its dimensions at fractional indices `i...`.

"""
function interpolate(A::AbstractArray{T,N},
                     i::Vararg{N,AbstractVector{<:Real}}) where {T,N}
    return interpolate(A, Val(1), i...)
end

function interpolate(A::AbstractArray{T,N},
                     ::Val{d},
                     xd::AbstractVector{<:Real},
                     xs::AbstractVector{<:Real}...) where {T,N,d}
    B = Array{float(T),N}(undef, )
    unsafe_interpolate!(B, A, xd, Val(d))
    return interpolate(A, Val(d+1), xs...)
end
function interpolate(A::AbstractArray{T,N},
                     ::Val{N},
                     xd::AbstractVector{<:Real},
                     xs::AbstractVector{<:Real}...) where {T,N,d}
    B = Array{float(T),N}(undef, )
    unsafe_interpolate!(B, A, xd, d)
    return interpolate(A, Val(d+1), xs...)
end


function interpolate!(dst::AbstractArray{T,N},
                      src::AbstractArray{T,N},
                      grd::AbstractVector{<:Real},
                      dim::Integer) where {T,N}
    @assert 1 ≤ dim ≤ N
    I_dst = axes(dst)
    I_src = axes(src)
    I_head = I_src[1:dim-1]
    I_tail = I_src[dim+1:N]
    I_run_src = I_src[dim]
    I_run_dst = I_dst[dim]

    I_dst[1:dim-1] == I_head || throw(DimensionMismatch(
        "leading axes are not the same"))
    I_dst[dim+1:N] == I_tail || throw(DimensionMismatch(
        "trailing axes are not the same"))
    axes(grd) == (I_run_src,) || throw(DimensionMismatch(
        "vector of positions has incompatible indices"))

    unsafe_foo!(dst, CartesianIndices(I_head), I_run_dst, CartesianIndices(I_tail),
                src, I_run_src, grd)
    return dst
                  ,
    const dim1 = size(src, 1)
    dst = Array(Float64, length(xvals))
    k = 0
    for x in xvals
        (i0, i1, u0, u1) = interpcoefs(x, dim1)
        k += 1
        dst[k] = src[i0]*u0 + src[i1]*u1
    end
    return dst
end

@Inline function interp_coefs(x::T, I::AbstractUnitRange{Int}) where {T<:AbstractFloat}
    if x < first(I)
        i1 = first(I)
        i2 = i2
        w1 = zero(T)
        w2 = one(T)
    elseif x ≥ last(I)
        i1 = last(I)
        i2 = i1
        w1 = one(T)
        w0 = zero(T)
    else
        i1 = floor(Int, x)
        i2 = i1 + 1
        w2 = x - floor(x)
        w1 = one(T) - w2
    end
    return i1, w1, i2, w2
end

function unsafe_interpolate(dst::AbstractArray{<:Any,N},
                            A::AbstractArray{<:Any,N},
                            x::AbstractVector{<:Real},
                            ::Val{d}) where {N,d}
    inds = axes(A)
    I = CartesianIndices(inds[1:d-1])
    J = axes(dst)[d]
    K = CartesianIndices(inds[d+1:N])
    Jp = inds[d]

    zerofill!(dst)

    # Assume A and dst are in column-major order.
    for k ∈ K
        for j ∈ J
            j1, w1, j2, w2 = interp_coefs(x[j], Jp)
            for i ∈ I
                dst[i,j,k] = w1*A[i,j1,k] + w2*A[i,j2,k]
            end
        end
    end
    return dst
end


zerofill!(A::AbstractArray) = fill!(A, zero(eltype(A)))

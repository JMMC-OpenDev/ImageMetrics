"""
    ImageMetrics.GammaCorrection{T}(γ) -> Γ
    ImageMetrics.GammaCorrection(γ) -> Γ

build a functor object `Γ` such that `Γ(x)` yields `sign(x)*abs(x)^γ`.

"""
struct GammaCorrection{T<:Real}
    exponent::T
end
(f::GammaCorrection)(x::Real) = sign(x)*abs(x)^f.exponent

absdif(x, y) = abs(x - y)
abs2dif(x, y) = abs2(x - y)

"""
    ImageMetrics.distance(x, y; kwds...) -> dist

yields the image distance between a restored image `x` and a reference image `y` computed
as:

    dist = Σᵢ Ψ(Γ(α⋅x[i] + β), Γ(y[i]))

where the parameters `α` and `β` and functions `Γ` and `Ψ` are specified by keywords:

| Keyword         | Symbol | Default    | Description            |
|:----------------|:-------|:-----------|:-----------------------|
| `scale`         | `α`    | `1`        | Brightness scale       |
| `bias`          | `β`    | `0`        | Brightness bias        |
| `enhancement`   | `Γ`    | `identity` | Brightness enhancement |
| `cost`          | `Ψ`    | `absdif`   | Pixelwise cost         |

"""
function distance(x::AbstractArray,
                  y::AbstractArray;
                  enhancement::Function = identity,
                  cost::Function = abs,
                  scale::Real = one(real(typeof(x))),
                  bias::Real = zero(real(typeof(x))))
    @assert_same_axes x y
    α, β, Γ, Ψ = scale, bias, enhancement, cost
    dist = zero(Ψ(Γ(α*zero(eltype(x)) + β) - Γ(zero(eltype(y)))))
    if axes(x) == axes(y)
        @inbounds @simd for i in eachindex(x, y)
            dist += Ψ(Γ(α*x[i] + β), Γ(y[i]))
        end
    else
        X = CartesianIndices(x)
        Y = CartesianIndices(y)
        Γz = Γ(β)
        # Distance of background w.r.t. `x`.
        let dx = zero(dist)
            @inbounds @simd for i in X
                Γx = Γ(α*x[i] + β)
                dx += Ψ(Γx, Γz)
            end
            dist += dx
        end
        # Distance of background w.r.t. `y`.
        let dy = zero(dist)
            @inbounds @simd for i in Y
                Γy = Γ(y[i])
                dy += Ψ(Γz, Γy)
            end
            dist += dy
        end
        # Distance of `x` w.r.t. `y` in common part minus respective distance to background.
        let dxy = zero(dist)
            @inbounds @simd for i in @range X ∩ Y
                Γx = Γ(α*x[i] + β)
                Γy = Γ(y[i])
                dxy += Ψ(Γx, Γy) - (Ψ(Γx, Γz) + Ψ(Γz, Γy))
            end
            dist += dxy
        end
    end
    return dist
end

"""
    ImageMetrics.GammaCorrection{T}(γ) -> Γ
    ImageMetrics.GammaCorrection(γ) -> Γ

build a functor object `Γ` such that `Γ(x)` yields `sign(x)*abs(x)^γ`.

"""
struct GammaCorrection{T<:Real}
    exponent::T
end
(f::GammaCorrection)(x::Real) = sign(x)*abs(x)^f.exponent

"""
    ImageMetrics.distance(x, y; kwds...) -> dist

yields the image distance between a restored image `x` and a reference image
`y` computed as:

    dist = Σᵢ Ψ(Γ(α⋅x[i] + β) - Γ(y[i]))

where the parameters `α` and `β` and functions `Γ` and `Ψ` are specified by
keywords:

| Keyword         | Symbol | Default    | Description            |
|:----------------|:-------|:-----------|:-----------------------|
| `scale`         | `α`    | `1`        | Brightness scale       |
| `bias`          | `β`    | `0`        | Brightness bias        |
| `enhancement`   | `Γ`    | `identity` | Brightness enhancement |
| `cost`          | `Ψ`    | `abs`      | Pixelwise cost         |

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
    @inbounds @simd for i in eachindex(x, y)
        dist += Ψ(Γ(α*x[i] + β) - Γ(y[i]))
    end
    return dist
end

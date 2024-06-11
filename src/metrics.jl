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
    ImageMetrics.absdif(x, y) -> abs(x - y)

yields the absolute difference between `x` and `y`.

"""
absdif(x::Number, y::Number) = abs(x - y)

"""
    ImageMetrics.abs2dif(x, y) -> abs2(x - y)

yields the squared absolute difference between `x` and `y`.

"""
abs2dif(x::Number, y::Number) = abs2(x - y)

"""
    ImageMetrics.distance(x, y; kwds...) -> dist

yields the image distance between a restored image `x` and a reference image `y` computed
as:

    dist = Σᵢ Ψ(Γ(α⋅x[i] + β), Γ(y[i]))

for pixels `i` in the common region of `x` and `y` and where the parameters `α` and `β`
and functions `Γ` and `Ψ` are specified by keywords:

| Keyword         | Symbol | Default                 | Description                  |
|:----------------|:-------|:------------------------|:-----------------------------|
| `scale`         | `α`    | `one(real(eltype(x)))`  | Brightness scale             |
| `bias`          | `β`    | `zero(real(eltype(x)))` | Brightness bias              |
| `out`           | `η`    | `zero(real(eltype(y)))` | Value of out-of-field pixels |
| `enhancement`   | `Γ`    | `identity`              | Brightness enhancement       |
| `cost`          | `Ψ`    | `absdif`                | Pixelwise cost               |

"""
function distance(x::AbstractArray{Tx,N},
                  y::AbstractArray{Ty,N};
                  enhancement::Function = identity,
                  cost::Function = abs,
                  scale::Number = one(real(Tx)),
                  bias::Number = zero(real(Tx)),
                  out::Number = zero(real(Ty))) where {Tx,Ty,N}
    # Convert floating-point type of α, β, and η whitout changing their units if any.
    T = floating_point_type(Tx, Ty) # floating-point type for computations
    α = convert_floating_point_type(T, scale)
    β = convert_floating_point_type(T, bias)
    η = convert_floating_point_type(T, out)
    Γ = enhancement
    Ψ = cost
    Γz = Γ(η)
    dist = zero(Ψ(Γz, Γz))
    if axes(x) == axes(y)
        @inbounds @simd for i in eachindex(x, y)
            dist += Ψ(Γ(α*x[i] + β), Γ(y[i]))
        end
    else
        X = CartesianIndices(x)
        Y = CartesianIndices(y)
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
                Γx, Γy = promote(Γ(α*x[i] + β), Γ(y[i]))
                dxy += Ψ(Γx, Γy) - (Ψ(Γx, Γz) + Ψ(Γz, Γy))
            end
            dist += dxy
        end
    end
    return dist
end

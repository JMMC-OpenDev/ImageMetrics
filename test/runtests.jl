using ImageMetrics
using Test

@testset "ImageMetrics.jl" begin
    @testset "Interpolation" begin
        x = -4.3:0.1:5.3
        y = -1.3:0.1:8.3
        i = firstindex(x):0.3:lastindex(x)
        j = firstindex(y)-2.1:0.3:lastindex(y)+3.2
        # 1-D interpolation
        xp = @inferred interpolate(x, i)
        @test xp ≈ first(xp):step(x)*step(i):last(xp)

        # 2-D interpolation
        f(x::Real, y::Real) = sin(x*x + y)
        z = f.(x, y')
        @test z === @inferred interpolate(z, :, :)
        @test z === @inferred interpolate(z, :)
        @test z === @inferred interpolate(z)
        r1 = @inferred interpolate(z, i, :)
        @test r1 == @inferred interpolate(z, i)
        r2 = @inferred interpolate(r1, :, j)
        r3 = @inferred interpolate(z, i, j)
        @test r2 == r3
    end
end

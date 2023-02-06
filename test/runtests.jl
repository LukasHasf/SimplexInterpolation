using Test
include("../src/main.jl")

@testset "point_inside_simplex" begin
    x = 1.0
    y = 1.0
    simplex_x = [0.0 0.0 3.0]
    simplex_y = [0.0 3.0 0.0]
    # Inside
    @test point_inside_simplex(x, y, simplex_x, simplex_y) == true
    # Outside
    @test point_inside_simplex(3.0, 3.0, simplex_x, simplex_y) == false
    # On the edge
    @test point_inside_simplex(1.7, 0.0, simplex_x, simplex_y) == true
end

@testset "Compare interpolation methods" begin
    xs = rand(10)
    ys = rand(10)
    zs = rand(10)
    xi = rand(2)
    yi = rand(2)
    zi1 = simplex_interpolate(xs, ys, zs, xi, yi)
    itp = SimplexInterpolator(xs, ys, zs)
    zi2 = interpolate(itp, xi, yi)
    @test zi1 == zi2
    xi = rand(1)
    yi = rand(1)
    zi1 = simplex_interpolate(xs, ys, zs, xi, yi)
    zi2 = interpolate(itp, xi, yi)
    @test zi1 == zi2
end

@testset "Compare barycentric coordinate computation" begin
    x = 1.0
    y = 1.0
    xs = [0.0 0.0 3.0]
    ys = [0.0 3.0 0.0]
    λ1 = compute_barycentric_coordinates(x, y, xs, ys)
    λ2 = compute_generalized_barycentric_coordinates(x, y, vec(xs), vec(ys))
    @test all(λ1 .≈ λ2)
end
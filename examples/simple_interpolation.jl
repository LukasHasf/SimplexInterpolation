include("../src/SimplexInterpolation.jl")
using Plots

function exc()
    n=4
    xs = rand(n)
    ys = rand(n)
    xs = [0.8045237241795106, 0.4469802438876931, 0.6968852120465394, 0.2406133811248723]#, 0.6]
    ys = [0.8826124018484488, 0.3882666291837913, 0.2999087734670568, 0.6276019990899697]#, 0.48]
    zs = Float64.(collect(1:length(xs)))
    Ny = 200
    Nx = 200
    itp_simplex = NaturalNeighborInterpolation.NaturalNeighborsInterpolator(xs, ys, zs)
    plot_summary(Ny, Nx, itp_simplex, xs, ys, zs, "example_interpolation_simplex.png")
    itp_natural = NaturalNeighborInterpolation.NaturalNeighborsInterpolator(xs, ys, zs; mode=NaturalNeighborInterpolation.MODE_NATURAL)
    plot_summary(Ny, Nx, itp_natural, xs, ys, zs, "example_interpolation_natural.png") 
end

function plot_summary(Ny, Nx, itp, xs, ys, zs, filename)
    zi = zeros(Ny, Nx)
    for c in CartesianIndices(zi)
        zi[c.I...] = only(NaturalNeighborInterpolation.interpolate(itp, (c.I ./ (Ny, Nx))...))
    end
    heatmap(zi)
    scatter!(ys .* Ny, xs .* Nx; label="", zcolor=zs)
    plot_triang(itp.tri, Nx, Ny)  
    savefig(filename)
end

function plot_triang(triangulation, Nx, Ny)
    for i in axes(triangulation.simplices, 1)
        simplex_ind = triangulation.simplices[i, :]
        ys = triangulation.points[simplex_ind, 1] .* Nx
        xs = triangulation.points[simplex_ind, 2] .* Ny
        plot!([xs..., xs[1]], [ys..., ys[1]];label="")
    end    
end

exc()
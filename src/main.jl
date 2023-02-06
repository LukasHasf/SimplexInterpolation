using Delaunay

const MODE_NATURAL = "natural"
const MODE_SIMPLEX = "simplex"

"""    compute_generalized_barycentric_coordinates(x, y, xs::AbstractVector, ys::AbstractVector)

Compute the barycentric coordinates `[λ₁, λ₂, ..., λₖ]` of a point `(x,y)` with respect to a set of points
with x-coordinates `xs` and y-coordinates `ys`.
"""
function compute_generalized_barycentric_coordinates(x, y, xs::AbstractVector, ys::AbstractVector)
    A = [ones(length(xs))'; xs'; ys']
    b = [length(xs), x, y]
    λ = A\b
    return λ
end
"""    compute_barycentric_coordinates(x, y, simplex_x, simplex_y)

For a simplex with vertices at x-coordinates `simplex_x` and y-coordinates `simplex_y`,
calculate the barycentric coordinates `[λ1, λ2, λ3]` of a point `(x,y)`.
"""
function compute_barycentric_coordinates(x, y, simplex_x, simplex_y)
    x1, x2, x3 = simplex_x
    y1, y2, y3 = simplex_y
    λ1 = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / ((y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3))
    λ2 = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / ((y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3))
    λ3 = 1 - λ1 - λ2
    return λ1, λ2, λ3 
end

"""    norm(v)

Return the Euclidean norm of the vector `v`.
"""
function norm(v)
    return compute_distance(v, zero.(v))
end

"""    compute_distance(a, b)

Compute the Euclidean distance between two points `a` and `b`.
"""
function compute_distance(a, b)
    return sqrt(sum(abs2, a .- b))
end

function distance_point_line(a, b, c)
    bc = c .- b
    ab = b .- a
    line_length = compute_distance(b, c)
    dist = abs(bc[1] * ab[2] - ab[1] * bc[2]) / line_length
    return dist
end

function nearest_neighbor_interpolation(points, values, query_point)
    nearest_value = nothing
    nearest_distance = nothing
    for i in axes(points, 1)
        point = points[i, :]
        dist = compute_distance(point, query_point)
        if isnothing(nearest_distance) || dist < nearest_distance
            nearest_distance = dist
            nearest_value = values[i]
        end
    end
    return nearest_value
end

function nearest_simplex_interpolation(triangulation, values, query_point)
    nearest_simplex_index = -1
    nearest_distance = nothing
    nearest_simplex_x = nothing
    nearest_simplex_y = nothing
    nearest_simplex_vertex_indices = nothing
    nearest_edge = nothing
    for i in axes(triangulation.simplices, 1)
        simplex_indices = triangulation.simplices[i, :]
        simplex_xs = triangulation.points[simplex_indices, 1]
        simplex_ys = triangulation.points[simplex_indices, 2]
        simplex_xs = [simplex_xs..., simplex_xs[1]]
        simplex_ys = [simplex_ys..., simplex_ys[1]]
        dist = nothing
        for j in 1:(length(simplex_xs)-1)
            b = [simplex_xs[j], simplex_ys[j]]
            c = [simplex_xs[j+1], simplex_ys[j+1]]
            d = distance_point_line(query_point, b, c)
            if isnothing(dist) || d < dist
                dist = d
                nearest_edge = [b, c]
            end
        end
        if isnothing(nearest_distance) || dist < nearest_distance
            nearest_simplex_index = i
            nearest_distance = dist
            nearest_simplex_x = simplex_xs
            nearest_simplex_y = simplex_ys
            nearest_simplex_vertex_indices = simplex_indices
        end
    end
    # Project query point onto nearest edge
    edge = nearest_edge[2] .- nearest_edge[1]
    projected_point = nearest_edge[1] .+ sum((query_point .- nearest_edge[1]) .* (edge)) / norm(edge)^2 .* edge
    λ = compute_barycentric_coordinates(projected_point..., nearest_simplex_x, nearest_simplex_y)
    return sum(λ .* values[nearest_simplex_vertex_indices])
end

"""    point_inside_simplex(x, y, simplex_x, simplex_y)

Return whether the point `(x,y)` is contained in the simplex with
x-coordinates `simplex_x` and y-coordinates `simplex_y`.
"""
function point_inside_simplex(x, y, simplex_x, simplex_y)
    # Compute the barycentric coordinates of the query point
    λ1, λ2, λ3 = compute_barycentric_coordinates(x,y, simplex_x, simplex_y)
    # Check if the query point is inside the simplex (i.e. if all barycentric coordiantes are non-negative)
    if λ1 >= 0 && λ2 >= 0 && λ3 >= 0
        return true
    end
    return false
end

"""    find_containing_simplex(x, y, triangulation)

Find the simplex in `triangulation` that contains the point `(x,y)`.

Returns the index of the containing simplex.

If no simplex is found, returns `-1`.
"""
function find_containing_simplex(x, y, triangulation)
    # Intiialized index of simplex containing the query point
    simplex_index = -1
    # Loop over all simplices
    for i in axes(triangulation.simplices, 1)
        # Construct current simplex
        simplex_ind = triangulation.simplices[i, :]
        xs = triangulation.points[simplex_ind, 1]
        ys = triangulation.points[simplex_ind, 2]
        # Check if query point is inside current simplex
        if point_inside_simplex(x, y, xs, ys)
            simplex_index = i
        end
    end
    return simplex_index
end

function simplex_weights(triangulation, x, y, simplex_index)
    # Get neighbors of vertex
    simplex_vertex_indices = triangulation.simplices[simplex_index, :]
    xs = triangulation.points[simplex_vertex_indices, 1]
    ys = triangulation.points[simplex_vertex_indices, 2]
    # Weight the value by the barycentric coordinates
    λ = compute_generalized_barycentric_coordinates(x, y, xs, ys)
    weights = Dict{Int64, Float64}(simplex_vertex_indices .=> λ)
    return collect(keys(weights)), collect(values(weights))
end

"""    calc_circumcenter(xs, ys)
Calculate the circumcenter of a triangle with x-coordinates `xs` and y-coordinates `ys`.
Returns `[x_cc, y_cc]`
"""
function calc_circumcenter(xs, ys)
    # Formula from https://en.wikipedia.org/wiki/Circumscribed_circle
    D = 2*(xs[1]*(ys[2]-ys[3]) + xs[2]*(ys[3]-ys[1]) + xs[3]*(ys[1]-ys[2]))
    x = 1/D * ((xs[1]^2+ys[1]^2)*(ys[2]-ys[3]) + (xs[2]^2+ys[2]^2)*(ys[3]-ys[1]) + (xs[3]^2+ys[3]^2)*(ys[1]-ys[2]))
    y = 1/D * ((xs[1]^2+ys[1]^2)*(xs[3]-xs[2]) + (xs[2]^2+ys[2]^2)*(xs[1]-xs[3]) + (xs[3]^2+ys[3]^2)*(xs[2]-xs[1]))
    return [x,y]
end

function calc_circumradius(xs, ys)
    x1, x2, x3 = xs
    y1, y2, y3 = ys
    a = norm([x2, y2] .- [x1, y1])
    b = norm([x3, y3] .- [x2, y2])
    c = norm([x1, y1] .- [x3, y3])
    s = (a+b+c)/2
    return a*b*c/(4*sqrt(s*(s-a)*(s-b)*(s-c)))
end

function point_inside_circumcircle(x, y, xs, ys)
    circumcenter = calc_circumcenter(xs, ys)
    circumradius = calc_circumradius(xs, ys)
    distance = norm([x,y] .- circumcenter)
    return distance < circumradius
end

function compute_natural_neighbors(triangulation, x, y)
    natural_neighbor_indices = []
    for i in axes(triangulation.simplices, 1)
        vertex_indices = triangulation.simplices[i,:]
        xs = triangulation.points[vertex_indices, 1]
        ys = triangulation.points[vertex_indices, 2]
        if point_inside_circumcircle(x, y, xs, ys)
            push!(natural_neighbor_indices, vertex_indices...)
        end
    end
    return unique(natural_neighbor_indices)
end

function natural_weights(triangulation, x, y)
    # Get natural neighbors of query point
    vertex_neighbors_indices = compute_natural_neighbors(triangulation, x, y)
    xs = triangulation.points[vertex_neighbors_indices, 1]
    ys = triangulation.points[vertex_neighbors_indices, 2]
    # Weight the value by the barycentric coordinates
    λ = compute_generalized_barycentric_coordinates(x, y, xs, ys)
    weights = Dict{Int64, Float64}(vertex_neighbors_indices .=> λ)
    return collect(keys(weights)), collect(values(weights))
end

function simplex_interpolate(tri::Triangulation, z::Vector{T}, xi::Vector{T}, yi::Vector{T}; mode=MODE_SIMPLEX) where {T <: Number}
    global MODE_NATURAL, MODE_SIMPLEX
    # Array to hold interpolated values
    zi = zeros(length(xi))
    # Indices for the convex hull
    convex_hull_indices = unique(tri.convex_hull)
    # Loop over query points
    for i in eachindex(xi)
        # Add point to delaunay triangulation
         simplex_index = find_containing_simplex(xi[i], yi[i], tri)
        if simplex_index==-1
            zi[i] = nearest_neighbor_interpolation(tri.points[convex_hull_indices, :], z[convex_hull_indices], [xi[i], yi[i]])
            #zi[i] = nearest_simplex_interpolation(tri, z, [xi[i], yi[i]])
            continue
        end
        if mode==MODE_SIMPLEX
            natural_neighbors, weights = simplex_weights(tri, xi[i], yi[i], simplex_index)
        elseif mode==MODE_NATURAL
            #new_points = cat(tri.points, [xi[i] yi[i]]; dims=1)
            #new_tri = delaunay(new_points)
            natural_neighbors, weights = natural_weights(tri, xi[i], yi[i])
        end
        weights ./= sum(weights)
        # Perform interpolation using the weights
        zi[i] = sum(z[natural_neighbors] .* weights)
    end
    return zi
end

function simplex_interpolate(x::Vector{T}, y::Vector{T}, z::Vector{T}, xi::Vector{T}, yi::Vector{T}) where {T <: Number}
    # Perform Delaunay triangulation of input data points
    points = collect([x'; y']')
    tri = delaunay(points)
    return simplex_interpolate(tri, z, xi, yi)
end

function simplex_interpolate(x::Vector{T}, y::Vector{T}, z::Vector{T}, xi::T, yi::T) where {T <: Number}
    return simplex_interpolate(x, y, z, [xi], [yi])
end

struct SimplexInterpolator{T}
    x::Vector{T}
    y::Vector{T}
    z::Vector{T}
    tri::Triangulation
    mode::String
end

function SimplexInterpolator(xs, ys, zs; mode=MODE_SIMPLEX)
    global MODE_NATURAL, MODE_SIMPLEX
    @assert mode in [MODE_NATURAL, MODE_SIMPLEX] "mode needs to be one of [$MODE_NATURAL, $MODE_SIMPLEX], but is $mode."
    points = collect([xs'; ys']')
    if eltype(points) == Float32
        points = Float64.(points)
    end
    if eltype(zs)==Float32
        zs = Float64.(zs)
    end
     if eltype(xs)==Float32
        xs = Float64.(xs)
    end
    if eltype(ys)==Float32
        ys = Float64.(ys)
    end
    tri = delaunay(points)
    return SimplexInterpolator(xs, ys, zs, tri, mode)
end

function interpolate(itp, xi::AbstractVector, yi::AbstractVector)
    return simplex_interpolate(itp.tri, itp.z, xi, yi; mode=itp.mode)
end

function interpolate(itp, xi::T, yi::T) where {T <: Number}
    return interpolate(itp, [xi], [yi])
end
# Standard libraries
using Random, LinearAlgebra, SparseArrays
# Geospatial libraries
import GeoInterface as GI, GeometryOps as GO, LibGEOS as LG
import Proj, ArchGDAL
import GeoFormatTypes: EPSG
# File reading
using Shapefile, DataFrames, GeoDataFrames
using Rasters # and raster processing
# Spatially binning the geometries
using SortTileRecursiveTree, AbstractTrees


# You can cache a Raster using an LRU cache like this,
# but in this case it doesn't help performance much because the chunks are too small.
# cached_dem = Rasters.rebuild(dem; data = DiskArrays.cache(parent(dem); maxsize = 3000))

"""
    minimum_elevation(masked_dem)

Find the minimum elevation in a masked DEM, and return a 2-tuple of the point at which it occurs `(long, lat)`.
"""
function minimum_elevation(masked_dem)
    # First, check if the masked DEM is empty
    if isempty(masked_dem)
        ext = GI.extent(masked_dem)
        return (ext.X[2] - ext.X[1], ext.Y[2] - ext.Y[1])
    end
    # Otherwise, construct an iterator that skips missing values
    value_iterator = skipmissing(masked_dem)
    if isnothing(iterate(value_iterator)) # check that the iterator is not empty
        ext = GI.extent(masked_dem)
        return (ext.X[2] - ext.X[1], ext.Y[2] - ext.Y[1])
    end
    # finally, find the minimum value and its index
    # this works only because Rasters.jl has a custom skipmissing that allows for this
    minval, minidx = findmin(value_iterator)
    # return the point at which the minimum value occurs
    return DimPoints(masked_dem)[minidx]
end

"""
    is_extent_suitable(raster, node; max_n_chunks = 100^2, chunk_size = 128^2)

Check if the extent of a node is suitable for loading into memory.
"""
function is_extent_suitable(raster, node; max_n_chunks=100^2, chunk_size=128^2)
    # check if the extent of this particular node is less than 100x100 chunks.
    # To do so, we convert the extent to indices within the raster,
    # and then check if the product of the lengths of the indices is less than `max_n_chunks * chunk_size`.
    # TODO: in future, maybe simply `crop` -> `DiskArrays.eachchunk(parent(cropped))` -> count chunks?
    # That seems like a cleaner approach...
    indices = Rasters.DD.dims2indices(raster, GI.extent(node))
    return prod(length.(indices)) < max_n_chunks * chunk_size
end
# NOTE: 
# I tried tuning the max number of chunks here,
# and 100x100 seemed like a sweet spot.
# 200x200 almost doubled the runtime for some reason,
# so I presume there's a balance between memory usage / download speed / decoding speed
# and runtime.
# Plus, since glaciers are nicely clustered in space,
# we can afford to do this relatively easily.

"""
    extract_idxs(node)

Extract the indices of the geometries that an `STRNode` covers.

Returns a vector of integers that are indices of the geometries that made the tree.
"""
function extract_idxs(node)::Vector{Int}
    return mapreduce(Base.Fix2(getproperty, :indices), vcat, AbstractTrees.Leaves(node))
end


function get_minimum_elevations(dem, geometries; progress=true)
    prog = progress ? ProgressMeter.Progress(length(geometries); desc="Getting minimum elevations") : nothing

    # Create a (default) STRtree with the geometries
    # This usually means a leaf node size of around 10 geometries,
    # which is not too onerous.
    tree = SortTileRecursiveTree.STRtree(geometries)

    # Preallocate the minimum elevations array so we can assign to it
    minimum_elevations = fill((NaN, NaN), length(geometries))
    # Now, we construct an iterator over the tree (using AbstractTrees.jl),
    # performing a depth first search over the tree from parent to child.
    # We also use the predicate to stop evaluating once we find a node that is suitable in extent.

    # Note that the predicate causes iteration to assume a leaf node if it's **`false`**,
    # otherwise it will continue iterating over children.
    # So, the predicate is that if the extent is unsuitable, we keep iterating (true),
    # but if it's suitable, we stop iterating over children (false).
    iterator = AbstractTrees.PreOrderDFS(
        Base.Fix1(!is_extent_suitable, dem), # predicate to stop the iterator
        tree.rootnode                        # the root node of the tree
    )

    # Loop over the tree in pre-order depth-first search
    for node::Union{SortTileRecursiveTree.STRNode,SortTileRecursiveTree.STRLeafNode} in iterator
        # We perform rasterization if the node's extent is sufficiently small,
        # or if the node is a leaf node.
        loaded_dem = if is_extent_suitable(dem, node)
            @debug("Loading DEM for node with extent $(GI.extent(node))")
            # Here, we know the node is smaller than the maximum size we're willing to load,
            # so we load the section of the DEM covering the node into memory,
            # and use that for the zonal operation.
            Rasters.read(Rasters.crop(dem; to=GI.extent(node)::GI.Extents.Extent)::Raster)::Raster
        elseif node isa SortTileRecursiveTree.STRLeafNode
            @debug("Using full DEM for node with extent $(GI.extent(node))")
            # Here, we know the node is a leaf node, and that the extent of the node is larger than
            # the maximum size we're willing to load, 
            # so we just use the full (lazy) DEM for the zonal operation.
            dem
        else
            # The node is neither suitable nor a leaf node,
            # so we skip it.
            @debug("Skipping node with extent $(GI.extent(node))")
            continue
        end
        @debug("Zonalling")
        # Extract the indices of the geometries that this node covers
        indices = extract_idxs(node)
        # Perform the zonal operation
        try
        minimum_elevations[indices] .= Rasters.zonal(
            minimum_elevation, loaded_dem;
            of=view(geometries, indices),
            skipmissing=false,  # provide a full masked Raster to the zonal operation
            threaded=true,      # thread over the geometries
            progress=false,     # don't show progress since we are already showing global progress
        )::Vector{Union{Missing,Tuple{Float64,Float64}}}
        catch e
            println("Error zonalling node with extent $(GI.extent(node))")
            display(Rasters.read(Rasters.crop(dem; to=GI.extent(node))))
            rethrow(e)
        end



        # Update the progress meter
        progress && ProgressMeter.update!(prog, prog.counter + length(indices))
    end

    return minimum_elevations
end


struct UnitCartesianFromGeographic <: CoordinateTransformations.Transformation
end

# median time 19 ns
"""
    (::UnitCartesianFromGeographic)(geographic_point)

Convert a geographic point (longitude, latitude) to a unit Cartesian 3D point on the sphere.

# Arguments
- `geographic_point`: A point with longitude (x) and latitude (y) in degrees.

# Returns
- `GeometryBasics.Point3`: The corresponding 3D point on the unit sphere.

# Notes
- Longitude is mapped to the azimuthal angle θ (in degrees).
- The polar angle ϕ is computed as 90° minus the latitude.
- Uses `sincosd` for efficient and numerically stable computation.
"""
function (::UnitCartesianFromGeographic)(geographic_point)
    # Longitude is directly translatable to a spherical coordinate
    # θ (azimuth)
    θdeg = GI.x(geographic_point)
    # The polar angle is 90 degrees minus the latitude
    # ϕ (polar angle)
    ϕdeg = 90 - GI.y(geographic_point)
    # Here, we use `sincosd` for two reasons:
    # 1. It's faster than calling `sind` and `cosd` separately
    # 2. `sind` and `cosd` are more numerically stable than 
    # `sin(deg2rad(...))` and `cos(deg2rad(...))` are.
    # If we have coordinates available in radians, using them directly is better.
    sinϕ, cosϕ = sincosd(ϕdeg)
    sinθ, cosθ = sincosd(θdeg)
    # Since this is the unit sphere, the radius is assumed to be 1,
    # and we don't need to multiply by it.
    return GeometryBasics.Point3(
        sinϕ * cosθ,
        sinϕ * sinθ,
        cosϕ
    )
end

struct GeographicFromUnitCartesian <: CoordinateTransformations.Transformation
end

# median time 36 ns
function (::GeographicFromUnitCartesian)(xyz::AbstractVector)
    @assert length(xyz) == 3 "GeographicFromUnitCartesian expects a 3D Cartesian vector"
    sph = SphericalFromCartesian()(xyz)
    return GeometryBasics.Point2(
        rad2deg(sph.θ),
        rad2deg(sph.ϕ)
    )
end

"""
    UnitSphericalCap{T}

A spherical cap, defined by a point on the sphere and a radius.
"""
struct UnitSphericalCap{T}
    point::GeometryBasics.Point3{T}
    radius::T
end

function UnitSphericalCap(point::GeometryBasics.Point3{T1}, radius::T2) where {T1,T2}
    return UnitSphericalCap{promote_type(T1, T2)}(point, radius)
end

UnitSphericalCap(point, radius) = UnitSphericalCap(GI.trait(point), point, radius)

UnitSphericalCap(::GI.PointTrait, point, radius) = UnitSphericalCap(GOC.Spherical(), GI.PointTrait(), point, radius)

function UnitSphericalCap(m::GOC.Spherical{T}, ::GI.PointTrait, point, radius) where {T}
    radius_on_unit_sphere = radius / m.radius
    return UnitSphericalCap(UnitCartesianFromGeographic()(point), radius_on_unit_sphere)
end


function spherical_cap_edge_points(cap::UnitSphericalCap{T}, num_points::Int=100, R::T=one(T)) where {T}
    # The problem with perpendicular-vector-based approaches is that they are vulnerable to oversensitivity, and the normalization actually makes the problem worse.
    # Instead, we construct points on the edge of the cap by rotating the points on the edge of the North Pole,
    # then rotating them to the correct orientation.

    edge_points = Vector{GeometryBasics.Point3{T}}(undef, num_points)
    # NOTE: this is the physical convention of spherical coordinates,
    # where θ is the polar angle (colatitude) and φ is the azimuthal angle.
    center_theta = acos(cap.point[3])
    center_phi = atan(cap.point[2], cap.point[1])
    # see https://math.stackexchange.com/questions/275134/rotation-matrix-sending-the-north-pole-to-a-specified-point-on-the-sphere 
    # for the rotation matrix
    # we directly convert to a static matrix to pre-compute and materialize the matrix,
    rotation = Rotations.RotZX{T}(center_phi + pi / 2, center_theta) |> Rotations.RotMatrix{3}

    sinphi, cosphi = sincos(cap.radius)
    for (i, theta) in enumerate(Base.range(zero(T), stop=2π, length=num_points))
        sintheta, costheta = sincos(theta)
        # Generate a series of Cartesian points on the edge of the North Pole
        cartesian_point = GeometryBasics.Point3{T}(sinphi * costheta, sinphi * sintheta, cosphi)
        # Rotate the points to the correct orientation, such that the "new" north pole is actually the cap's point
        edge_points[i] = LinearAlgebra.normalize(rotation * cartesian_point) * R
    end

    return edge_points
end


"""
    to_latlong_polygon(cap::UnitSphericalCap, n_points::Int)

Convert a `UnitSphericalCap` to a geographic (latitude/longitude) polygon.

# Arguments
- `cap::UnitSphericalCap`: The spherical cap to convert.
- `n_points::Int`: Number of points to use for approximating the polygon edge.

# Returns
- `GI.Polygon`: A polygon in latitude/longitude (EPSG:4326) representing the cap.

The function samples `n_points` along the edge of the cap, converts them to geographic coordinates, and returns a polygon.
"""
function to_latlong_polygon(cap::UnitSphericalCap, n_points::Int)
    points = GeographicFromUnitCartesian().(spherical_cap_edge_points(cap, n_points))
    return GI.Polygon([GI.LinearRing(points)]; crs=GFT.EPSG(4326))
end

#=
# Disabled: importing MakieCore here makes the package unprecompilable. Makie and MakieCore both
# define `deepcopy(::Observables.Observable)`, and the resulting method overwriting is an error
# during module precompilation. Plot a cap by converting it with `to_latlong_polygon` instead.
#
import MakieCore
function MakieCore.convert_arguments(::Type{MakieCore.Mesh}, cap::UnitSphericalCap)
    # get the edge points
    points = spherical_cap_edge_points(cap, 50, 1.0)
    # add the center point
    push!(points, cap.point)
    # create the faces
    faces = [GeometryBasics.TriangleFace(i, i + 1, length(points)) for i in 1:(length(points)-2)]
    return (GeometryBasics.normal_mesh(points, faces),)
end
=#

function trace_downstream(start_id::Integer, ids, next_ids; maxiters=length(ids))
    # prog = ProgressUnknown(; dt = 0.1, desc = "Tracing downstream...")
    if start_id == 0
        return Base.nonmissingtype(eltype(ids))[]
    end
    visited_ids = Base.nonmissingtype(eltype(ids))[]
    current_id = start_id
    current_idx = searchsortedfirst(ids, start_id)
    iter_count = 0
    while current_id != 0 && iter_count <= maxiters
        current_idx = searchsortedfirst(ids, current_id)
        push!(visited_ids, current_id)
        current_id = next_ids[current_idx]
        # next!(prog)
        iter_count += 1
    end
    # finish!(prog)
    return visited_ids
end

function haversine_distance(x, y)
    return haversine_distance(GOC.Spherical(), x, y)
end

function haversine_distance(manifold::GOC.Spherical, x, y)
    return haversine_distance(manifold, GI.trait(x), x, GI.trait(y), y)
end

haversine_distance(manifold::GOC.Spherical, ::Union{GI.AbstractTrait, Nothing}, x, ::Union{GI.AbstractTrait, Nothing}, y) = error("Not implemented yet for $(GI.trait(x)) and $(GI.trait(y)), please file an issue if you need this.")

haversine_distance(manifold::GOC.Spherical, ::GI.PointTrait, x, ::GI.PointTrait, y) = Distances.Haversine(manifold.radius)((GI.x(x), GI.y(x)), (GI.x(y), GI.y(y)))

haversine_distance(manifold::GOC.Spherical, ::GI.PointTrait, x, ::GI.AbstractCurveTrait, y) = GO.applyreduce(min, GI.PointTrait(), y; init = Inf) do yp
    haversine_distance(manifold, x, yp)
end



"""
    flux_accumulate!(river_inputs, id, nextdown_id, headbasin, majorbasin_id)

Accumulate river fluxes by propagating flows downstream through a river network.

# Arguments
- `river_inputs`: DimensionalArray of basin/river inputs [m s⁻¹] with dimensions [Ti, ID]
- `id`: Vector of river/basin IDs 
- `nextdown_id`: Vector of next downstream river/basin IDs (0 indicates no downstream reach)
- `headbasin`: Boolean vector indicating headwater basins (true = no upstream reaches)
- `majorbasin_id`: Vector of major/main basin IDs used for parallel processing

# Details
For each major basin, this function:
1. Identifies headwater basins and their downstream connections
2. Iteratively propagates flows downstream by adding upstream flows to downstream reaches
3. Continues until all basins in the network have been processed
4. Performs accumulation in-place, modifying the input `river_inputs` array

The function processes each major basin in parallel using multiple threads.
flux_accumulate! is insanely fast.
"""
function flux_accumulate!(river_inputs, id, nextdown_id, headbasin, majorbasin_id; accumulated_length = false)
    # sum river inputs going downstream to get total river flux     

    Threads.@threads for majorbasin in unique(majorbasin_id) #[2 min]
        # majorbasin = unique(majorbasin_id)[1]
        index = majorbasin_id .== majorbasin

        id0 =  id[index]
        nextdown_id0 = nextdown_id[index]
        headbasin0 = headbasin[index]
        
        # identify head basins that do not flow into other basins
        touched = nextdown_id0 .== 0

        # identify head basins that flow into other basins
        flowdown_now = headbasin0 .& (.!touched) # these basins are ready to flow down
        flowdown_later = copy(flowdown_now) # these basins are waiting for other basins catch fill up
        
        # exclude basins that will later be touched by flowdown
        flowdown_now_ids = setdiff(unique(nextdown_id0[flowdown_now]), unique(nextdown_id0[.!(touched .| flowdown_now)]))
        flowdown_now[flowdown_now] = in.(nextdown_id0[flowdown_now], Ref(flowdown_now_ids)) 
        flowdown_later = flowdown_later .& (.!flowdown_now)

        while any(flowdown_now)
            id1 = id0[flowdown_now]
            nextdown_id1 = nextdown_id0[flowdown_now]

            # a loop here seems to be faster than a view
            if accumulated_length
                 for to_id in unique(nextdown_id1)
                    from_id = id1[nextdown_id1.==to_id]
                    river_inputs[At(to_id)] += maximum(river_inputs[At(from_id)])
                end
            else
                for (from_id, to_id) in zip(id1, nextdown_id1) 
                    river_inputs[:, At(to_id)] .+= parent(river_inputs[:, At(from_id)])
                end
            end  

            touched[:] = flowdown_now .| touched

            if all(touched)
                break
            else
                flowdown_now[:] = (in.(id0, Ref(unique(nextdown_id0[flowdown_now]))) .| flowdown_later) .& (.!touched)
                flowdown_later[:] = flowdown_now

                flowdown_now_ids = setdiff(unique(nextdown_id0[flowdown_now]), unique(nextdown_id0[.!(touched .| flowdown_now)]))
                flowdown_now[flowdown_now] = in.(nextdown_id0[flowdown_now], Ref(flowdown_now_ids))
                flowdown_later[:] = flowdown_later .& (.!flowdown_now)
            end
        end

        # this is a sanity check that all basins were touched
        if !all(touched)
            @warn "for major basin $(majorbasin) $(sum(touched)) basins were not touched"
        end
    end
    return river_inputs
end

function linear_reservoir_impulse_response_monthly(Tb)
    #Tb = 45 #days
    Vb = 1 #mm

    Qb =[]
    for i in 1:10000
        if i < 30
            Qsb = 1 #mm/day
        else
            Qsb = 0 #mm/day
        end
        Qb = push!(Qb, Vb/Tb)
        Vb += (Qsb - Qb[i])
    end


    c = 1
    dd = 30
    Qfrac = []
    for _ in 1:10
        push!(Qfrac, round(sum(Qb[c:c+dd-1]) ./ sum(Qb), digits=2))
        c += dd
    end

    # Trim the zero tail. For a long enough residence time none of the ten monthly bins rounds to
    # zero, and `findfirst` then returns `nothing` -- the old `findfirst(Qfrac .== 0) - 1` threw
    # `MethodError: no method matching -(::Nothing, ::Int64)` for Tb >= ~60 days. Keep everything
    # when there is no zero bin to trim at.
    first_zero = findfirst(iszero, Qfrac)
    if !isnothing(first_zero)
        Qfrac = Qfrac[1:first_zero-1]
    end

    return Qfrac
end


function apply_vector_impulse_resonse!(M, impulse_resonse)
    nir = length(impulse_resonse)
    m, n, k = size(M)
    Q1 = zeros(eltype(M), (k + nir - 1))

    for i in 1:m
        for j in 1:n
            Q1 .= 0
            Q0 = M[i, j, :]
            for r in 1:k
                Q1[r:r+nir-1] .+= Q0[r] .* impulse_resonse
            end
            M[i, j, :] = Q1[1:k]
        end
    end
end



# load in river reaches
function river_reaches(rivers_paths; col_names=nothing)
    rivers = fill(DataFrame(),length(rivers_paths))

    Threads.@threads for i in eachindex(rivers_paths)
        if isnothing(col_names)
            rivers[i] = GeoDataFrames.read(rivers_paths[i])
        else
            rivers[i] = GeoDataFrames.read(rivers_paths[i])[:, col_names]
        end
    end
    
    rivers = reduce(vcat, rivers)
    sort!(rivers, [:COMID])
    return rivers
end


function compute_population(gpw_ras, rivers, country_polygons, gmax, runoff, buffer_radii; progress=true)
    m2deg = 0.0000089932
    
    countries = unique(country_polygons[!, :country])
    filter!(!=(""), countries)

    dgmax = Dim{:gmax}(gmax; metadata=Dict(:longname => "maximum monthly glacier fraction of total river flux", :units => "fraction"))
    drunoff = Dim{:runoff}(runoff; metadata=Dict(:longname => "minimum runoff flux during gmax", :units => "m³/s"))

    dbuffer = Dim{:buffer}(buffer_radii; metadata=Dict(:longname => "buffer radius", :units => "m"))
    dcountry = Dim{:country}(countries; metadata=Dict(:longname => "country", :units => ""))

    population = zeros(dcountry, dgmax, drunoff,dbuffer; name="population")

    source_crs1 = GFT.EPSG(4326)
    target_crs1 = GFT.EPSG(3857)

    country_polygons[!, :geometry] = _validate_convert.(country_polygons[:, :geometry])

    prog = progress ? ProgressMeter.Progress(length(population); desc="Extracting population for each county...") : nothing
    
    for country in dcountry # having issues with PROJError: pipeline: Pipeline: too deep recursion when using Threads.@threads

        country_polygon = country_polygons[country_polygons.country.==country, :geometry]
        country_polygon = reduce(LibGEOS.union, country_polygon)

        for gmax in dgmax
            for runoff in drunoff

                river_idx = (rivers[:, :gmax_avg] .>= gmax) .& (rivers[:, :country] .== country) .& (rivers[:, :runoff_max_avg] .>= runoff)

                if !any(river_idx)
                    ProgressMeter.update!(prog, prog.core.counter + length(dbuffer))
                    continue
                end

                geom1 = rivers[river_idx, :geometry]

                # extract points
                geom1 = reduce(vcat, _get_point_subset.(geom1; subsample=10))
    
                # project to web-mercador for equal area projection
                geom1 = GO.reproject(geom1, source_crs1, target_crs1; always_xy=true)

                for buffer_radius in dbuffer
                    # printstyled("\ncountry: $country, gmax: $gmax, buffer_radius: $buffer_radius \n"; color=:red)
                    
                    # buffer
                    geom = GO.buffer(geom1, buffer_radius)

                    # project back to geographic coordinates
                    geom = GO.reproject(geom, target_crs1, source_crs1; always_xy=true)

                    # GO.union is hopeless
                    #GO.union(geom[1], geom[2]; target=GI.PolygonTrait())
                    #GI.Polygon(GI.LinearRing(collect(GI.getpoint(geom[1]))))

                    # convert geometry to LibGEOS to use thier union function
                    geom = _validate_convert.(geom)
                    geom = filter(!isnothing, geom)

                    # union all the polygons
                    try
                        geom = reduce(LibGEOS.union, geom)
                    catch
                        printstyled("\nerror unioning geometry for country: $country, gmax: $gmax, buffer_radius: $buffer_radius\n"; color=:red)
                        progress ? next!(prog) : nothing
                        continue
                    end

                    # trim to country boundaries
                    geom = LibGEOS.intersection(country_polygon, geom)

                    total_pop = Rasters.zonal(sum, gpw_ras; of=geom, skipmissing=true)
                    population[At(country), At(gmax), At(runoff), At(buffer_radius)] += total_pop

                    progress ? next!(prog) : nothing
                end
            end
        end
    end
    finish!(prog)

    return population
end


 function _get_point_subset(geom; subsample=10)

    pts = GI.getpoint(geom)
    if length(pts) < subsample
        pts = GO.centroid(geom)
    else
        pts = collect(pts)
        idx = collect(1:10:length(pts))
        idx[end] = length(pts)
        pts = pts[idx]
    end

    return pts
end 

function _validate_convert(geom)
    try
        geom = GI.convert(LibGEOS, geom)
        if !LibGEOS.isValid(geom)
            geom = nothing
        end
    catch
        geom = nothing
    end
    return geom
end


# ============================================================================
# Glacier Routing Module
#
# This module handles the routing of glacier meltwater through river networks.
#
# Key Functionality:
# - Imports necessary packages for glacier routing analysis
# - Sets up paths and constants for data processing
# - Finds lowest elevation points for each glacier to determine meltwater entry points
# - Maps glaciers to drainage basins and river networks
# - Traces downstream flow paths for glacier meltwater
# - Calculates and accumulates glacier runoff through river networks
# - Exports results to NetCDF format for further analysis
#
# Data Sources:
# - Glacier outlines and attributes
# - Digital elevation models (30m resolution)
# - River network data from MERIT Hydro
# - GLDAS flux data (corrected and uncorrected versions)
#
# Workflow:
# 1. Identify lowest points in glacier boundaries
# 2. Map glaciers to drainage basins and nearest rivers
# 3. Trace downstream flow paths
# 4. Calculate glacier runoff inputs to river segments
# 5. Accumulate runoff through river networks
# 6. Export results to NetCDF format
# ============================================================================

# Sets up paths and constants needed for glacier routing analysis.
#
# Configuration includes:
# - File paths for input/output data
# - Physical constants for unit conversions
# - River network data paths (corrected and uncorrected versions)
# - GLDAS flux data configuration
#
# This section establishes all necessary file paths and constants required for the
# glacier routing workflow, including paths to glacier outlines, river networks,
# and output locations for processed data.
"""
    glacier_routing_file(surface_mask) -> String

Path of the per-glacier routing table (`RiverID`, `RiverIDTrace`, ...) that [`glacier_routing`](@ref)
writes for the glaciers of `surface_mask`, next to that mask's individual-glacier outlines.
"""
glacier_routing_file(surface_mask) =
    replace(pathlocal[Symbol("$(surface_mask)_individual")], ".gpkg" => "_routing.arrow")

"""
    glacier_rivers_file(surface_mask) -> String

Path of the river reaches fed by the glaciers of `surface_mask`, written by [`glacier_routing`](@ref).
The land-surface fluxes routed onto those reaches (`land_surface_model_routing.jl`) are saved beside it
with the same stem, so every glacier-river product carries the mask in its name.
"""
glacier_rivers_file(surface_mask) =
    joinpath(pathlocal[:river], "riv_pfaf_MERIT_Hydro_v07_Basins_v01_$(surface_mask).arrow")

function glacier_routing(; surface_mask="glacier", route_rgi2ocean=["19"])
    Unitful.register(MyUnits)
    glacier_summary_riverflux_file = replace(pathlocal[:glacier_summary], ".nc" => "_riverflux.nc")
   
    # Load local configuration paths
    paths = pathlocal

    # Derived data paths
    glacier_lowest_point_path = replace(paths[Symbol("$(surface_mask)_individual")], ".gpkg" => "_lowest_point.gpkg")
    glacier_routing_path = glacier_routing_file(surface_mask)

    # River network paths
    # Note: Two options for river data - using corrected or uncorrected data

    # Uncorrected river paths (currently used)
    rivers_paths = allfiles(paths[:river]; fn_endswith="MERIT_Hydro_v07_Basins_v01.shp", fn_startswith="riv_pfaf")
    glacier_rivers_path = glacier_rivers_file(surface_mask)

    # GLDAS flux data configuration
    use_corrected_flux = false
    if use_corrected_flux 
        # Paths for corrected flux data using in-situ observations
        river_Q_path = joinpath(paths.data_dir, "rivers/Collins2024/Qout_pfaf_ii_GLDAS_COR_M_1980-01_2009-12_utc/")
        river_Q_files = readdir(river_Q_path; join=true)
        lev2_basins = [splitpath(fn)[end][11:12] for fn in river_Q_files]

        # Backfill missing data with ensemble simulations
        river_Q_path_ens = joinpath(paths.data_dir, "rivers/Collins2024/Qout_pfaf_ii_GLDAS_ENS_M_1980-01_2009-12_utc/")
        river_Q_files_ens = readdir(river_Q_path_ens; join=true)

        for fn in river_Q_files_ens
            if !any(splitpath(fn)[end][11:12] .== lev2_basins)
                river_Q_files = push!(river_Q_files, fn)
            end
        end
    else
        # Use uncorrected ensemble simulations only
        river_Q_path_ens = joinpath(paths.data_dir, "rivers/Collins2024/Qout_pfaf_ii_GLDAS_ENS_M_1980-01_2009-12_utc/")
        river_Q_files = readdir(river_Q_path_ens; join=true)
    end


# Find the lowest elevation point for each glacier
#
# Identifies the lowest elevation point within each glacier boundary by overlaying
# glacier geometries on a 30m Copernicus DEM. These points represent where glacier
# meltwater likely enters the river network.
#
# The function:
# 1. Reads glacier outline geometries from GeoPackage
# 2. Overlays geometries on the 30m DEM to find minimum elevation points
# 3. Extracts longitude and latitude coordinates of these points
# 4. Saves the augmented glacier dataset with minimum elevation coordinates
#
# Only runs if the output file doesn't already exist to avoid redundant processing.
if !isfile(glacier_lowest_point_path)
    glaciers = GeoDataFrames.read(paths[Symbol("$(surface_mask)_individual")])
    dem = Raster("/mnt/devon-r2/shared_data/copernicus-dem-30m/DEM.vrt"; lazy=true)
    mes = get_minimum_elevations(dem, glaciers.geometry)
    glaciers.min_elev_long = first.(mes)
    glaciers.min_elev_lat = last.(mes)
    GeoDataFrames.write(glacier_lowest_point_path, glaciers)
end


# Find drainage basins and rivers for each glacier
#
# Maps glaciers to their drainage basins and river networks by:
# 1. Identifying which drainage basin contains each glacier's lowest point
# 2. Finding the closest river to each glacier's lowest point
# 3. Tracing the downstream flow path for glacier meltwater
# 4. Mapping glaciers to all downstream river segments
# 5. Identifying which rivers terminate at the ocean
# 6. Assigning continent and country information to each river segment
#
# The function uses spatial indexing for efficient querying and processes data in parallel
# where possible. Results are saved to disk to avoid redundant processing in future runs.
#
# Only executes if the output file doesn't already exist.
if !isfile(glacier_rivers_path)
#begin
    glaciers = GeoDataFrames.read(glacier_lowest_point_path)

    # read the file containing the basins
    basins = GeoDataFrames.read(paths[:river_basins])

    # construct points from minimum latitude and longitude
    pts = GI.Point.(glaciers.min_elev_long, glaciers.min_elev_lat)

    # Convert the `basin` geometries to Julia geometries,
    # so that access is super fast
    basin_geometries = GO.tuples(basins.geometry)

    # Build a spatial index for the basins
    tree = SortTileRecursiveTree.STRtree(basins.geometry; nodecapacity=3)

    # Iterate over the points and find the first basin that intersects with each point
    # Prefilter by the basins that could possibly intersect with each point
    basin_indices = zeros(Union{Int64,Missing}, size(pts))
    @showprogress dt=1 desc="Locating drainage basins..." Threads.@threads :greedy for (idx, pt) in enumerate(pts)
        # Query the tree for the point
        potential_basin_idxs = SortTileRecursiveTree.query(tree, pt)
        # Find the first basin that intersects with the point
        first_intersecting_basin = findfirst(Base.Fix1(GO.intersects, pt), view(basin_geometries, potential_basin_idxs))
        # Store the basin index, or missing if no intersection was found
        basin_indices[idx] = isnothing(first_intersecting_basin) ? missing : potential_basin_idxs[first_intersecting_basin]
    end

    valid = .!ismissing.(basin_indices)
    glaciers[!, :BasinID] = Vector{Union{Missing,Int64}}(fill(missing, length(pts)))
    glaciers[valid, :BasinID] = basin_indices[valid]

    # identify the river that each glacier is closest to
    # load in river reaches
    col_names = [ "geometry",  "COMID", "NextDownID", "up1", "lengthkm"]
    riv = GeoDataFrames.read.(rivers_paths)
    rivers = reduce(vcat, r[:,col_names] for r in riv)
    tree = SortTileRecursiveTree.STRtree(rivers.geometry; nodecapacity=3)

    river_indices = zeros(Union{Int64,Missing}, size(pts))
    @showprogress dt=1 desc="Locating closest river..." Threads.@threads :greedy for (idx, pt) in enumerate(pts)
        # Query the tree for the point
        potential_river_idxs = Int64[];
        buff = 0;
        while isempty(potential_river_idxs)
            buff += 1_000 # distance in meters
            potential_river_idxs = SortTileRecursiveTree.query(tree, to_latlong_polygon(UnitSphericalCap(pt, buff), 16))
        end

        if length(potential_river_idxs) == 1
            river_indices[idx] = potential_river_idxs[1]
        else
            river_distances = [haversine_distance(pt, rivers.geometry[river_idx]) for river_idx in potential_river_idxs]
            river_indices[idx] = first(potential_river_idxs[river_distances .== minimum(river_distances)])
        end
    end

    glaciers[!, :RiverID] = rivers.COMID[river_indices]

    # set the river id to 0 for glaciers that flow directly to the ocean
    for reg in route_rgi2ocean
        idx = glaciers.O1Region .== reg
        if any(idx)
            glaciers[idx, :RiverID] .= 0
        end
    end

    glaciers[!, :RiverIDTrace] .= [Int64[]]
    @showprogress dt=1 desc="Tracing downstream..." for r in eachrow(glaciers)
    #r = eachrow(glaciers)[1]
        r.RiverIDTrace = trace_downstream(r.RiverID, rivers.COMID, rivers.NextDownID; maxiters=length(rivers.COMID))
    end

    # save the glaciers with the routing information
    GeoDataFrames.write(glacier_routing_path, glaciers)

    # identify just those rivers that receive glacier meltwater
    glacier_melt_rivers = unique(reduce(vcat, glaciers.RiverIDTrace))
    rivers = rivers[in.(rivers.COMID, Ref(glacier_melt_rivers)), :]

    # map glaciers to each downstream rivers
    rivers[!,  :glacier_table_row] .= [Vector{Int64}()]

    @showprogress dt=1 desc="Mapping glaciers to downstream rivers..." Threads.@threads :greedy for r in eachrow(rivers)
        r.glacier_table_row =  findall(in.(Ref(r.COMID), glaciers.RiverIDTrace))
    end

    rivers[!, :RGIId] = [glaciers.RGIId[idx] for idx in rivers.glacier_table_row]

    # determine which rivers flow to the ocean
    river2ocean_path = joinpath(paths.data_dir, "rivers/Collins2024/riv_pfaf_ii_MERIT_Hydro_v07_Basins_v01_coast")
    river2ocean_files = readdir(river2ocean_path; join=true)
    river2ocean_files = filter(fn -> Base.contains(fn, ".shp"), river2ocean_files)
    rivers[!, :ocean_terminating] .= false

    river2ocean = GeoDataFrames.read(river2ocean_files[1])
    for fn in river2ocean_files[2:end]
        try
            river2ocean = reduce(vcat, (river2ocean , GeoDataFrames.read(fn)))
        catch
        end
    end

    idx = [any(rid .== river2ocean.COMID) for rid in rivers.COMID]
    rivers[idx, :ocean_terminating] .= true

    rivers = GeoDataFrames.read(glacier_rivers_path)
    rivers[!, :continent] .= ""
    rivers[!, :country] .= ""
    rivers[!, :centroid] = GO.centroid.(rivers.geometry)
    tree = SortTileRecursiveTree.STRtree(rivers[!, :centroid])

    # Assign continent to each river segment
    continents = GeoDataFrames.read(paths.continents)

    @time Threads.@threads for r in eachrow(continents) # [4.5 min]
        query_result = SortTileRecursiveTree.query(tree, GI.extent(r.geometry))
        index = GO.intersects.(rivers.centroid[query_result], Ref(r.geometry))
        rivers[query_result[index], :continent] .= r.CONTINENT
    end

    # Assign country to each river segment
    countries = GeoDataFrames.read(paths[:countries])
    @time Threads.@threads for r in eachrow(countries) # [2.5 min]
        query_result = SortTileRecursiveTree.query(tree, GI.extent(r.geometry))
        index = GO.intersects.(rivers.centroid[query_result], Ref(r.geometry))
        rivers[query_result[index], :country] .= r.SOVEREIGNT
    end

    GeoDataFrames.write(glacier_rivers_path, rivers[:, Not(:centroid)])
end


# Route glacier runoff through river networks.
#
# This function processes glacier runoff data and routes it through downstream river networks:
# 1. Loads glacier runoff data from a NetCDF file
# 2. Converts runoff from Gt to kg/s
# 3. Interpolates values to monthly time steps
# 4. Maps glacier runoff to river segments
# 5. Accumulates runoff downstream through the river network
# 6. Saves the resulting river flux data to a NetCDF file
#
# Arguments
# - `glacier_routing_path`: Path to glacier routing information
# - `pathlocal[:glacier_summary]`: Path to input glacier runoff data
# - `glacier_summary_riverflux_file`: Path for output river flux data
# - `glacier_rivers_path`: Path to river network data
# - `paths`: Dictionary of project paths

begin # [4 minutes]
    vars2route=["runoff"]

    glaciers = copy(GeoDataFrames.read(glacier_routing_path)) # getting a strange issue with scope later on so using copy

    glacier_flux_nc = NCDataset(pathlocal[:glacier_summary])
    summary_mask = get(glacier_flux_nc.attrib, "surface_mask", nothing)
    summary_mask == String(surface_mask) ||
        error("$(pathlocal[:glacier_summary]) was built for surface mask $(repr(summary_mask)), not " *
              "$(repr(String(surface_mask))): rerun glacier_summary_file with surface_mask=$(repr(String(surface_mask)))")
        
    glacier_flux0 = Dict()

    for tsvar in vars2route
    #tsvar = "runoff"
        var0 = nc2dd(glacier_flux_nc[tsvar])

        # convert from Gt to kg³/s
        dTi = dims(var0, :Ti)
        Δdate = Second.(unique(Day.(val(dTi)[2:end] .- val(dTi)[1:end-1])))

        # Gt per month
        glacier_flux0[tsvar] = var0[:, 2:end]
        glacier_flux0[tsvar][:,:] =  diff(parent(var0), dims=2)
        glacier_flux0[tsvar] = uconvert.(u"kg", glacier_flux0[tsvar]) ./ Δdate
    end

    # interpolate to the 1st of every month
    glacier_flux = Dict()
    for tsvar in vars2route # [30 seconds]
        var0 = glacier_flux0[tsvar]
        dTi = dims(var0, :Ti)
        drgi = dims(var0, :rgiid)
        
        yrange = extrema(Dates.year.(dTi))
        t0 = [Dates.DateTime(y,m,1) for m in 1:12, y in yrange[1]:yrange[2]][:];
        t0 = t0[(t0 .>= dTi[1]) .& (t0 .<= dTi[end])]
        dTi0 = Ti(t0)
        glacier_flux[tsvar] = zeros(drgi,dTi0)* u"kg/s"

        tms = Dates.datetime2epochms.(val(dTi));
        dt = unique(tms[2:end] .- tms[1:end-1]);
        @assert length(dt) == 1

        t0ms = Dates.datetime2epochms.(t0)

        #scale for efficiency of interpolation
        t0ms = (t0ms .- tms[1]) ./  dt[1] .+ 1

        Threads.@threads for rgi in drgi
            A = var0[At(rgi), :]
            valid = .!isnan.(A)
            if any(valid)
                    itp = Interpolations.extrapolate(Interpolations.interpolate(A[valid], Interpolations.BSpline(Interpolations.Quadratic())), NaN * u"kg/s")
                foo = itp.(t0ms .- findfirst(valid) .+ 1)
                (tsvar == "runoff") ? (foo[foo .< 0* u"kg/s"] .= 0 * u"kg/s") : nothing
                glacier_flux[tsvar][At(rgi), :] = foo
            else
                glacier_flux[tsvar][At(rgi), :] = fill(NaN, length(t0ms))* u"kg/s"
            end
        end
    end

    # this takes 1 min to load
    rivers = GeoDataFrames.read(glacier_rivers_path)
    rivers[!, :basin02] = floor.(Int16, rivers.COMID./1000000)
    rivers[!, :headbasin] = rivers.up1 .== 0
    rivers[!, :longitude] = getindex.(GO.centroid.(rivers.geometry), 1)
    rivers[!, :latitude] = getindex.(GO.centroid.(rivers.geometry), 2)

    # if there is no input from upstream then it is a head basin
    # This is needed as only a subset of the full river network is routed
    rivers[!, :headbasin] .|= .!in.(rivers.COMID, Ref(unique(rivers.NextDownID)))

    dcomid = Dim{:COMID}(rivers.COMID)
    dTi = dims(glacier_flux[first(vars2route)], :Ti)
    river_inputs = Dict()
    river_flux = Dict()

    for tsvar in vars2route
        tsvar = "runoff"
        river_inputs[tsvar] = zeros(eltype(glacier_flux[tsvar]), (dTi, dcomid))

        # need to sum as multiple glaciers can flow into the same river
        for r in eachrow(glaciers[glaciers.RiverID .!= 0, :])
            # r = eachrow(glaciers[glaciers.RiverID .!= 0, :])[1]
            river_inputs[tsvar][:, At(r.RiverID)] .+= glacier_flux[tsvar][At(r.RGIId), :]
        end

        # accululate inputs downstream
        river_flux[tsvar] = flux_accumulate!(river_inputs[tsvar], rivers.COMID, rivers.NextDownID, rivers.headbasin, rivers.basin02) #[< 1 second]

        # save data as netcdf
        # data needs to be Float32 (not Real) for saving as netcdf
        river_flux[tsvar] = Float32.(river_flux[tsvar])
    end

    dstack = DimStack(getindex.(Ref(river_flux), vars2route); name=Symbol.(vars2route))
    
    # save to netcdf [having permission issues]
    if isfile(glacier_summary_riverflux_file)
        rm(glacier_summary_riverflux_file)
    end

    NCDataset(glacier_summary_riverflux_file, "c") do ds

        data_dims = dims(dstack)

        # First all dimensions need to be defined
        for dim in data_dims
            dname = string(DimensionalData.name(dim));
            defDim(ds, dname, length(dim))
        end

        # now add the variables
        for dim in data_dims
            dname = string(DimensionalData.name(dim));
            d = defVar(ds, dname, val(val(dim)), (dname,))
            if DateTime <: eltype(dim)
                d.attrib["cf_role"] = "timeseries_id";
            end
        end

        for vaname in keys(dstack)
            v = defVar(ds, "$vaname", ustrip.(parent(dstack[vaname])), string.(DimensionalData.name.(data_dims)))
            v.attrib["units"] = string(Unitful.unit(dstack[vaname][1]))
        end

        # add latitude and longitude [This is a hack until I can get geometry into a NetCDF]
        defVar(ds, "latitude", rivers.latitude, string.(DimensionalData.name.(data_dims[2:2])))
        defVar(ds, "longitude", rivers.longitude, string.(DimensionalData.name.(data_dims[2:2])))

        # add global attributes
        ds.attrib["title"] = "river flux from glacier model"
        ds.attrib["version"] = "beta - " * Dates.format(now(), "yyyy-mm-dd")
        ds.attrib["featureType"] = "timeSeries";

        println("glacier output saved to $glacier_summary_riverflux_file")
    end;
end
end
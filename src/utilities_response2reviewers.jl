# Import required packages
using CSV
using DataFrames
import GlacierHydrology
using DimensionalData
using Dates
using CairoMakie
using Statistics
using ProgressMeter
using NCDatasets
using FileIO

"""
    load_wgms_data(paths)

Load and preprocess WGMS mass balance and glacier data.

# Arguments
- `paths`: Named tuple or struct with keys :wgms_mb, :wgms_glacier, :path2glacier, :project_dir

# Returns
- Tuple (wgms_mb, wgms_glacier, ds, glaciers) where wgms_mb and wgms_glacier are DataFrames,
  ds is a DimStack from the project NetCDF, and glaciers is a GeoDataFrame.

# Examples
```julia
julia> paths = (wgms_mb="...", wgms_glacier="...", path2glacier="...", project_dir="...")
julia> wgms_mb, wgms_glacier, ds, glaciers = load_wgms_data(paths)
```
"""
function load_wgms_data(paths)
    # Load WGMS mass balance and glacier attribute tables
    wgms_mb = CSV.read(paths.wgms_mb, DataFrame)
    wgms_glacier = CSV.read(paths.wgms_glacier, DataFrame)

    # Augment WGMS data table with latitude, longitude, and geometry (point)
    ia, ib = GlacierHydrology.intersectindices(collect(wgms_mb.glacier_id), collect(wgms_glacier.id))
    wgms_mb[!, :latitude] = wgms_glacier[ib, :latitude]
    wgms_mb[!, :longitude] = wgms_glacier[ib, :longitude]
    wgms_mb[!, :geometry] = GlacierHydrology.GI.Point.(wgms_mb.longitude, wgms_mb.latitude)

    # Load model data (netcdf) and individual glacier polygons
    path2outfile = joinpath(GlacierHydrology.pathlocal[:project_dir], "Gardner2025_glacier_2deg.nc")

    println("glacier model data last modified: $(Dates.unix2datetime(mtime(path2outfile))): $path2outfile")
    ds = GlacierHydrology.netcdf2dimstack(path2outfile)
    glaciers = GlacierHydrology.GeoDataFrames.read(paths.path2glacier)
    
    return wgms_mb, wgms_glacier, ds, glaciers
end

"""
    preallocate_model_columns!(wgms_mb)

Pre-allocate columns for modeled quantities in the WGMS DataFrame (in-place).

# Arguments
- `wgms_mb`: DataFrame of WGMS mass balance data; new columns are added for runoff, dm, smb, and seasonal variants.

# Returns
- The modified DataFrame (in-place).

# Examples
```julia
julia> preallocate_model_columns!(wgms_mb)
julia> # wgms_mb now has :runoff, :dm, :smb, etc.
```
"""
function preallocate_model_columns!(wgms_mb)
    wgms_mb[!, :runoff] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
    wgms_mb[!, :dm] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
    wgms_mb[!, :smb] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
    wgms_mb[!, :runoff_winter] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
    wgms_mb[!, :dm_winter] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
    wgms_mb[!, :smb_winter] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
    wgms_mb[!, :runoff_summer] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
    wgms_mb[!, :dm_summer] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
    wgms_mb[!, :smb_summer] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
    wgms_mb[!, :geotile_centroid] = Vector{Union{Missing,GlacierHydrology.GI.Point}}(missing, nrow(wgms_mb))
    wgms_mb[!, :geotile_glacier_area] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
    wgms_mb[!, :geotile_area] = Vector{Union{Missing,Float64}}(missing, nrow(wgms_mb))
end

"""
    setup_geotile_grid(ds)

Setup gridded longitude and latitude intervals and cell areas for geotile calculations.

# Arguments
- `ds`: DimStack or similar with X and Y dimensions defining the grid spacing

# Returns
- Raster of cell areas (e.g. from Rasters.cellarea) in the same CRS as the grid.

# Examples
```julia
julia> cell_area = setup_geotile_grid(ds)
julia> # Use for area-weighted aggregation
```
"""
function setup_geotile_grid(ds)
    geotile_width_out = abs(dims(ds, :X)[1] - dims(ds, :X)[2])
    lond = X(GlacierHydrology.DimensionalData.Dimensions.Lookups.Sampled(
        -180+geotile_width_out/2:geotile_width_out:180-geotile_width_out/2,
        sampling=GlacierHydrology.DimensionalData.Dimensions.Lookups.Intervals(GlacierHydrology.DimensionalData.Dimensions.Lookups.Center());
        metadata=Dict("long_name" => "longitude", "units" => "degrees")))
    latd = Y(GlacierHydrology.DimensionalData.Dimensions.Lookups.Sampled(
        -90+geotile_width_out/2:geotile_width_out:90-geotile_width_out/2,
        sampling=GlacierHydrology.DimensionalData.Dimensions.Lookups.Intervals(GlacierHydrology.DimensionalData.Dimensions.Lookups.Center());
        metadata=Dict("long_name" => "latitude", "units" => "degrees")))
    cell_area = GlacierHydrology.Rasters.cellarea(GlacierHydrology.Rasters.Raster(zeros(lond, latd); crs=GlacierHydrology.GeoFormatTypes.EPSG(4326)))
    return cell_area
end

"""
    extract_time_series(ds, longitude, latitude, begin_date, end_date)

Extract time series for runoff, dm, and smb over the specified time interval.
Returns a tuple of (runoff_cumulative, dm_cumulative, smb_cumulative).
"""
function extract_time_series(ds, longitude, latitude, begin_date, end_date)
    runoff = ds.runoff[X=Near(longitude), Y=Near(latitude), Ti=Touches(begin_date, end_date)]::DimVector{Float64}
    runoff_cumulative = cumsum(runoff)
    dm = ds.dm[X=Near(longitude), Y=Near(latitude), Ti=Touches(begin_date, end_date)]::DimVector{Float64}
    dm_cumulative = cumsum(dm)
    smb = ds.smb[X=Near(longitude), Y=Near(latitude), Ti=Touches(begin_date, end_date)]::DimVector{Float64}
    smb_cumulative = cumsum(smb)
    return runoff_cumulative, dm_cumulative, smb_cumulative
end

"""
    calculate_metrics(runoff_cumulative, dm_cumulative, smb_cumulative, geotile_glacier_area)

Calculate runoff, dm, and smb metrics from cumulative time series.
Returns a tuple of (runoff, dm, smb).
"""
function calculate_metrics(runoff_cumulative, dm_cumulative, smb_cumulative, geotile_glacier_area)
    runoff = (runoff_cumulative[end] - runoff_cumulative[1]) / geotile_glacier_area / 1000
    dm = (dm_cumulative[end] - dm_cumulative[1]) / geotile_glacier_area / 1000
    smb = (smb_cumulative[end] - smb_cumulative[1]) / geotile_glacier_area / 1000
    return runoff, dm, smb
end

"""
    process_observation_row!(r, ds, cell_area)

Process a single observation row to extract synthetic mass balance data.
"""
function process_observation_row!(r, ds, cell_area)
    # Only process if observation has valid date interval
    if !ismissing(r.begin_date) && !ismissing(r.end_date)
        # Retrieve geotile glacier area and total cell area
        r.geotile_glacier_area = ds.area[X=Near(r.longitude), Y=Near(r.latitude)]::Float64
        r.geotile_area = cell_area[X=Near(r.longitude), Y=Near(r.latitude)]::Float64

        # Extract time series for the full observation interval
        runoff_cumulative, dm_cumulative, smb_cumulative = extract_time_series(ds, r.longitude, r.latitude, r.begin_date, r.end_date)

        # Store synthetic outputs if non-empty
        if length(runoff_cumulative) > 0
            r.runoff, r.dm, r.smb = calculate_metrics(runoff_cumulative, dm_cumulative, smb_cumulative, r.geotile_glacier_area)
            r.geotile_centroid = GlacierHydrology.GI.Point(refdims(runoff_cumulative)[1][1], refdims(runoff_cumulative)[2][1])

            # If a midseason date exists, also compute winter/summer metrics
            if !ismissing(r.midseason_date)
                r.midseason_date = DateTime(r.midseason_date)

                # Compute winter values (start to midseason)
                runoff_cum_winter, dm_cum_winter, smb_cum_winter = extract_time_series(ds, r.longitude, r.latitude, r.begin_date, r.midseason_date)
                if length(runoff_cum_winter) > 0
                    r.runoff_winter, r.dm_winter, r.smb_winter = calculate_metrics(runoff_cum_winter, dm_cum_winter, smb_cum_winter, r.geotile_glacier_area)
                end

                # Compute summer values (midseason to end)
                runoff_cum_summer, dm_cum_summer, smb_cum_summer = extract_time_series(ds, r.longitude, r.latitude, r.midseason_date, r.end_date)
                if length(runoff_cum_summer) > 0
                    r.runoff_summer, r.dm_summer, r.smb_summer = calculate_metrics(runoff_cum_summer, dm_cum_summer, smb_cum_summer, r.geotile_glacier_area)
                end
            end
        end
    end
end

"""
    add_validity_flags!(wgms_mb)

Add validity flags for comparing synthetic and observed values.
"""
function add_validity_flags!(wgms_mb)
    # Calculate differences between synthetic and observed values and add validity flags
    wgms_mb[!, :synth_minus_obs] = wgms_mb[:, :smb] .- wgms_mb[:, :annual_balance]
    wgms_mb[!, :synth_minus_obs_winter] = wgms_mb[:, :smb_winter] .- wgms_mb[:, :winter_balance]
    wgms_mb[!, :synth_minus_obs_summer] = wgms_mb[:, :smb_summer] .- wgms_mb[:, :summer_balance]
    wgms_mb[!, :both_valid] = (.!ismissing.(wgms_mb.annual_balance) .& .!ismissing.(wgms_mb.dm)) .& .!isnan.(wgms_mb.synth_minus_obs)
    wgms_mb[!, :both_valid_winter] = (.!ismissing.(wgms_mb.winter_balance) .& .!ismissing.(wgms_mb.dm_winter)) .& .!isnan.(wgms_mb.synth_minus_obs_winter)
    wgms_mb[!, :both_valid_summer] = (.!ismissing.(wgms_mb.summer_balance) .& .!ismissing.(wgms_mb.dm_summer)) .& .!isnan.(wgms_mb.synth_minus_obs_summer)
end

"""
    add_runoff_and_dm!(wgms_mb, ds)

Augment WGMS table with modeled synthetic mass balance, runoff, and other metrics.
This is the main function that orchestrates the data processing.
"""
function add_runoff_and_dm!(wgms_mb, ds)
    # Pre-allocate columns for modeled quantities
    preallocate_model_columns!(wgms_mb)
    
    # Setup gridded longitude and latitude intervals and cell areas
    cell_area = setup_geotile_grid(ds)

    # Main loop over each observation row in WGMS data
    for r = eachrow(wgms_mb)
        process_observation_row!(r, ds, cell_area)
    end

    # Add validity flags
    add_validity_flags!(wgms_mb)
end

"""
    calculate_wgms_summary_statistics(wgms_mb, index)

Calculate and print summary statistics for WGMS mass balance data.
Returns a tuple of (index_all, index_area, df2, wgms_annual_coverage).
"""
function calculate_wgms_summary_statistics(wgms_mb, index)
    # Count all valid annual balance observations, print summary
    index_all = .!ismissing.(wgms_mb.annual_balance)
    println("total number of annual_balance observations in the WGMS dataset: $(sum(index_all)) [$(length(unique(wgms_mb.glacier_name[index_all]))) named glaciers]")

    # Count all valid study-period observations, print summary
    println("total number of observations with dates in the WGMS dataset within the study period: $(sum(index)) [$(length(unique(wgms_mb.glacier_name[index]))) named glaciers]")

    # Calculate coverage as fraction of global glacier area observed each year
    index_area = index .& .!ismissing.(wgms_mb.area)
    gdf = groupby(wgms_mb[index_area, :], :year)
    df2 = DataFrames.combine(gdf, :area => sum => :area_sum, :annual_balance => std => :annual_balance_std, :synth_minus_obs => std => :synth_minus_obs_std)
    
    return index_all, index_area, df2
end

"""
    calculate_annual_coverage(df2, ds)

Calculate annual coverage statistics and return coverage percentages.
"""
function calculate_annual_coverage(df2, ds)
    wgms_annual_coverage = df2.area_sum ./ sum(ds.area[:]) * 100
    println("average annual coverage: $(round(mean(wgms_annual_coverage), digits=2))%")
    return wgms_annual_coverage
end

"""
    plot_area_coverage(df2, wgms_annual_coverage)

Create and save a plot of area coverage over time.
"""
function plot_area_coverage(df2, wgms_annual_coverage)
    f = GlacierHydrology._publication_figure(columns=1, rows=1)
    ax1 = CairoMakie.Axis(f[1, 1]; ylabel="global glacier area observed [%]")
    lines!(df2.year, wgms_annual_coverage)
    fill_between!(ax1, df2.year, wgms_annual_coverage, 0)
    xlims!(ax1, minimum(df2.year), maximum(df2.year))
    ylims!(ax1, (0, nothing))
    
    fname = joinpath(GlacierHydrology.pathlocal.figures, "wgms_mb_annual_coverage.png")
    CairoMakie.save(fname, f)
    return f
end

"""
    calculate_glacier_count_statistics(wgms_mb, index_area)

Calculate statistics for number of unique glaciers observed each year.
Returns df2 with glacier count data.
"""
function calculate_glacier_count_statistics(wgms_mb, index_area)
    gdf = groupby(wgms_mb[index_area, :], :year)
    df2 = DataFrames.combine(gdf, :glacier_name => (x -> length(unique(x))) => :unique_glaciers, :annual_balance => std => :annual_balance_std, :synth_minus_obs => std => :synth_minus_obs_std)
    println("average number of unique glaciers per year: $(round(Int, mean(df2.unique_glaciers)))")
    return df2
end

"""
    plot_glacier_count_timeseries(df2)

Create and save a plot of glacier count over time.
"""
function plot_glacier_count_timeseries(df2; rgi_count=274531)
    # Use RGI v6 and v7 glacier counts for global context
    #rgi6_count = CONFIG.rgi6_count
    #rgi7_count = CONFIG.rgi7_count

    f = GlacierHydrology._publication_figure(columns=1, rows=1)
    ax1 = CairoMakie.Axis(f[1, 1]; ylabel="glacier count")
    ax2 = CairoMakie.Axis(f[1, 1], yaxisposition=:right, ylabel="% of global glacier count")

    lines!(ax1, df2.year, df2.unique_glaciers)
    fill_between!(ax1, df2.year, df2.unique_glaciers, 0)
    # Optionally: percentage axis for glacier counts
    # lines!(ax2, df2.year, (df2.unique_glaciers ./ rgi_count) * 100)
    # fill_between!(ax2, df2.year, (df2.unique_glaciers ./ rgi7_count) * 100, 0)

    ymax = ceil(maximum(df2.unique_glaciers) / 10) * 10
    xlims!(ax1, minimum(df2.year), maximum(df2.year))
    ylims!(ax1, (0, ymax))
    xlims!(ax2, minimum(df2.year), maximum(df2.year))
    ylims!(ax2, (0, ymax / rgi_count * 100))
    
    fname = joinpath(GlacierHydrology.pathlocal.figures, "wgms_mb_glacier_count.png")
    CairoMakie.save(fname, f)
    return f
end

"""
    calculate_spatial_variogram(wgms_mb, mincount=3)

Calculate spatial variogram of annual balance as a function of distance.
Returns variogram data and distance bins.
"""
function calculate_spatial_variogram(wgms_mb; mincount=3)
    wgms_mb[!, :geometry] = GlacierHydrology.GI.Point.(wgms_mb.longitude, wgms_mb.latitude) # Ensure all rows have point geometry

    # Define distance bins (in meters) for variogram calculation
    delta_distance = vcat(0, 10_000 .* cumsum(collect(1:2:16)) .+ 10_000)
    ddistance = Dim{:distance}(delta_distance)
    index_variogram = .!ismissing.(wgms_mb.annual_balance) .& (wgms_mb.year .> 1999)
    variogram = zeros(ddistance)

    # Calculate mean absolute annual-balance difference as a function of distance
    for buffer_radius in ddistance[2:end]
        intermediate_year = Vector{Float64}()
        gdf = groupby(wgms_mb[index_variogram, :], :year)

        buffer_ind = findfirst(delta_distance .== buffer_radius)
        if buffer_ind == 1
            buffer_inner = 0
        else
            buffer_inner = delta_distance[buffer_ind - 1]
        end

        # Group by year to calculate statistics per year
        @showprogress desc = "Calculating variogram of annual balance for buffer radius $buffer_radius m..." Threads.@threads for df2 in gdf
            
            intermediate_pt = Vector{Float64}()

            for r in eachrow(df2)
                # Donut-shaped spatial window: within buffer_radius and outside buffer_inner
                cap1 = GlacierHydrology.UnitSphericalCap(r.geometry, buffer_radius)
                cap2 = GlacierHydrology.UnitSphericalCap(r.geometry, buffer_inner)
                polygon1 = GlacierHydrology.to_latlong_polygon(cap1, 30)
                polygon2 = GlacierHydrology.to_latlong_polygon(cap2, 30)

                index_within_radius = GlacierHydrology.GO.intersects.(Ref(polygon1), df2.geometry) .& .!GlacierHydrology.GO.intersects.(Ref(polygon2), df2.geometry)
                if sum(index_within_radius) > mincount
                    absdiff = abs.(df2.annual_balance[index_within_radius] .- r.annual_balance)
                    mean_absdiff = mean(absdiff[absdiff .!= 0])

                    try
                        push!(intermediate_pt, mean_absdiff)
                    catch
                        display(mean_absdiff)
                        error("Error calculating mean absolute difference for year $(df2.year[1]) and buffer radius $buffer_radius m")
                    end
                end
            end
            if length(intermediate_pt) > 0
                push!(intermediate_year, mean(intermediate_pt))
            end
        end
        variogram[distance=At(buffer_radius)] = mean(intermediate_year)
    end
    
    return variogram, ddistance
end

"""
    calculate_geotile_misfit_statistics(wgms_mb, index)

Calculate geotile misfit statistics and return summary data.
"""
function calculate_geotile_misfit_statistics(wgms_mb, index)
    gdf = groupby(wgms_mb[index, :], :geotile_centroid)
    df2 = DataFrames.combine(gdf,
        :synth_minus_obs => std => :synth_minus_obs_std,
        :synth_minus_obs => length => :nobs,
        :geotile_area => first => :geotile_area,
        :synth_minus_obs => mean => :synth_minus_obs_mean)
    
    return df2
end

"""
    calculate_misfit_summary_statistics(wgms_mb, index, df2)

Calculate and print summary statistics for misfit analysis.
Returns radius and MAD values.
"""
function calculate_misfit_summary_statistics(wgms_mb, index, df2)
    # Compute mean radius of geotiles with at least 65 observations
    min_observations = 1
    index_count = df2.nobs .>= min_observations
    radius = round(Int, sqrt(mean(df2.geotile_area[index_count]) / pi) / 1000)
    println("average radius of geotiles with at least $min_observations observations: $(radius) km")
    
    synth_minus_obs_mad = round(mean(abs.(wgms_mb[index, :synth_minus_obs])), digits=2)
    println("mean absolute difference between in situ and this study: $(synth_minus_obs_mad)")
    
    return radius, synth_minus_obs_mad
end

"""
    plot_variogram(ddistance, variogram, radius, synth_minus_obs_mad)

Create and save a plot of the spatial variogram.
"""
function plot_variogram(ddistance, variogram, radius, synth_minus_obs_mad; plot_point=true)
    f = GlacierHydrology._publication_figure(columns=1, rows=1)
    ax1 = CairoMakie.Axis(f[1, 1]; ylabel="mean absolute difference [m w.e. yr⁻¹]", xlabel="distance to other observations [km]")
    scatterlines!(ax1, ddistance.val ./ 1000, variogram.data)
    plot_point ? plot!(ax1, Point(Float64(radius), synth_minus_obs_mad), color=:red) : nothing
    
    fname = joinpath(GlacierHydrology.pathlocal.figures, "wgms_mb_variogram.png")
    CairoMakie.save(fname, f)
    return f
end

"""
    create_2d_histogram(wgms_mb, index, bound_lower=-5, bound_upper=3)

Create a 2D histogram (heatmap) of observed vs modeled mass balance.
Returns the figure and density data.
"""
function create_2d_histogram(wgms_mb; bound_lower = -5, bound_upper=3, colormap=Reverse(:viridis))
    bins = bound_lower:0.25:bound_upper
    df = GlacierHydrology.binstats(wgms_mb, [:annual_balance, :dm], [bins, bins], :both_valid; missing_bins=true, col_function=sum)
    x = unique(df[:, 1])
    y = unique(df[:, 2])
    density = reshape(df.both_valid_sum, length(x), length(y))

    f = GlacierHydrology._publication_figure(columns=1, rows=1)
    ax1 = CairoMakie.Axis(f[1, 1]; ylabel="in situ [m w.e. yr⁻¹]", xlabel="this study [m w.e. yr⁻¹]")
    hm = heatmap!(x, y, density; colormap)
    Colorbar(f[1, 2], hm)
    lines!(ax1, [bound_lower, bound_upper], [bound_lower, bound_upper], color=:black) # 1:1 reference line
    display(f)
    
    return f, density
end

"""
    create_misfit_histograms(wgms_mb, index)

Create histograms for misfit statistics by season/period.
"""
function create_misfit_histograms(wgms_mb, index)
    hist(wgms_mb[index, :synth_minus_obs])
    hist(wgms_mb[wgms_mb[:, :both_valid_summer], :synth_minus_obs_summer])
    hist(wgms_mb[wgms_mb[:, :both_valid_winter], :synth_minus_obs_winter])
end

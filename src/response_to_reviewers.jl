# Load packages and data paths
begin
    using CSV
    using DataFrames
    import GlobalGlacierAnalysis as GGA
    import GlacierVolumeChange as GVC
    using DimensionalData
    using Dates
    using CairoMakie
    using Statistics
    using ProgressMeter
    using NCDatasets
    using FileIO
    using CFTime
    using GeoDataFrames
    import GeometryOps as GO
    
    # Import functions from utilities
    include("utilities_response2reviewers.jl")

    # data paths
    paths = (
        wgms_mb = "/mnt/bylot-r3/data/glacier_mb/WGMS/DOI-WGMS-FoG-2025-02b/data/mass_balance.csv",
        wgms_glacier = "/mnt/bylot-r3/data/glacier_mb/WGMS/DOI-WGMS-FoG-2025-02b/data/glacier.csv",

        path2glacier = GGA.pathlocal[Symbol("glacier_individual")],
    )

    # to include in uncertainty
    paths = merge(GGA.pathlocal, paths)


    # Analysis parameters
    study_start_year = 2000
    min_neighbors_variogram = 3
    geotile_thresholds = [1, 20, 65]
    histogram_bounds = (-5, 3)
    historic_analysis_mincount_range = 5:15

    # RGI glacier counts for context
    rgi6_count = 215547
    rgi7_count = 274531

    error_quantile=0.95
    error_scaling=1.5

    discharge_fractional_error = 0.15

    reference_ensemble_file = GVC.reference_ensemble_file
    ensemble_reference_file = replace(reference_ensemble_file, "_aligned.jld2" => "_synthesized.jld2")

    glacier_summary_file = GGA.pathlocal[:glacier_summary]

    # path2perglacier = replace(ensemble_reference_file, ".jld2" => "_perglacier.jld2")
    path2discharge = GGA.discharge_global_file(GGA.binned_filled_fileparts(reference_ensemble_file).surface_mask)
    path2rgi_regions = paths[:rgi6_regions_shp]

    path2river_flux = joinpath(paths[:river], "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier.arrow")


    path2runs_filled, params = GVC.binned_filled_filepaths(;
        project_id=:v01,
        surface_masks=["glacier", "glacier_rgi7"],
        dem_ids=["best", "cop30_v2"],
        curvature_corrects=[false, true],
        amplitude_corrects=[true],
        binning_methods=["median", "nmad3", "nmad5"],
        fill_params=[1, 2, 3, 4],
        binned_folders=[GVC.analysis_paths(; geotile_width=2).binned, replace(GVC.analysis_paths(; geotile_width=2).binned, "binned" => "binned_unfiltered")],
        include_existing_files_only=true
    )
    path2runs_synthesized = replace.(path2runs_filled, "aligned.jld2" => "synthesized.jld2")
end

# load glacier summary data for context 
glacier_summary_file = GGA.pathlocal[:glacier_summary]

glacier_flux_nc = NCDataset(glacier_summary_file)

glacier_points = Point.(glacier_flux_nc[:longitude][:], glacier_flux_nc[:latitude][:])

basins = vcat(GeoDataFrames.read(GGA.pathlocal[:river_level03_as_basins]), GeoDataFrames.read(GGA.pathlocal[:river_level03_si_basins]))
# find index to basin for each glacier point

index_basin = zeros(Int, length(glacier_points))

Threads.@threads for i in eachindex(basins.geometry)
    ind = GO.contains.(Ref(basins.geometry[i]), glacier_points)
    lock = ReentrantLock()
    if any(ind)
        Base.lock(lock) do
            index_basin[ind] .= i
        end
    end
end

table_basins = Dict()
table_basins["indus"] = [4030033640]
table_basins["ganges + brahmaputra"] = [4030025450, 4030022790]
table_basins["amu darya + syr darya"] = [4030050220, 4030050240]
table_basins["Tarim"] = [4030050210, 4030050980, 4030050270]
table_basins["Lake Balkash + Gobi Interior"] = [4030050230, 4030050290, 4030050430, 3030024310, 4030050320]

time_index = (glacier_flux_nc[:Ti] .> Date("2000-01-01")) .& (glacier_flux_nc[:Ti] .< Date("2025-01-01"))
decyear = collect(GGA.decimalyear.(glacier_flux_nc[:Ti]))

n = sum(time_index);

basin_results = DataFrame(
    basin_name = String[],
    basin_ids = Vector{Int}[],
    net_acc_gtyr = Float64[],
    dm_gtyr = Float64[],
    runoff_gtyr = Float64[]
)

# Dictionary to store raw time series for each basin and variable
basin_timeseries = Dict{String, Dict{Symbol, Vector{Float64}}}()

# Store the time vector for the time series
basin_time_vector = decyear[time_index]

for (basin_name, basin_ids) in table_basins
    println("Basin: $basin_name")

    trend_values = Dict{Symbol, Float64}()
    basin_timeseries[basin_name] = Dict{Symbol, Vector{Float64}}()

    for var in [:net_acc, :dm, :runoff]
        foo = zeros(n)
        for basin_id in basin_ids
            index = index_basin .== findfirst(basins.HYBAS_ID .== basin_id)
            if any(index)
                if :net_acc == var
                    foo += vec(sum(glacier_flux_nc[:acc][index, time_index], dims=1)) .- vec(sum(glacier_flux_nc[:ec][index, time_index], dims=1))
                else
                    foo += vec(sum(glacier_flux_nc[var][index, time_index], dims=1))
                end
            end
        end

        # Store the raw time series
        basin_timeseries[basin_name][var] = foo

        foo_fit = GGA.curve_fit(GVC.offset_trend_seasonal, decyear[time_index] .- mean(decyear[time_index]), foo, GGA.p3)
        trend_values[var] = foo_fit.param[2]

        println("  $var: $(round(foo_fit.param[2], digits=1)) Gt/yr")
    end

    push!(basin_results, (
        basin_name = basin_name,
        basin_ids = basin_ids,
        net_acc_gtyr = trend_values[:net_acc],
        dm_gtyr = trend_values[:dm],
        runoff_gtyr = trend_values[:runoff]
    ))
end





# load hugonnet data
# load in all regions
begin
    hugonnet = DataFrame()
    for rgi = 1:19
        rgi_str = lpad(string(rgi), 2, '0')
        path2hugonnet = "/mnt/bylot-r3/data/glacier_mb/Hugonnet2021/time_series/time_series_$(rgi_str)/dh_$(rgi_str)_rgi60_pergla_rates.csv"
        hugonnet = vcat(hugonnet, CSV.read(path2hugonnet, DataFrame))
    end

    # find annual balance and format similar to wgms_mb for spatial variogram
    hugonnet[!, "start_date"] = DateTime.(getindex.(split.(hugonnet.period, "_"), 1))
    hugonnet[!, "end_date"] = DateTime.(getindex.(split.(hugonnet.period, "_"), 2))

    hugonnet[!, "end_year"] = year.(hugonnet[!, "end_date"])
    hugonnet[!, "start_year"] = year.(hugonnet[!, "start_date"])
    hugonnet[!, "dt_years"] = hugonnet[!, "end_year"] - hugonnet[!, "start_year"]
    hugonnet[!, "year"] = hugonnet[!, "end_year"]
    index = hugonnet[!, "dt_years"] .== 1
    hugonnet = hugonnet[index, :]
    hugonnet[!, :annual_balance] = hugonnet.dhdt * 0.85


    rgi6 = GeoDataFrames.read(GGA.pathlocal.glacier_individual)
    index = indexin(hugonnet.rgiid, rgi6.RGIId)

    # it seems some ids do not match... I think this is because hugonnet updated glaciers that are not in rgi6
    index_not_match = isnothing.(index)
    hugonnet = hugonnet[.!index_not_match, :]

    hugonnet[!, :longitude] = rgi6.CenLon[index[.!index_not_match]]
    hugonnet[!, :latitude] = rgi6.CenLat[index[.!index_not_match]]

subset_random = rand(1:nrow(hugonnet), round(Int, nrow(hugonnet)/100))
variogram, ddistance = calculate_spatial_variogram(deepcopy(hugonnet[subset_random, [:annual_balance, :longitude, :latitude, :year]]); mincount=min_neighbors_variogram);

f_var_hugonnet = plot_variogram(ddistance, variogram, 0, 0)
display(f_var_hugonnet)
save(joinpath(GGA.pathlocal.figures, "wgms_mb_variogram_hugonnet.png"), f_var_hugonnet)

end


# link hugonnet glaciers to wgms_mb
index_into_wgms = indexin(hugonnet.rgiid, wgms_glacier.rgi60_ids)
index_valid = .!isnothing.(index_into_wgms)
hugonnet[!, :wgms_id] .= 0
hugonnet[index_valid, :wgms_id] = wgms_glacier.id[index_into_wgms[index_valid]]


wgms_mb[!, :annual_balance_hugonnet] .= NaN
for r in eachrow(wgms_mb)
    index = findfirst((hugonnet.dt_years .== 1) .& (hugonnet.wgms_id .== r.glacier_id) .& (hugonnet.year .== r.year))
    if !isnothing(index)
        r.annual_balance_hugonnet = hugonnet.dhdt[index] .* 0.85
    end
end

delta = wgms_mb.annual_balance_hugonnet .- wgms_mb.annual_balance
index_valid = (.!ismissing.(delta)) .& (.!isnan.(delta))
f = Figure();
ax = Axis(f[1, 1]; xlabel="wgms minus hugonnet [m w.e. yr⁻¹]", ylabel="count");
hist!(ax, delta[index_valid])
f
save(joinpath(GGA.pathlocal.figures, "wgms_mb_histogram_hugonnet.png"), f)


# Create precipitation raster of maximum monthly precipitation
begin
    path2precip = "/mnt/bylot-r3/data/glacier_mb/GPCP/precip.mon.ltm.1991-2020.nc"

    precip = GGA.Rasters.Raster(path2precip, name=:precip);
    precip_max = dropdims(getindex.(argmax(precip, dims=3), 3), dims=:Ti)

    output_path = joinpath(splitpath(path2precip)[1:end-1]..., "precip.mon_of_max.ltm.1991-2020.tif")
    precip_max = reverse(precip_max, dims=:Y)

    precip_max = rebuild(precip_max, data = UInt8.(precip_max.data))
    precip_max = set(precip_max, X => dims(precip_max, X).val .-180)

    GGA.Rasters.write(output_path, precip_max, force=true)
end;

# =============================================================================
# This script generates the Extended Data Figures for the manuscript.
#
# It runs the discharge and model calibration analyses, using the GlacierHydrology
# package and the glacier volume products of GlacierVolumeChange.jl.
#
# The script is organized by Extended Data Figure number, with each section
# corresponding to a figure or analysis step in the manuscript.
#
# Key steps include:
#   - Setting up analysis parameters and file paths
#   - Calculating global glacier discharge
#   - Calibrating the GEMB model
#
# The figures of the altimetry chain (binning, filling, synthesis, ensemble spread and
# sub-sampling) are scripts in the glacier-volume-paper repository, built on GlacierVolumeChange.jl.
# =============================================================================
begin
    using CSV
    using DataFrames
    import GlacierHydrology
    import GlacierVolumeChange as GVC
    using DimensionalData
    using Dates
    using CairoMakie
    using Statistics
    using ProgressMeter
    using NCDatasets
    using FileIO
    using GeoDataFrames
    
    # Import functions from utilities
    include("utilities_response2reviewers.jl")

    force_remake = false # it is best practice to delete files than to use force_remake
    geotile_width = 2
    force_remake_before = nothing

    single_geotile_test = GlacierHydrology.geotiles_golden_test[15] # nothing  "lat[+78+80]lon[-080-078]" #

    #single_geotile_test = #"lat[+60+62]lon[-148-146]"; #"lat[+62+64]lon[-148-146]"; #"lat[+58+60]lon[-136-134]"; #"lat[+58+60]lon[-138-136]"
    #single_geotile_test = "lat[-44-42]lon[+170+172]" # New Zealand
    downscale_to_glacier_method = "area"
    gemb_run_id = 9
    gemb_physics_member = GlacierHydrology.gemb_reference_physics_member

    ens = GlacierHydrology.ensemble_paths(; geotile_width, physics_member=gemb_physics_member)
    path2runs_synthesized = [ens.synthesized_ref]
    ensemble_reference_file = ens.synthesized_ref
    surface_mask_ref = ens.surface_mask_ref

    # data paths
    paths = (
        wgms_mb="/mnt/bylot-r3/data/glacier_mb/WGMS/DOI-WGMS-FoG-2025-02b/data/mass_balance.csv",
        wgms_glacier="/mnt/bylot-r3/data/glacier_mb/WGMS/DOI-WGMS-FoG-2025-02b/data/glacier.csv",

        path2glacier=GlacierHydrology.pathlocal[Symbol("glacier_individual")],
    )

    # Analysis parameters
    study_start_year = 2000
    min_neighbors_variogram = 3
    geotile_thresholds = [1, 20, 65]
    histogram_bounds = (-5, 3)
    historic_analysis_mincount_range = 5:15

    # RGI glacier counts for context
    rgi6_count = 215547
    rgi7_count = 274531


    discharge_fractional_error = 0.15

    glacier_summary_file = GlacierHydrology.pathlocal[:glacier_summary]

    # path2perglacier = replace(ensemble_reference_file, ".jld2" => "_perglacier.jld2")
    path2discharge = GlacierHydrology.discharge_global_file(surface_mask_ref)

   
    path2rgi_regions = GlacierHydrology.pathlocal[:rgi6_regions_shp]

end;

gemb = GlacierHydrology.gemb_ensemble_dv(; gemb_run_id, surface_mask=surface_mask_ref, physics_member=gemb_physics_member);

# 12. Calculate global glacier solid ice discharge, filled Antarctic data with modeled SMB
discharge = GlacierHydrology.global_discharge_filled(;
    surface_mask=surface_mask_ref,
    discharge_global_fn=GlacierHydrology.discharge_global_file(surface_mask_ref),
    gemb,
    discharge2smb_max_latitude=-60,
    discharge2smb_equilibrium_period=(Date(1979), Date(2000)),
    pscale=1,
    ΔT=0,
    geotile_width=2,
    force_remake_before=GlacierHydrology.gemb_ensemble_built(GlacierHydrology.gemb_info(; gemb_run_id, physics_member=gemb_physics_member),
                                                [surface_mask_ref]),
    force_remake_before_hypsometry=nothing
);

gemb = GlacierHydrology.dv_adjust4discharge!(gemb, discharge)

# [Extended Data Figure 8 & 11d]
#single_geotile_test = GlacierHydrology.geotiles_golden_test[1];

calibrate_to = [:trend, :five_year, :annual, :all];
for calibrate_to in calibrate_to
    (cost, dv_altim, geotiles_in_group) = GlacierHydrology.gemb_calibration(
        path2runs_synthesized,
        gemb[:dv];
        geotile_width,
        geotile_grouping_min_feature_area_km2=100,
        single_geotile_test,
        seasonality_weight=GlacierHydrology.seasonality_weight,
        ΔT_to_pscale_weight=GlacierHydrology.ΔT_to_pscale_weight,
        force_remake_before=GlacierHydrology.gemb_ensemble_built(GlacierHydrology.gemb_info(; gemb_run_id, physics_member=gemb_physics_member),
                                                    [surface_mask_ref]),
        calibrate_to
    );

    print("when calibrate to $(calibrate_to), the model fit is: ")

    f = GlacierHydrology.plot_model_fit(gemb, discharge, dv_altim, cost, geotiles_in_group)
    display(f)

    fname = "$(replace(splitpath(first(path2runs_synthesized))[end], ".jld2" => ""))_cost_$(calibrate_to).png"
    fpath = joinpath(GlacierHydrology.pathlocal[:figures], splitpath(first(path2runs_synthesized))[5], fname)
    save(fpath, f)
end


# [Extended Data Figure 9 & 10]
begin
    outfile_prefix = "Gardner2025_geotiles_rates"
    fname = joinpath(GlacierHydrology.pathlocal[:data_dir], "project_data", "$(outfile_prefix)_myr.gpkg")
    geotiles0 = GeoDataFrames.read(fname)


    f = GlacierHydrology.plot_ref_pscale_ΔT_summary(path2runs_synthesized, ensemble_reference_file, gemb_physics_member;
        rgi2plot=[99], show_title=false, bin_width=0.5,
        distance_from_origin_penalty=GlacierHydrology.gemb_info(; gemb_run_id, physics_member=gemb_physics_member).distance_from_origin_penalty)
    fname = joinpath(GlacierHydrology.pathlocal[:figures], "Extended_Data_Figure_9.png")
    display(f[1])
    CairoMakie.save(fname, f[1])

    # plot a histogram of the rates
    f = GlacierHydrology.plot_hist_gemb_altim_trend_amplitude(geotiles0)
    display(f)

    outfile = joinpath(GlacierHydrology.pathlocal[:figures], "Extended_Data_Figure_10.png")
    CairoMakie.save(outfile, f)    
end

# Compare modeled runoff with the Rounce et al. PyGEM projections
begin
    date_range = DateTime(2000, 1, 1).. DateTime(2025, 1, 1)
    scenario = "ssp126"
    var2extract = "glac_runoff_monthly"
    path2regional = joinpath(GlacierHydrology.pathlocal[:project_dir], "Gardner2025_regional_timseries.nc")
    println("Gardner2025_regional_timseries.nc last modified: $(Dates.unix2datetime(mtime(path2regional)))")
    ds_gardner = GlacierHydrology.netcdf2dimstack(path2regional)
    rain_gardner_include = true
    drgi = dims(ds_gardner[:fac],:rgi);
    dstudy = Dim{:study}(["gardner", "rounce"])
    runoff = zeros(drgi, dstudy)


    for rgi_region in drgi
    #rgi_region = 1
    # Load Rounce data and Gardner output for comparison
        begin
            folder = "/mnt/bylot-r3/data/glacier_mb/PyGEM/CMIP_global_2000_2100/$(lpad(string(rgi_region), 2, '0'))/"

            if !isdir(folder)
                continue
            end

            path2rounce = GlacierHydrology.allfiles(folder; 
                fn_startswith="R$(lpad(string(rgi_region), 2, '0'))_glac_runoff_monthly_$(scenario)",
                fn_endswith="2000_2100_all.nc"
            )

            ds_example = NCDataset(path2rounce[1])
            dTi = Dim{:Ti}(reinterpret.(DateTime, collect(ds_example[:time])))
            dmodel = Dim{:model}(collect(ds_example[:Climate_Model]))

            runoff_rounce =zeros(dTi)

            for path in path2rounce
            #path = path2rounce[1]
                NCDataset(path) do ds_rounce0
                #ds_rounce0 = NCDataset(path)
                    drgi = Dim{:rgi}(collect(ds_rounce0[:RGIId]))
                    runoff_rounce0 = ds_rounce0[var2extract][:, :, :]::Array{Union{Missing,Float64},3}
                    runoff_rounce0 = DimArray(runoff_rounce0, (dTi, drgi, dmodel)) * 1E-9 # m^3/day to km^3
                    runoff_rounce0 = mean(runoff_rounce0, dims=:model)
                    runoff_rounce0 = coalesce.(runoff_rounce0, 0.0)
                    runoff_rounce0 = dropdims(sum(runoff_rounce0, dims=:rgi), dims=(:rgi, :model))
                    runoff_rounce .+= runoff_rounce0
                end
            end

            runoff_rounce = runoff_rounce[Ti = date_range]

            if rain_gardner_include
                runoff_gardner = ds_gardner[:runoff][date=date_range, rgi=At(rgi_region)] .+ ds_gardner[:rain][date=date_range, rgi=At(rgi_region)]
            else
                runoff_gardner = ds_gardner[:runoff][date=date_range, rgi=At(rgi_region)]
            end
            runoff_gardner .-= runoff_gardner[1]

            runoff_gardner = diff(runoff_gardner)

            # Plot comparison
            f = Figure();
            ax = f[1, 1] = Axis(f, ylabel="runoff [Gt]", title="RGI $(rgi_region)");
            lines!(ax, val(dims(runoff_rounce, :Ti)), collect(runoff_rounce); label="Rounce 2022")
            lines!(ax, val(dims(runoff_gardner, :date)), collect(runoff_gardner); label="Gardner 2025")
            f[1, 2] = Legend(f, ax)

            println("mean fractional difference rgi $(rgi_region): Rounce 2022 minus Gardner 2025: $(round(Int, ((sum(runoff_rounce) .- sum(runoff_gardner)) ./ sum(runoff_gardner)) * 100))%")

            runoff[rgi = At(rgi_region), study = At("gardner")] = sum(runoff_gardner)
            runoff[rgi = At(rgi_region), study = At("rounce")] = sum(runoff_rounce)
        end; 
        display(f)
    end

    runoff0 = deepcopy(runoff);
    # convert from total runoff in Gt to Gt/yr
    years = Dates.Day(maximum(date_range) - minimum(date_range)).value / 365.25
    runoff[study=At("gardner")] ./= years
    runoff[study=At("rounce")] ./= years

    # polulate HMA and Global runoff
    runoff[rgi = At(99), study = At("gardner")] = sum(runoff[study=At("gardner")])
    runoff[rgi = At(99), study = At("rounce")] = sum(runoff[study=At("rounce")])

    runoff[rgi = At(98), study = At("gardner")] = sum(runoff[study=At("gardner")][rgi = At(13:15)])
    runoff[rgi = At(98), study = At("rounce")] = sum(runoff[study=At("rounce")][rgi = At(13:15)])


    f = Figure();
    ax = Axis(f[1, 1]; ylabel="PyGEMv1.0.2", xlabel="This study");
    scatter!(ax, vec(runoff[study = At("gardner")]), vec(runoff[study = At("rounce")]));
    xmax = ceil(maximum(runoff[study=At("gardner")]), digits=-2)
    xlims!(ax, 0, xmax)
    ylims!(ax, 0, xmax)
    lines!(ax, [0, xmax], [0, xmax]);
    f

    delta = (runoff[study=At("rounce")] .- runoff[study=At("gardner")]) ./ runoff[study=At("gardner")] * 100

    rgis2plot = [19, 18, 17, 16, 98, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]

    f = GlacierHydrology._publication_figure(columns=1, rows=1);
    ax = Axis(f[1, 1]; yticks = (1:length(rgis2plot), GlacierHydrology.rginum2label.(rgis2plot)), xtickformat = "{:.0f}%", title = "difference in runoff between PyGEMv1.0.2 and this study")

    barplot!(ax, delta[rgi = At(rgis2plot)].data; direction=:x)
    display(f)
    save(joinpath(GlacierHydrology.pathlocal.figures, "Extended_Data_Figure_18.png"), f)

    println("PyGEMv1.0.2 produces $(round(Int,delta[rgi = At(99)]))% runoff compared to this study")
    gardner_minus_rounce = (runoff[study=At("gardner"), rgi = At(99)].-runoff[study=At("rounce"), rgi = At(99)]) ./ runoff[study=At("rounce"), rgi = At(99)] * 100
    println("This study produces $(round(Int,gardner_minus_rounce))% runoff compared to PyGEMv1.0.2")

    rounce_gardner_total_flux = (runoff[study=At("rounce"), rgi=At(99)] .- runoff[study=At("gardner"), rgi=At(99)]) ./ (runoff[study=At("gardner"), rgi=At(99)] + 155) * 100
    println("PyGEMv1.0.2 produces $(round(Int,rounce_gardner_total_flux))% total freshwater flux ompared to this study")
    gardner_rounce_total_flux = (runoff[study=At("gardner"), rgi=At(99)] .- runoff[study=At("rounce"), rgi=At(99)]) ./ (runoff[study=At("rounce"), rgi=At(99)] + 155) * 100
    println("This study produces $(round(Int,gardner_rounce_total_flux))% total freshwater flux ompared to PyGEMv1.0.2")
end

# Compare to WGMS data
begin
    # Load WGMS data
    wgms_mb, wgms_glacier, ds, glaciers = load_wgms_data(paths)

    # Apply function to fill WGMS table with synthetic results
    add_runoff_and_dm!(wgms_mb, ds)

    # Select valid observations within the study period (after 2000)
    index = wgms_mb.both_valid .& (wgms_mb.begin_date .> DateTime(study_start_year, 1, 1))

    # Fill in missing glacier areas by joining to glacier polygons based on intersection
    @showprogress desc = "Adding missing areas..." Threads.@threads for r in eachrow(wgms_mb)[index.&ismissing.(wgms_mb.area)]
        ind = findfirst(GlacierHydrology.GO.intersects.(glaciers.geometry, Ref(r.geometry)))
        if !isnothing(ind)
            r.area = round(glaciers[ind, :Area] * 1e6)
        end
    end

    # Write one validation point per unique glacier to geopackage (point locations only)
    begin
        gdf = groupby(wgms_mb[index, :], :glacier_name)
        df2 = DataFrames.combine(gdf, :geometry => first => :geometry, :geometry => length => :count)
        GlacierHydrology.GeoDataFrames.write(joinpath(GlacierHydrology.pathlocal[:project_dir], "wgms_mb_validation_points.gpkg"),
            df2; crs=GlacierHydrology.EPSG(4326))
    end

    # Attach glacier geometry to validation points and write to geopackage
    begin
        rename!(df2, :geometry => :point_location)
        df2[!, :geometry] .= repeat([glaciers.geometry[1]], nrow(df2))

        @showprogress desc = "Adding glacier geometry to validation points..." Threads.@threads for r in eachrow(df2)
            ind = findfirst(GlacierHydrology.GO.intersects.(glaciers.geometry, Ref(r.point_location)))
            if !isnothing(ind)
                r.geometry = glaciers.geometry[ind]
            end
        end

        # Exclude points with no valid geometry and write result
        no_intersection = df2[:, :geometry] .== Ref(glaciers.geometry[1])
        GlacierHydrology.GeoDataFrames.write(joinpath(GlacierHydrology.pathlocal[:project_dir], "wgms_mb_validation.gpkg"),
            df2[.!no_intersection, Not(:point_location)]; crs=GlacierHydrology.EPSG(4326))
    end

    # --- Peyto Glacier investigation (example, disabled) ---
    #=
    # Identify rows corresponding to Peyto Glacier after 1999
    index_peyto = (wgms_mb.glacier_id .== 57) .& (wgms_mb.year .> 1999) 

    # Compute modeled and observed cumulative mass balance and visualize
    area = ds.area[X=Near(wgms_mb.longitude[findlast(index_peyto)]), Y=Near(wgms_mb.latitude[findlast(index_peyto)])]::Float64
    dm = ds.dm[X=Near(wgms_mb.longitude[findlast(index_peyto)]), Y=Near(wgms_mb.latitude[findlast(index_peyto)]), Ti=Touches(DateTime(2000, 1, 1), DateTime(2025, 1, 1))]::DimVector{Float64}
    dm_cumulative = cumsum(dm) ./ area / 1000
    f = Figure();
    ax1 = Axis(f[1, 1]);
    scatterlines!(wgms_mb[index_peyto, :year], cumsum(wgms_mb[index_peyto, :annual_balance]); label="in situ");
    scatterlines!(GlacierHydrology.decimalyear.(dims(dm_cumulative, :Ti).val), dm_cumulative.data; label="modeled");
    f
    =#

    wgms_mean = mean(wgms_mb[index, :].annual_balance)
    model_mean = mean(wgms_mb[index, :].dm)

    println("FoG glaciers have a mean of $(round(wgms_mean, digits=2)) m w.e. and the extracted model values have a mean of $(round(model_mean, digits=2)) m w.e. ")

    # ---- Summary statistics and coverage for WGMS mass balance data ----

    # Calculate summary statistics
    index_all, index_area, df2 = calculate_wgms_summary_statistics(wgms_mb, index)

    # Calculate and plot area coverage
    wgms_annual_coverage = calculate_annual_coverage(df2, ds)
    f_area = plot_area_coverage(df2, wgms_annual_coverage)
    save(joinpath(GlacierHydrology.pathlocal.figures, "Extended_Data_Figure_12b.png"), f_area)

    # Calculate and plot glacier count statistics
    df2_glacier = calculate_glacier_count_statistics(wgms_mb, index_area)
    f_count = plot_glacier_count_timeseries(df2_glacier)
    save(joinpath(GlacierHydrology.pathlocal.figures, "Extended_Data_Figure_12a.png"), f_count)

    # ---- Calculate spatial variogram of annual balance ----
    variogram, ddistance = calculate_spatial_variogram(wgms_mb; mincount=min_neighbors_variogram)

    # Compute and summarize geotile misfit statistics
    df2 = calculate_geotile_misfit_statistics(wgms_mb, index)
    radius, synth_minus_obs_mad = calculate_misfit_summary_statistics(wgms_mb, index, df2)

    # Plot variogram
    f_var = plot_variogram(ddistance, variogram, radius, synth_minus_obs_mad)
    display(f_var)
    save(joinpath(GlacierHydrology.pathlocal.figures, "Extended_Data_Figure_14.png"), f_var)

    println("at 100km observations have a mean absolute difference of $(round(variogram[distance = At(100_000)], digits=2)) m w.e. yr⁻¹ and the mean absolute difference between in situ and this study is $(round(synth_minus_obs_mad, digits=2)) m w.e. yr⁻¹")

    # Coverage calculation for each year (repeat, probably for reporting)
    index_area = index .& .!ismissing.(wgms_mb.area)
    gdf = groupby(wgms_mb[index_area, :], :year)
    df2 = DataFrames.combine(gdf, :area => sum => :area_sum, :annual_balance => std => :annual_balance_std, :synth_minus_obs => std => :synth_minus_obs_std)
    wgms_annual_coverage = df2.area_sum ./ sum(ds.area[:]) * 100
    println("average annual coverage: $(round(mean(wgms_annual_coverage), digits=2))%")

    # --- 2D histogram (heatmap) of observed vs modeled mass balance ---
    f_hist, density = create_2d_histogram(wgms_mb[index, :]; bound_lower=histogram_bounds[1], bound_upper=histogram_bounds[2], colormap=:tempo)
    save(joinpath(GlacierHydrology.pathlocal.figures, "Extended_Data_Figure_13.png"), f_hist)

    # --- Histograms for misfit statistics by season/period ---
    create_misfit_histograms(wgms_mb, index)
    save(joinpath(GlacierHydrology.pathlocal.figures, "wgms_mb_misfit_histograms.png"), f_hist)
end




# SI Table for comparing runoff from other studies
begin

    # load discharge for each RGI [<1s]
    discharge_rgi = GlacierHydrology.discharge_rgi(path2discharge, path2rgi_regions; fractional_error=discharge_fractional_error);

    # load results for all runs for each RGI [18s]
    runs_rgi = GlacierHydrology.runs2rgi(path2runs_synthesized);

    # fit trends for overlaping periods with RACMO studies

    # Arctic Canada North and South: 2000-2015
    dates4trend = [DateTime(2000, 3, 1), DateTime(2015, 12, 15)]
    runs_rgi_fits = GlacierHydrology.rgi_trends(runs_rgi, discharge_rgi, dates4trend);
    region_fits = GlacierHydrology.region_fit_ref_and_err(runs_rgi_fits, ensemble_reference_file; error_quantile, error_scaling, discharge=discharge_rgi)

    r = region_fits[rgi=At(3), varname=At("runoff"), parameter=At("trend"), error=At(false)];
    println("Arctic Canada North Runoff 2000-2015: $(round(r, digits=2))) Gt/yr")

    r = region_fits[rgi=At(4), varname=At("runoff"), parameter=At("trend"), error=At(false)];
    println("Arctic Canada South Runoff 2000-2015: $(round(r, digits=2))) Gt/yr")

    # Iceland: 2000-2019
    dates4trend = [DateTime(2000, 3, 1), DateTime(2019, 12, 15)]
    runs_rgi_fits = GlacierHydrology.rgi_trends(runs_rgi, discharge_rgi, dates4trend);
    region_fits = GlacierHydrology.region_fit_ref_and_err(runs_rgi_fits, ensemble_reference_file; error_quantile, error_scaling, discharge=discharge_rgi)

    r = region_fits[rgi=At(6), varname=At("runoff"), parameter=At("trend"), error=At(false)];
    println("Iceland Runoff 2000-2015: $(round(r, digits=2))) Gt/yr")

    # Svalbard: 2000-2018
    dates4trend = [DateTime(2000, 3, 1), DateTime(2018, 12, 15)]
    runs_rgi_fits = GlacierHydrology.rgi_trends(runs_rgi, discharge_rgi, dates4trend);
    region_fits = GlacierHydrology.region_fit_ref_and_err(runs_rgi_fits, ensemble_reference_file; error_quantile, error_scaling, discharge=discharge_rgi)

    r = region_fits[rgi=At(7), varname=At("runoff"), parameter=At("trend"), error=At(false)];
    println("Svalbard Runoff 2000-2018: $(round(r, digits=2))) Gt/yr")

    # Southern Andes: 2000-2023
    dates4trend = [DateTime(2000, 3, 1), DateTime(2023, 12, 15)]
    runs_rgi_fits = GlacierHydrology.rgi_trends(runs_rgi, discharge_rgi, dates4trend);
    region_fits = GlacierHydrology.region_fit_ref_and_err(runs_rgi_fits, ensemble_reference_file; error_quantile, error_scaling, discharge=discharge_rgi)

    r = region_fits[rgi=At(17), varname=At("runoff"), parameter=At("trend"), error=At(false)];
    println("Southern Andes Runoff 2000-2023: $(round(r, digits=2))) Gt/yr")
end

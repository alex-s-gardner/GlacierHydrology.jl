# run_all.jl
#
# Glacier mass fluxes, river routing and population analyses, built on the glacier volume-change
# dataset of GlacierVolumeChange.jl. Each major step is numbered below.
#
#  1-10. Glacier volume change: GlacierVolumeChange.jl `create_dataset` (archives, ancillary data,
#        binning, filling, synthesis) and its volume product for every ensemble member. Run it first;
#        this script reads the products.
# 11.  Create glacier model (GEMB) classes for each geotile
# 12.  Calculate global glacier discharge, including filled Antarctic data with modeled SMB
# 13.  Calibrate GEMB model to altimetry data for grouped geotiles
# 14.  Calibrate GEMB model to altimetry data for each synthesized geotile dataset
# 15.  Generate glacier-level summary files with key statistics and metrics
# 16.  Export geotile-level glacier change trends and amplitudes to GIS-compatible files
# 17.  Route glacier runoff through river networks
# 18.  Route land surface model runoff through the glacier-fed river reaches found in 17
# 19.  Calculate gmax (maximum glacier contribution to river flux)
# 20.  Generate point-based figures for gmax visualization
# 21.  Produce regional results for analysis and sharing
# 22.  Generate extended data figures for manuscript
# 23.  Analyze population affected by glacier-fed river changes
# 24.  Create gridded output files for the manuscript
#
# Implements checkpoint logic to avoid redundant processing, allowing efficient restarts.
begin
    import GlacierHydrology
    import GlacierVolumeChange as GVC
    using Unitful
    using CairoMakie
    Unitful.register(GlacierHydrology.MyUnits)
    using Dates
    using DimensionalData

    geotile_width = 2
    single_geotile_test = nothing #GlacierHydrology.geotiles_golden_test[1]
    gemb_run_id = 9
    # The GEMB physics member steps 11-24 run against. Each member's calibration (steps 13-14) is written
    # to its own folder, `gemb_physics/<member>/`, beside the altimetry members; another member is a
    # rerun of steps 11-14 with this changed.
    gemb_physics_member = GlacierHydrology.gemb_reference_physics_member

    # Glacier volume change: one product per ensemble member, written by GVC.create_dataset
    ens = GlacierHydrology.ensemble_paths(; geotile_width, physics_member=gemb_physics_member)
    path2runs_filled = ens.filled
    reference_ensemble_file = ens.filled_ref
    path2runs_synthesized = ens.synthesized
    binned_synthesized_dv_files = ens.gembfit_dv
    binned_synthesized_dv_file_ref = ens.gembfit_dv_ref
    # GEMB is binned, adjusted and calibrated per surface mask; the glacier-level products (steps 15-24)
    # spread the reference member's geotile values over its own inventory's glaciers.
    surface_masks = ens.surface_masks
    surface_mask_ref = ens.surface_mask_ref
end

# 11. Create glacier model (GEMB) classes for each geotile
include("gemb_tiles_binning.jl")

# Steps 12-14 all descend from the GEMB ensemble step 11 writes, so they share one cutoff: anything
# older than the newest of its per-mask files is stale and is rebuilt. Step 11 itself has no cutoff --
# it always rebuilds -- and the cutoff throws if any of those files is missing.
force_remake_before_gemb = GlacierHydrology.gemb_ensemble_built(GlacierHydrology.gemb_info(; gemb_run_id, physics_member=gemb_physics_member),
                                                   String.(surface_masks))

# 12. Calculate global glacier solid ice discharge, filled Antarctic data with modeled SMB
# GEMB, its discharge adjustment and its calibration are per surface mask: each altimetry member is
# fitted against GEMB area-weighted by its own mask's hypsometry.
gemb = Dict{String,Any}()
discharge = Dict{String,Any}()
for m in String.(surface_masks)
    gemb[m] = GlacierHydrology.gemb_ensemble_dv(; gemb_run_id, surface_mask=m, physics_member=gemb_physics_member)
    discharge[m] = GlacierHydrology.global_discharge_filled(;
        surface_mask=m,
        discharge_global_fn=GlacierHydrology.discharge_global_file(m),
        gemb=gemb[m],
        discharge2smb_max_latitude = -60,
        discharge2smb_equilibrium_period=(Date(1979), Date(2000)),
        pscale=1,
        ΔT=0,
        geotile_width=2,
        force_remake_before=force_remake_before_gemb,
        force_remake_before_hypsometry=nothing
    )
    GlacierHydrology.dv_adjust4discharge!(gemb[m], discharge[m])
end

# 13. Calibrate GEMB (Glacier Energy and Mass Balance) model to altimetry data by finding optimal fits for grouped geotiles 
# [4 hrs on 128 cores, 150GB of memory]
# The forcing prior and its penalty weight are the GEMB run's (`GlacierHydrology.GEMB_RUNS`); a run whose
# prior is not derived yet is fitted unpenalized.
GlacierHydrology.gemb_calibration(
    path2runs_synthesized,
    Dict(m => g[:dv] for (m, g) in gemb);
    geotile_width,
    geotile_grouping_min_feature_area_km2=100,
    single_geotile_test,
    seasonality_weight=GlacierHydrology.seasonality_weight,
    ΔT_to_pscale_weight=GlacierHydrology.ΔT_to_pscale_weight,
    force_remake_before=force_remake_before_gemb,
    calibrate_to=:all
)

# 14. Calibrate GEMB model to altimetry data for each synthesized geotile dataset [3.5 min, 150GB of memory]
GlacierHydrology.geotile_synthesis_gembfit_dv(
    path2runs_synthesized, 
    discharge,
    gemb;
    geotile_grouping_min_feature_area_km2 = 100,
    geotile_width,
    force_remake_before=force_remake_before_gemb,
)


vars2plot = ["dm", "dm_altim", "dv"]
rgi2plot = [99]

f = GlacierHydrology.plot_ensemble_rgi_timeseries(vars2plot, path2runs_synthesized;  center_on_dates=DateTime(2002, 5, 1):Month(1):DateTime(2005, 12, 15), rgi2plot)
display.(f);

distance_from_origin_penalty = GlacierHydrology.gemb_info(; gemb_run_id, physics_member=gemb_physics_member).distance_from_origin_penalty
f = GlacierHydrology.plot_ref_pscale_ΔT_summary(path2runs_synthesized, ens.synthesized_ref, gemb_physics_member; rgi2plot, seasonality_weight=GlacierHydrology.seasonality_weight, distance_from_origin_penalty, ΔT_to_pscale_weight=GlacierHydrology.ΔT_to_pscale_weight);
display.(f);

for i in eachindex(rgi2plot)
    fname = "RGI$(rgi2plot[i])_Ws$(GlacierHydrology.seasonality_weight)_Wd$(distance_from_origin_penalty)_Wm2p$(GlacierHydrology.ΔT_to_pscale_weight).png"
    CairoMakie.save(joinpath(GlacierHydrology.pathlocal[:figures], fname), f[i])
end

# 15. Generate glacier level summary nc file with key statistics and metrics for further analysis and sharing  [8 min]
glacier_summary_file = GlacierHydrology.glacier_summary_file(
    binned_synthesized_dv_files,
    binned_synthesized_dv_file_ref;
    error_quantile=0.95,
    geotile_width,
    error_scaling=1.5, 
    reference_period=(DateTime(2000, 4, 1), DateTime(2024, 12, 31)), 
    surface_mask=surface_mask_ref,
    force_remake_before=DateTime(2028, 1, 1),
)
    

# 16. Export geotile-level glacier change trends and amplitudes to GIS-compatible files [13 s]
GlacierHydrology.gembfit_dv2gpkg(binned_synthesized_dv_file_ref; 
    outfile_prefix="Gardner2025_geotiles_rates", 
    datelimits=(DateTime(2000, 3, 1), DateTime(2025, 1, 1))
)

# 17. Route glacier runoff through river networks for glacier flux calculation [5 min]
# The first run for an inventory also finds each glacier's lowest point on the 30 m DEM and traces its
# river reaches; both are cached per inventory (`GlacierHydrology.glacier_routing_file`, `GlacierHydrology.glacier_rivers_file`).
GlacierHydrology.glacier_routing(; surface_mask=surface_mask_ref);

# 18. Route land surface model runoff through the glacier-fed river reaches of step 17
include("land_surface_model_routing.jl")

# 19. Calculate gmax (maximum glacier contribution to river flux) [3.4 min]
include("gmax_global.jl")

# 20. Generate point-based figures for gmax visualization [1 min]
include("gmax_point_figure.jl")

# 21. Produce regional results for analysis and sharing [7 min]
include("regional_results.jl")

# 22. Generate extended data figures for manuscript
include("manuscript_extended_data_figures.jl")

# 23. Analyze population affected by glacier-fed river changes using buffer analysis [~5hrs]
include("river_buffer_population.jl")

# 24. Create gridded output files for the manuscript
include("gridded_output_files.jl")

# run on local machine to copy over figures to my computer
#=
rsync -rav devon:/mnt/bylot-r3/altim_figs/ /Users/gardnera/Research/24_01_GlobalGlacier/figures/
rsync -rav devon:/mnt/bylot-r3/data/project_data/ ~/data/GlobalGlacier/project_data/
=#
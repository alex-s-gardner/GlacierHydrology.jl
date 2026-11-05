#=
gridded_output_files.jl

Creates gridded NetCDF output files of glacier mass flux data at multiple resolutions.

This script:
1. Loads synthesized geotile data from multiple analysis runs
2. Calculates ensemble means and errors across different processing configurations
3. Downscales data from geotiles to regular lat/lon grids (0.5° and 2° resolution)
4. Converts mass balance variables (runoff, melt, refreeze, SMB, etc.) from kg/m² to kg
5. Saves gridded monthly time series with metadata to NetCDF format
6. Performs sanity checks to verify data consistency

Output files: Gardner2025_glacier_{resolution}deg.nc
=#

# Import required modules and packages for analysis and plotting
begin
    import GlobalGlacierAnalysis as GGA
    import GlacierVolumeChange as GVC
    using Unitful                # Physical units support
    Unitful.register(GGA.MyUnits) # Register custom units from GGA
    using Dates                  # Date and time handling
    using DimensionalData        # For working with labeled dimensions
    using Rasters                # Raster data support
    using Statistics             # Statistical functions
    using CairoMakie             # Plotting library
    using LsqFit                 # Least-squares fitting


    error_quantile = 0.95
    error_scaling = 1.5

    project_id = :v01
    geotile_width = 2

end

# User-defined parameters

# Synthesis parameters of input data
for geotile_width_out in [0.5, 2]
begin
    path2runs_filled, params = GVC.binned_filled_filepaths(;
        project_id,
        surface_masks=["glacier", "glacier_rgi7"],
        dem_ids=["best", "cop30_v2"],
        curvature_corrects=[false, true],
        amplitude_corrects=[true],
        binning_methods=["median", "nmad3", "nmad5"],
        fill_params=[1, 2, 3, 4],
        binned_folders=[GVC.analysis_paths(; geotile_width).binned, replace(GVC.analysis_paths(; geotile_width).binned, "binned" => "binned_unfiltered")],
        include_existing_files_only=true
    )

    path2runs_synthesized = replace.(path2runs_filled, "aligned.jld2" => "synthesized.jld2")

    reference_ensemble_file = GVC.reference_ensemble_file
        
    binned_synthesized_dv_files = GGA.gemb_physics_path.(replace.(path2runs_synthesized, ".jld2" => "_gembfit_dv.jld2"),
                                                         GGA.gemb_reference_physics_member)
    binned_synthesized_dv_file_ref = GGA.gemb_physics_path(replace(reference_ensemble_file, "_aligned.jld2" => "_synthesized_gembfit_dv.jld2"),
                                                           GGA.gemb_reference_physics_member)
        
    path2ref = binned_synthesized_dv_file_ref
    path2files = setdiff(binned_synthesized_dv_files, [binned_synthesized_dv_file_ref])


    # set to be start of record for instantanious (not cumulative) error calculation
    reference_period = (DateTime(1979, 1, 1), DateTime(1979, 12, 31))

    if geotile_width < geotile_width_out
        error("geotile_width_out must be smaller or equal to geotile_width")
    end

    # load in raw synthesis data with errors 
    geotiles0 = GGA.geotiles_mean_error(path2ref, path2files; error_quantile, error_scaling, reference_period) # returns variables in units of kg/m2 [90s]

    geotiles_original = GGA.geotiles_mask_hyps("glacier", geotile_width)
    geotiles_original[!, :glacier_area_km2] = sum.(geotiles_original[!, :glacier_area_km2])
end

# rearrange from geotiles to gridded variables, and multiply by glacier area to give units of monthly mass flux [kg]
begin
    # create global raster for mask coverage calculation
    out_var = :glacier_frac
    geotiles_out = GVC.geotiles_w_mask(geotile_width_out)

    lond = X(DimensionalData.Dimensions.Lookups.Sampled(-180+geotile_width_out/2:geotile_width_out:180-geotile_width_out/2, sampling=DimensionalData.Dimensions.Lookups.Intervals(DimensionalData.Dimensions.Lookups.Center()); metadata =  Dict("long_name" => "longitude", "units" => "degrees")))
    latd = Y(DimensionalData.Dimensions.Lookups.Sampled(-90+geotile_width_out/2:geotile_width_out:90-geotile_width_out/2, sampling=DimensionalData.Dimensions.Lookups.Intervals(DimensionalData.Dimensions.Lookups.Center()); metadata = Dict("long_name" => "latitude", "units" => "degrees")))
    tid = dims(geotiles0, Ti)
    geotiled = dims(geotiles0, :geotile)

    # Calculate cell area for each cell [m^2]
    glacier_area = Raster(zeros(lond, latd); crs = GGA.GeoFormatTypes.EPSG(4326))
    cell_areas = Rasters.cellarea(glacier_area)

    for gt in eachrow(geotiles_out)
        glacier_area[X=gt.extent.X[1] .. gt.extent.X[2], Y=gt.extent.Y[1] .. gt.extent.Y[2]] .= gt.glacier_frac
    end
    
    glacier_area = glacier_area .* cell_areas # area in m^2

    # check that the glacier area is the same as the original geotiles
    area_original = sum(geotiles_original[!, :glacier_area_km2])
    area_downscaled = sum(glacier_area) * 1E-6
    @assert (abs(area_downscaled - area_original) / area_original) < 2E-3

    vars = collect(val(dims(geotiles0, :varname)))

    names_var = Symbol[]
    values_var = [] # area in kg
    
    # create variable level metadata dictionary
    var_metadata = Dict()
    for var0 in vars
        if var0 == "runoff"
            var_metadata[var0] = Dict("long_name" => "runoff", "description" => "total monthly glacier runoff", "units" => "kg")
        elseif var0 == "melt"
            var_metadata[var0] = Dict("long_name" => "melt", "description" => "total monthly glacier melt", "units" => "kg")
        elseif var0 == "refreeze"
            var_metadata[var0] = Dict("long_name" => "meltwater refreeze", "description" => "total monthly glacier meltwater refreeze", "units" => "kg")
        elseif var0 == "rain"
            var_metadata[var0] = Dict("long_name" => "rain", "description" => "total monthly rain on glacier", "units" => "kg")
        elseif var0 == "smb"
            var_metadata[var0] = Dict("long_name" => "surface mass balance", "description" => "total monthly glacier surface mass balance", "units" => "kg")
        elseif var0 == "acc"
            var_metadata[var0] = Dict("long_name" => "accumulation", "description" => "total monthly accumulation on glacier", "units" => "kg")
        elseif var0 == "dm"
            var_metadata[var0] = Dict("long_name" => "mass change", "description" => "total monthly glacier mass change", "units" => "kg")
        elseif var0 == "ec"
            var_metadata[var0] = Dict("long_name" => "evaporation and condensation", "description" => "total monthly glacier evaporation, condensation, sublimation, and deposition", "units" => "kg")
        elseif var0 == "fac"
            var_metadata[var0] = Dict("long_name" => "firn air content", "description" => "change in monthly firn air content", "units" => "mm")
        end
    end

    # add error variable metadata
    for var0 in vars
        meta_copy = copy(var_metadata[var0])
        meta_copy["long_name"] = meta_copy["long_name"] * " error"
        meta_copy["description"] = "error (2σ, fully correlated) in " * meta_copy["description"]
        var_metadata[var0*"_error"] = meta_copy
    end

    
    var_metadata["area"] = Dict("long_name" => "glacier area", "description" => "glacier area", "units" => "m^2")

    global_metadata = Dict(
        "source" => "Alex Gardner/JPL, alex.s.gardner@jpl.nasa.gov",
        "title" => "gridded glacier monthly mass flux",
        "citation" => "Gardner, A. S., Schlegel, N.-J., Greene, C. A., Hugonnet, R., Menounos, B., Wiese, D. N., & Berthier, E. (in review). Glacier contributions to rivers and oceans in the early twenty first century.",
        "version" => "final - " * Dates.format(now(), "yyyy-mm-dd")
    )

    for var0 in vars
        _eltype = eltype(ustrip(geotiles0[varname = At(var0), error = At(false)]))

        push!(names_var, Symbol(var0))
        push!(values_var, zeros(_eltype, lond, latd, tid; metadata = var_metadata[var0]))

        push!(names_var, Symbol(var0 * "_error"))
       
        push!(values_var, zeros(_eltype, lond, latd; metadata = var_metadata[var0*"_error"]))

    end

    # add glacier area variable
    push!(names_var, Symbol("area"))
    push!(values_var, DimArray(parent(glacier_area), (lond, latd); metadata = var_metadata["area"]))
    
    # Bulid DimStack
    ds =  DimStack((; zip(names_var, values_var)...); metadata = global_metadata)
    N = length(tid)

    # downscale data to new resolution
    for var0 in vars
        #var0 = first(vars)
        for gtid in geotiled
            #gtid = geotiled[300]

            ext = GGA.geotile_extent(gtid)

            foo = ustrip(geotiles0[varname=At(var0), geotile=At(gtid), error=At(false)])

            # is is equivelent to foo[2] = foo[2] - foo[1]
            foo[2:end] = diff(foo, dims=Ti)
            
            # set first value equal to second value as we have one less value with diff
            foo[1] = foo[2]

            f =  @view ds[Symbol(var0)][X=ext.X[1] .. ext.X[2], Y=ext.Y[1] .. ext.Y[2]]
            for sl in  eachslice(f, dims=(X,Y))
                sl[:] = foo
            end
            
            # assume error is fully correlated with the signal (after some exploration this seems to be the best error model)
            ds[Symbol(var0 * "_error")][X=ext.X[1] .. ext.X[2], Y=ext.Y[1] .. ext.Y[2]] .= ustrip(geotiles0[varname=At(var0), geotile=At(gtid), error=At(true)])[end] / N
        
        end
    end

    for k in keys(ds)
        if k == :area
            continue
        elseif hasdim(ds[k], :Ti)
            for var in eachslice(ds[k], dims=Ti)
                var .*= glacier_area # [kg/m2 * m2 = kg]
            end
        else
            ds[k] .*= glacier_area
        end
    end
end

begin
    # save output
    path2outfile = joinpath(GGA.pathlocal[:project_dir], "Gardner2025_glacier_$(geotile_width_out)deg.nc")

    # save to disk
    GGA.dimstack2netcdf(ds, path2outfile; deflatelevel=2, chunksizes=[180, 90, 64])
end


# sanity check
begin
    path2outfile = joinpath(GGA.pathlocal[:project_dir], "Gardner2025_glacier_$(geotile_width_out)deg.nc")
    ds = GGA.netcdf2dimstack(path2outfile)

    runoff = ds[:runoff][:,:,:] :: DimArray{Float64,3}
    global_runoff = dropdims(sum(runoff, dims=(:X, :Y)), dims=(:X, :Y), )

    # slope is per year
    kg2Gt = 1e-12
    runoff_rate_downscaled = GGA.calculate_slope(cumsum(global_runoff, dims=:Ti)) * kg2Gt # Gt/yr

    gts = GGA.geotiles_mask_hyps("glacier", geotile_width)
    (ia, ib) = GGA.intersectindices(collect(dims(geotiles0, :geotile)), gts.id)

    foo = geotiles0[varname=At("runoff"), error=At(false)]*1E6 .* sum.(gts[ib, :glacier_area_km2])    
    foo = ustrip(dropdims(sum(foo; dims=:geotile), dims=:geotile))

    runoff_rate_original = GGA.calculate_slope(foo) * kg2Gt # Gt/yr

    if (abs(runoff_rate_downscaled - runoff_rate_original) / runoff_rate_original) < 1.5E-3
        printstyled("Average runoff rate of downscaled and original data are within 0.15% of each other\n"; color=:light_green)
    else
        error("Average runoff rate of downscaled and original data are not within 0.15% of each other")
    end
end
end


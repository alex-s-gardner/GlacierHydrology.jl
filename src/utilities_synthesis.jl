"""
    geotile_grouping!(geotiles0, glaciers, min_area_km2; geotile_groups_manual=nothing)

Group geotiles that share large glaciers and apply manual grouping overrides.

# Arguments
- `geotiles0`: DataFrame with geotile information
- `glaciers`: DataFrame with glacier geometries and areas
- `min_area_km2`: Minimum glacier area threshold for grouping
- `geotile_groups_manual`: Optional manual grouping specifications

# Returns
- Modified DataFrame with group assignments and rectangular geometries

# Examples
```julia
julia> geotile_grouping!(geotiles0, glaciers, 100; geotile_groups_manual=geotile_groups_forced())
julia> # geotiles0 now has :group column
```
"""
function geotile_grouping!(geotiles0, glaciers, min_area_km2; geotile_groups_manual=nothing)
    # Initialize columns to store discharge indices and geotile intersections
    geotiles0[!, :glacier_overap_ind] .= [Int64[]]
    geotiles0[!, :geotile_intersect] .= [Int64[]]

    if "Area" in names(glaciers)
        glaciers_large = glaciers[glaciers.Area.>min_area_km2, :]
    else
        glaciers_large = glaciers[glaciers.area_km2.>min_area_km2, :]
    end
    geometry_column = first(GI.geometrycolumns(glaciers_large))

    # Build spatial index tree for efficient intersection testing
    tree = SortTileRecursiveTree.STRtree(glaciers_large[:, geometry_column]; nodecapacity=3)

    # For each geotile, find intersecting glaciers
    @showprogress dt = 1 desc = "Match glaciers_large with geotiles ..." Threads.@threads for gt in eachrow(geotiles0)
        # Query tree for potential intersecting glaciers
        potential_idxs = SortTileRecursiveTree.query(tree, gt.extent)

        # Find glaciers that actually intersect the geotile
        intersecting = findall(Base.Fix1(GO.intersects, extent2rectangle(gt.extent)),
            view(glaciers_large[:, geometry_column], potential_idxs))

        gt.glacier_overap_ind = potential_idxs[intersecting]
    end

    # Group geotiles that share glaciers
    @showprogress dt = 1 desc = "Group geotiles by glaciers_large ..." for (i, gt) in enumerate(eachrow(geotiles0))
        if .!isempty(gt.glacier_overap_ind)
            # Find other geotiles that share glaciers with this one
            intersecting_geotiles = .!isdisjoint.(Ref(gt.glacier_overap_ind), geotiles0.glacier_overap_ind)
            intersecting_geotiles[i] = false

            if any(intersecting_geotiles)
                gt.geotile_intersect = findall(intersecting_geotiles)
            end
        end
    end

    # Identify connected groups of geotiles
    connectivity = geotiles0.geotile_intersect
    geotiles0[!, :group] = connected_groups(geotiles0.geotile_intersect)


    # Assign new group numbers to manual overrides
    if !isnothing(geotile_groups_manual)
        group0 = maximum(geotiles0.group)
        for grp in geotile_groups_manual
            for g in grp
                geotiles0[findfirst(geotiles0.id .== g), :group] = group0
            end
            group0 += 1
        end
    end

    # Convert geotile extents to rectangles 
    geometry_column = first(GI.geometrycolumns(geotiles0))
    geotiles0[!, geometry_column] = extent2rectangle.(geotiles0.extent)

    return geotiles0[:, Not(:extent, :glacier_overap_ind, :geotile_intersect)]
end

"""
    geotile_grouping(; surface_mask="glacier", min_area_km2=100, geotile_width=2, force_remake_before=nothing)

Group geotiles based on glacier connectivity and save results to file.

# Arguments
- `surface_mask`: Surface mask type to use for filtering (default: "glacier")
- `min_area_km2`: Minimum glacier area threshold for grouping (default: 100)
- `geotile_width`: Width of geotiles in degrees (default: 2)
- `force_remake_before`: Date before which to force recalculation (default: nothing)

# Returns
- `geotiles_out`: GeoDataFrame containing geotile groups with columns :geometry, :id, :group

# Description
This function identifies groups of geotiles that are connected by large overlapping glaciers.
It either loads existing grouping results from file or calculates new groupings by:
1. Filtering geotiles to only include those with glacier coverage
2. Loading glacier geometry data
3. Calling geotile_grouping! to perform the actual grouping logic
4. Saving results to an Arrow file for future use

The function uses caching to avoid recalculating groupings unless the output file
doesn't exist or is older than the force_remake_before date.

# Examples
```julia
julia> geotiles_out = geotile_grouping(; surface_mask="glacier", min_area_km2=100)
julia> # GeoDataFrame with :geometry, :id, :group
```
"""
function geotile_grouping(; surface_mask="glacier", min_area_km2=100, geotile_width=2, force_remake_before=nothing)

    geotile_groups_fn = joinpath(pathlocal[:data_dir], "project_data", "geotile_groups_$(surface_mask)_$(min_area_km2).arrow")
    if !isfile(geotile_groups_fn) || (!isnothing(force_remake_before) && Dates.unix2datetime(mtime(geotile_groups_fn)) < force_remake_before)
        geotiles0 = geotiles_mask_hyps(surface_mask, geotile_width)
        geotiles0[!, :area_km2] = sum.(geotiles0[:, "$(surface_mask)_area_km2"])

        geometry0 = GeoDataFrames.read(pathlocal[Symbol("$(surface_mask)_individual")])

        # identify groups of geotiles that are connected by large overlapping glaciers
        geotiles0 = geotile_grouping!(geotiles0, geometry0, min_area_km2; geotile_groups_manual=geotile_groups_forced())

        geotiles_out = geotiles0[:, [:geometry, :id, :group]]

        # Write geotile groups to file, excluding some columns
        GeoDataFrames.write(geotile_groups_fn, geotiles_out; crs=GFT.EPSG(4326))
    else
        geotiles_out = GeoDataFrames.read(geotile_groups_fn)
    end

    return geotiles_out

end


"""
    discharge2smb(glaciers; discharge2smb_max_latitude=-60, discharge2smb_equilibrium_period=(Date(1979), Date(2000)))

Calculate discharge from surface mass balance (SMB) trends for glaciers below a specified latitude.

# Arguments
- `glaciers`: DataFrame containing glacier data with columns :CenLat, :CenLon, :smb, and :area_km2.
  `smb` is a cumulative surface mass balance DimArray over `:date` in metres of ice equivalent (at
  `δice`); `area_km2` holds the glacier's area in km², summed over its elements.
- `discharge2smb_max_latitude`: Maximum latitude threshold for calculating discharge (default: -60)
- `discharge2smb_equilibrium_period`: Time period for equilibrium calculation (default: 1979-2000)

# Returns
- `discharge0`: DataFrame with columns :latitude, :longitude, :discharge_gtyr, :discharge_err_gtyr, :frontal_ablation_gtyr

# Description
This function estimates glacier discharge by calculating SMB trends during an equilibrium period.
For glaciers below the specified latitude threshold, it fits a linear trend to the SMB data
and, assuming the glacier is in balance, takes discharge as that trend times the glacier area,
converted from km³ ice equivalent to Gt with `δice`.

# Examples
```julia
julia> discharge0 = discharge2smb(glaciers; discharge2smb_max_latitude=-60)
julia> # DataFrame with latitude, longitude, discharge_gtyr
```
"""
function discharge2smb(glaciers; discharge2smb_max_latitude=-60, discharge2smb_equilibrium_period=(Date(1979), Date(2000)))
    index_discharge2smb = glaciers.CenLat .< discharge2smb_max_latitude

    discharge0 = DataFrame(latitude=glaciers[index_discharge2smb, :CenLat], longitude=glaciers[index_discharge2smb, :CenLon], discharge_gtyr=NaN, discharge_err_gtyr=NaN, frontal_ablation_gtyr=NaN)

    ddate = dims(glaciers[1, :smb], :date)
    decyear = decimalyear.(ddate)
    Δdecyear = decyear .- decyear[1]
    index_date = (ddate .>= discharge2smb_equilibrium_period[1]) .& (ddate .<= discharge2smb_equilibrium_period[2])

    for (i, glacier) in enumerate(eachrow(glaciers[index_discharge2smb, :]))
        y = glacier.smb[index_date]
        fit = curve_fit(offset_trend, Δdecyear[index_date], y .- mean(y), offset_trend_p)
        # m i.e. yr⁻¹ × km² / 1000 = km³ i.e. yr⁻¹; × δice / 1000 = Gt yr⁻¹
        discharge0[i, :discharge_gtyr] = fit.param[2] * sum(glacier.area_km2) / 1000 * δice / 1000
    end
    return discharge0
end



"""
    discharge_global_file(surface_mask) -> String

Path of the per-glacier discharge table (Gt/yr) that [`global_discharge_filled`](@ref) writes for
`surface_mask`. South of 60°S the table is built from that mask's glaciers.
"""
discharge_global_file(surface_mask) = replace(pathlocal[:discharge_global], ".jld2" => "_$(surface_mask).jld2")

"""
    global_discharge_filled(;
        surface_mask="glacier",
        discharge_global_fn=discharge_global_file(surface_mask),
        gemb=DimStack(),
        discharge2smb_max_latitude=-60,
        discharge2smb_equilibrium_period=(Date(1979), Date(2000)),
        pscale=1,
        ΔT=0,
        geotile_width=2,
        force_remake_before=nothing,
        force_remake_before_hypsometry=nothing
    )

Compute global glacier discharge by merging observed and modeled values, and return it per geotile.

This function assembles a best-estimate, comprehensive global discharge dataset by combining:
1. Measured glacier discharge (from direct observations where available)
2. Estimated discharge for all other glaciers based on GEMB-modeled surface mass balance (SMB)

The routine performs the following steps:
- Loads a saved discharge file if it exists and is up to date;
- Otherwise, loads geotiles and identifies those with glacier coverage,
- Loads SMB results from a GEMB run, applying requested precipitation and temperature scaling,
- Aggregates GEMB SMB timeseries for each glacier, area-weighted by the glacier's share of its geotile,
- Converts SMB to discharge for unmeasured glaciers using a specified empirical relationship,
- Concatenates the observed and estimated discharge values,
- Ensures no negative values (clips at zero) and that Antarctic total discharge is within a plausible range,
- Saves the result for future fast loading.

# Arguments
- `surface_mask`: (String) Geotile/glacier mask type (e.g. "glacier"). Default: "glacier".
- `discharge_global_fn`: (String) Path to read/write the combined discharge file. Default:
  `discharge_global_file(surface_mask)`. The file records `surface_mask`, and is rejected if it was
  written for another.
- `gemb`: The GEMB output DimStack from `gemb_ensemble_dv(; gemb_run_id, surface_mask)`; it must be built
  on `surface_mask`.
- `discharge2smb_max_latitude`: (Float) Only glaciers south of this latitude have discharge estimated via SMB-to-discharge when not observed. Default: -60.
- `discharge2smb_equilibrium_period`: (Tuple{Date,Date}) Date range for SMB equilibrium (long-term mean). Default: (1979, 2000).
- `pscale`: (Real) Precipitation scaling factor used for GEMB data extraction. Default: 1.
- `ΔT`: (Real) Elevation offset (m) to sample GEMB data at when extracting timeseries. Default: 1.
- `geotile_width`: (Integer) Geotile side length in degrees. Default: 2.
- `force_remake_before`: (Union{Nothing,DateTime}) Force regeneration if file is older than this. Default: nothing.
- `force_remake_before_hypsometry`: (Union{Nothing,DateTime}) Like above, but only for internal DEM/hypsometry products.

# Returns
- DimStack over the GEMB geotiles with one layer per discharge column (`discharge`, `discharge_err`,
  `frontal_ablation`), in km³ ice equivalent per year (`discharge2geotile(...; mass2volume=true)`),
  with `"surface_mask"` in its metadata.
  The file at `discharge_global_fn` holds the per-glacier table in Gt/yr.

# Details
- For glaciers with direct discharge measurements, those values are used.
- For unmeasured glaciers, discharge is estimated by downscaling GEMB-modeled SMB based on glacier area and converting SMB to discharge with a region-dependent empirical relation.
- Results are clipped so that discharge cannot be negative.
- Antarctic discharge sanity check: final total must be 80–115 Gt/yr, otherwise an error is raised.

# Examples
```julia
julia> discharge = global_discharge_filled(; surface_mask="glacier", geotile_width=2)
julia> # DimStack or DataFrame with discharge per geotile
```

# Notes
- The per-glacier table is in Gt/yr; the returned geotile DimStack is in km³ ice equivalent per year.
- The function can be forced to remake its output if `force_remake_before` is newer than the last save time.
- Downscaling to individual glaciers uses area-weighting.
"""
function global_discharge_filled(;
    surface_mask="glacier",
    discharge_global_fn=discharge_global_file(surface_mask),
    gemb = DimStack(),
    discharge2smb_max_latitude=-60,
    discharge2smb_equilibrium_period=(Date(1979), Date(2000)),
    pscale=1,
    ΔT=0,
    geotile_width=2,
    force_remake_before=nothing,
    force_remake_before_hypsometry=nothing
)

    surface_mask = String(surface_mask)
    _check_surface_mask(gemb[:smb], surface_mask, "gemb")

    if file_is_current(discharge_global_fn, force_remake_before)
        file_mask, discharge = FileIO.load(discharge_global_fn, "surface_mask", "discharge")
        file_mask == surface_mask ||
            error("$(discharge_global_fn) holds discharge for surface mask $(file_mask), not $(surface_mask)")
    else
        # Get geotiles containing glaciers and initialize parameters:
        # - Filter to only tiles with glaciers
        # - Set precipitation scaling to 1.0 (no scaling)
        # - Set height adjustment to 0

        # Load geotiles with glacier mask and filter to only include tiles with glacier coverage
        geotiles = _geotile_load_align(;
            surface_mask,
            geotile_width,
            only_geotiles_w_area_gt_0=true
        )

        smb = gemb[:smb][pscale=At(pscale), ΔT=At(ΔT)]

        surface_mask_hypsometry = geotile_hypsometry(geotiles, surface_mask; dem_id=:cop30_v2, force_remake_before=force_remake_before_hypsometry)

        ddate = dims(smb, :date)
        surface_mask_hypsometry[!, "smb"] = [fill(0.0, ddate) for _ in 1:nrow(surface_mask_hypsometry)]

        Threads.@threads for geotile in geotiles.id
            gindex = findall(geotile .== surface_mask_hypsometry.geotile)
            
            if isempty(gindex)
                continue
            end

            gt_dv = smb[geotile=At(geotile)]

            if all(isnan.(gt_dv))
                continue
            end

            garea = sum.(surface_mask_hypsometry[gindex, :area_km2])
            if sum(garea) == 0
                continue
            end

            gweighting = garea ./ sum(garea)
            for (i, ig) in enumerate(gindex)
                if garea[i] == 0
                    continue
                end
                # km³ i.e. over the geotile → m i.e. over this glacier
                surface_mask_hypsometry[ig, "smb"][:] = gt_dv * gweighting[i] ./ garea[i] .* 1000
            end
        end

        # For unmeasured glaciers, estimate discharge using surface mass balance
        # Only process glaciers with non-zero area
        discharge0 = discharge2smb(
            surface_mask_hypsometry;
            discharge2smb_max_latitude,
            discharge2smb_equilibrium_period
        )

        # Combine estimated discharge with measured discharge data
        discharge = glacier_discharge()
        discharge = vcat(discharge, discharge0)

        # Set any negative discharge values to zero
        discharge[discharge.discharge_gtyr.<0, :discharge_gtyr] .= 0

        # Antarctic discharge is a plausibility gate, not a fitted quantity: south of 60S nothing is
        # measured, so every value here is the GEMB SMB trend over the equilibrium period and the total
        # tracks whatever forcing the ensemble was built from.
        antarctic_discharge = sum(discharge[discharge.latitude.<-60, :discharge_gtyr])
        if antarctic_discharge < 80 || antarctic_discharge > 115
            error("Antarctic discharge =  $antarctic_discharge Gt/yr, should be between 80 and 115 Gt/yr")
        end

        # Save the combined measured and estimated discharge data
        save(discharge_global_fn, Dict("discharge" => discharge, "surface_mask" => surface_mask))
    end

    discharge = discharge2geotile(discharge, dims(gemb, :geotile); mass2volume=true)
    DD.metadata(discharge)["surface_mask"] = surface_mask

    return discharge
end



"""
    geotile_synthesis_gembfit_dv(path2runs, discharge, gemb; geotile_width, force_remake_before, geotile_grouping_min_feature_area_km2=100)

Synthesize geotile-level data with GEMB fit parameters and compute calibrated timeseries.

This function processes multiple runs to create calibrated geotile-level datasets by:
1. Loading and aligning geotiles for each surface mask.
2. Calculating glacier hypsometry and identifying geotiles with glaciers.
3. Applying GEMB fit parameters (pscale, ΔT) to calibrate variables.
4. Computing derived variables (discharge, volume change, mass change).
5. Averaging certain variables by geotile groups.
6. Adding regional glacier inventory (RGI) identifiers.

# Arguments
- `path2runs::Vector{String}`: Vector of file paths to binned synthesized data files.
- `discharge::AbstractDict`: surface mask (`String` or `Symbol`) => geotile discharge from
  `global_discharge_filled(; surface_mask, ...)`, in km³ ice equivalent per year.
- `gemb::AbstractDict`: surface mask => GEMB DimStack from `gemb_ensemble_dv(; gemb_run_id, surface_mask)`,
  adjusted with `dv_adjust4discharge!`. Each run uses the entries for its own surface mask; a mask
  without an entry, or an entry built on another mask, throws. Every entry must be of the same GEMB
  physics member ([`gemb_physics_member`](@ref)): each run's fit is read from, and its output written
  to, that member's folder ([`gemb_physics_path`](@ref)), and a fit whose provenance names another
  member throws.
- `geotile_width::Int`: Width of geotiles in degrees.
- `force_remake_before::Union{Nothing, DateTime}`: Optional DateTime to force regeneration of files created before this date.
- `geotile_grouping_min_feature_area_km2::Int=100`: Minimum feature area (km²) for geotile grouping.

# Returns
- Nothing. Saves processed data to files with "_gembfit_dv.jld2" suffix.

# Variables Processed
- Input variables: dv_altim, runoff, fac, smb, rain, acc, melt, ec, refreeze
- Derived variables: discharge, dv (volume change), dm (mass change), dm_altim (altimetry mass change)

# Examples
```julia
julia> geotile_synthesis_gembfit_dv(path2runs, discharge, gemb; geotile_width=2, force_remake_before=nothing)
julia> # Writes _gembfit_dv.jld2 files for each run
```
"""
function geotile_synthesis_gembfit_dv(path2runs, discharge::AbstractDict, gemb::AbstractDict; geotile_width, force_remake_before, geotile_grouping_min_feature_area_km2 = 100)
    gemb = Dict(mask => g for (mask, _, g) in _select_by_surface_mask(gemb, path2runs; name="gemb"))
    discharge = Dict(mask => d for (mask, _, d) in _select_by_surface_mask(discharge, path2runs; name="discharge"))
    for mask in keys(gemb)
        _check_surface_mask(gemb[mask][:dv], mask, "gemb[$(repr(mask))]")
        _check_surface_mask(discharge[mask], mask, "discharge[$(repr(mask))]")
    end

    # Every mask's GEMB ensemble is one physics member's, and each run is read from and written to that
    # member's folder, against a fit whose provenance names the same member.
    members = unique(gemb_physics_member(g[:dv]) for g in values(gemb))
    length(members) == 1 ||
        error("the GEMB ensembles passed are of different physics members ($(join(repr.(sort(members; by=string)), ", "))); " *
              "synthesize one member at a time")
    physics_member = only(members)
   
    # Load elevation change data to get geotile dimensions

    # Initialize dictionaries to store geotiles and glaciers for each surface mask
    geotiles = Dict()
    has_glacier = Dict()
    area_km2 = Dict()

    params = binned_filled_fileparts.(path2runs)
    surface_masks = unique(getindex.(params, :surface_mask))


    # Process each surface mask to load aligned geotiles and calculate hypsometry
    for surface_mask in surface_masks

        # Load and align geotiles with the synthesized data dimensions
        geotiles[surface_mask] = _geotile_load_align(;
            surface_mask,
            geotile_width,
            only_geotiles_w_area_gt_0=true
        )

        area_km2[surface_mask] = _geotile_area_km2(; surface_mask, geotile_width)

        has_glacier[surface_mask] = dropdims(sum(area_km2[surface_mask], dims=:height), dims=:height) .> 0
    end

    # Load synthesized data and GEMB fit parameters
    @showprogress desc = "Computing calibrated geotile level data timeseries for all runs..." Threads.@threads for binned_synthesized_file in path2runs

        binned_synthesized_dv_file = gembfit_dv_path(binned_synthesized_file, physics_member)

        if file_is_current(binned_synthesized_dv_file, force_remake_before)
            printstyled("    -> Skipping $(binned_synthesized_dv_file) because it was created after force_remake_before: $force_remake_before\n"; color=:light_green)
            continue
        else
            run_parameters = binned_filled_fileparts(binned_synthesized_file)
            surface_mask = run_parameters.surface_mask
            gembfit_file = gembfit_path(binned_synthesized_file, physics_member)
            _check_gemb_physics_provenance(get(TOML.parsefile(gembfit_provenance_file(gembfit_file)),
                                               "tile_run_parameters", Dict{String,Any}()),
                                           physics_member, gembfit_file)
            gemb_fit = GeoDataFrames.read(gembfit_file)

            geotiles0 = individual_geotile_synthesis_gembfit_dv(binned_synthesized_file, gemb[surface_mask], gemb_fit, discharge[surface_mask], geotiles[surface_mask], area_km2[surface_mask])

            FileIO.save(binned_synthesized_dv_file, Dict("geotiles" => geotiles0))
        end
    end
end

"""
    discharge2geotile(discharge, dgeotile::Dim; mass2volume=false)

Aggregate discharge from a DataFrame to geotiles using a dimension of geotile IDs.

# Arguments
- `discharge`: DataFrame with longitude, latitude, and discharge columns (e.g. discharge_gtyr)
- `dgeotile`: Dim{:geotile} of geotile identifiers
- `mass2volume`: If true, convert mass (Gt/yr) to volume (km³/yr) using ice density (default: false)

# Returns
- DimStack with one layer per discharge variable, dimensioned by geotile (units Gt/yr or km³/yr).

# Examples
```julia
julia> dgeotile = Dim{:geotile}(geotiles.id)
julia> discharge0 = discharge2geotile(discharge, dgeotile; mass2volume=true)
```
"""
function discharge2geotile(discharge, dgeotile::Dim; mass2volume = false)

    if mass2volume
        mass2volume = 1000 / δice
        units = "km3/yr"
    else
        mass2volume = 1.0
        units = "Gt/yr"
    end

    dvarname = setdiff(names(discharge), ["longitude", "latitude", "extent"])
    outvarname = Symbol.(replace.(dvarname, "_gtyr" => ""))
    

    discharge0 = DimStack([zeros(dgeotile; name) for name in outvarname]; metadata=Dict("units" => units))

    for geotile in dgeotile
        geotile_extent0 = geotile_extent(geotile)
        index = within.(Ref(geotile_extent0), discharge.longitude, discharge.latitude)
        if any(index)
            for i in eachindex(dvarname)
                discharge0[outvarname[i]][geotile=At(geotile)] = sum(discharge[index, dvarname[i]]) * mass2volume
            end
        end
    end   
    return discharge0
end

"""
    individual_geotile_synthesis_gembfit_dv(binned_synthesized_file, gemb, gemb_fit, discharge, geotiles0, area_km2)

Compute calibrated geotile-level time series for a single run using GEMB fit parameters.

Combines altimetry-derived volume change with GEMB-scaled SMB, discharge, and derived variables
for each geotile, applying per-geotile pscale and ΔT from gemb_fit.

# Arguments
- `binned_synthesized_file`: Path to the synthesized altimetry file of an ensemble member; its volume
  change is read from the member's volume product (`GlacierVolumeChange.read_volume`)
- `gemb`: GEMB model output (e.g. DimStack from gemb_ensemble_dv)
- `gemb_fit`: DataFrame with per-geotile fit parameters (pscale, ΔT, group)
- `discharge`: Discharge rate aggregated to geotiles in km³ ice equivalent per year, as returned by
  `discharge2geotile(...; mass2volume=true)`
- `geotiles0`: DataFrame of geotiles to process
- `area_km2`: DimArray of glacier area per geotile and height, used to convert meters ice equivalent to km³

# Returns
- DataFrame of geotiles with time series columns (dv_altim, runoff, fac, smb, etc.) and metadata.
  Each column is in its reporting unit, recorded as its `"units"` metadata: `"km3"` for the volumes
  `dv`, `dv_altim` and `fac` ([`is_volume_varname`](@ref)) and `"Gt"` for every mass, the GEMB
  components, `discharge`, `dm` and `dm_altim`. GEMB volumes are km³ of ice equivalent and are
  converted with `δice/1000` Gt per km³, here only.

# Examples
```julia
julia> geotiles0 = individual_geotile_synthesis_gembfit_dv(binned_synthesized_file, gemb, gemb_fit, discharge, geotiles0, area_km2)
```
"""
function individual_geotile_synthesis_gembfit_dv(binned_synthesized_file, gemb, gemb_fit, discharge, geotiles0, area_km2)
    discharge_units = get(DD.metadata(discharge), "units", nothing)
    discharge_units == "km3/yr" ||
        error("discharge must be in km3/yr (discharge2geotile(...; mass2volume=true)), got units = $(repr(discharge_units))")

    # Filter to only include geotiles containing glaciers
    volume2mass = δice / 1000
    geotiles0 = deepcopy(geotiles0) # geotiles0 is modified within the function

    # Add GEMB fit parameters and mass conversion factors
    geotiles0[!, :pscale] .= 1.0
    geotiles0[!, :ΔT] .= 0.0
    geotiles0[!, :mie2cubickm] .= 0.0
    geotiles0[!, :group] .= 0

    dpscale = dims(gemb, :pscale)
    dΔT = dims(gemb, :ΔT)

    dv_gemb0 = gemb_dv_sample(dpscale[1], dΔT[1], gemb[:smb][geotile=At(geotiles0.id[1])])
    ddate_gemb = dims(dv_gemb0, :date)
    Δdecyear = decimalyear.(ddate_gemb) .- decimalyear.(ddate_gemb)[1]
    geotiles0[!, :discharge] .= [fill(0.0, length(ddate_gemb)) for _ in 1:nrow(geotiles0)]

    for geotile in eachrow(geotiles0)
        fit_index = findfirst(isequal(geotile.id), gemb_fit.id)
        if isnothing(fit_index)
            continue
        end
        
        geotile.pscale = gemb_fit[fit_index, :pscale]
        geotile.ΔT = gemb_fit[fit_index, :ΔT]
        geotile.group = gemb_fit[fit_index, :group]
        area_km20 = sum(area_km2[geotile = At(geotile.id)])
        geotile.mie2cubickm = area_km20 / 1000  # Convert meters ice equivalent to cubic kilometers

        # `discharge` is a rate in km³ ice equivalent per year; accumulate it from the first GEMB date
        geotile.discharge = discharge[:discharge][geotile = At(geotile.id)] * Δdecyear
    end

    # add discharge date metadata
    colmetadata!(geotiles0, "discharge", "date", collect(val(ddate_gemb)), style=:note)

    # Process each variable

    ## Populate with altimetry data
    # Special handling for height change data
   
    dv_altim = GVC.read_volume(GVC.volume_product_path(binned_synthesized_file)).dv

    # Initialize time series arrays
    ddate_altim = dims(dv_altim, :date)
    geotiles0[!, :dv_altim] = [fill(0.0, length(ddate_altim)) for _ in 1:nrow(geotiles0)]

    # Add date metadata for time series
    colmetadata!(geotiles0, :dv_altim, "date", collect(val(ddate_altim)), style=:note)

    # Apply GEMB scaling and convert units for each geotile
    dgeotile = dims(dv_altim, :geotile)

    for geotile in eachrow(geotiles0)
        if in(geotile.id, dgeotile)
            geotile[:dv_altim] = dv_altim[geotile=At(geotile.id)]
        end
    end

    ## Polulate with GEMB data
    # Initialize time series arrays
    for k in keys(gemb)

        geotiles0[!, Symbol(k)] = [zeros(length(ddate_gemb)) for _ in 1:nrow(geotiles0)]

        # Add date metadata for time series
        colmetadata!(geotiles0, Symbol(k), "date", collect(val(ddate_gemb)), style=:note)
    end

    for geotile in eachrow(geotiles0)
        # Apply GEMB scaling and convert units for each geotile
        if in(geotile.id, dgeotile)
            for k in keys(gemb)
                geotile[Symbol(k)][:] = gemb_dv_sample(geotile.pscale, geotile.ΔT, gemb[k][geotile=At(geotile.id)]).data
            end
        end
    end

    # modeled dm
    geotiles0[!, :dm] = (geotiles0[:, :dv] .- geotiles0[:, :fac]) .* volume2mass # to units of Gt
    colmetadata!(geotiles0, "dm", "date", collect(val(ddate_gemb)), style=:note)

    # altimetry dm
    geotiles0[!, :dm_altim] = copy(geotiles0[!, :dv_altim])

    decyear_fac = decimalyear.(colmetadata(geotiles0, "fac", "date"))
    decyear_altim = decimalyear.(colmetadata(geotiles0, "dv_altim", "date"))
    for r in eachrow(geotiles0)
        itp = DataInterpolations.LinearInterpolation(r.fac, decyear_fac)
        fac = itp(decyear_altim)
        r.dm_altim = (r.dv_altim .- fac) .* volume2mass # in units of Gt
    end

    colmetadata!(geotiles0, "dm_altim", "date", colmetadata(geotiles0, "dv_altim", "date"), style=:note)

    # `dichage`, `dm` and `dv` must be average by geotile groups as they are not valid for single geotiles
    vars2average = ["discharge", "dv", "dm", "dv_altim", "dm_altim"]
    gdf = groupby(geotiles0, :group)

    for g in gdf
        #g = gdf[4]
        if nrow(g) > 1
            for varn in vars2average
                #varn = "dv_altim"
                foo = zeros(eltype(g[1, varn]), length(g[1, varn]))
                for r in eachrow(g)
                    foo .+= r[varn]
                end
                foo = foo ./ sum(g[:, :mie2cubickm])

                for r in eachrow(g)
                    r[varn] = foo .* r[:mie2cubickm]
                end
            end
        end
    end

    # Every mass in Gt, every volume in km³, after the group averaging (which is linear and in the
    # stored units).
    for v in vcat("dv_altim", string.(collect(keys(gemb))), "discharge", "dm", "dm_altim")
        if is_volume_varname(v)
            colmetadata!(geotiles0, v, "units", "km3", style=:note)
        else
            v in ("dm", "dm_altim") || foreach(x -> x .*= volume2mass, geotiles0[!, v])
            colmetadata!(geotiles0, v, "units", "Gt", style=:note)
        end
    end

    #### now compute regional estimates
    geotiles0[!, :rgiid] .= Int8(0)
    for i = 1:19
        rgi = "rgi$i"
        geotiles0[geotiles0[:, rgi].>0, :rgiid] .= i
        geotiles0 = geotiles0[:, Not(rgi)]
    end


    varnames_all = vcat("dv_altim", String.(keys(gemb))..., "discharge", "dm", "dm_altim")
    # sanity check
    for varname in varnames_all
        if any([all(isnan.(v)) for v in geotiles0[:, varname]])
            println("$(varname) has all NaNs")
        end
    end

    return geotiles0
end
        

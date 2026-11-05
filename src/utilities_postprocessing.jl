"""
    discharge_rgi(path2discharge, path2rgi_regions; fractional_error = 0.15)

Aggregate glacier discharge data by RGI regions and calculate uncertainty.

# Arguments
- `path2discharge`: Path to file containing discharge data
- `path2rgi_regions`: Path to file containing RGI region geometries
- `fractional_error`: Fractional error to apply to discharge values (default: 0.15)

# Returns
- Array with dimensions [region, error] containing discharge values and uncertainties
  for each RGI region, including special aggregations for HMA (98) and global (99)

# Examples
```julia
julia> drgi = discharge_rgi(path2discharge, path2rgi_regions; fractional_error=0.15)
julia> global_discharge = drgi[At(99), At(false)]
```
"""
function discharge_rgi(path2discharge, path2rgi_regions; fractional_error = 0.15)

    drgi = Dim{:rgi}([1:19; 98; 99])
    derror = Dim{:error}([false, true])
    
    discharge = load(path2discharge, "discharge")

    # determine wich region discharge is within
    rgi_regions = GeoDataFrames.read(path2rgi_regions)

    discharge[!, :geometry] = tuple.(discharge.longitude, discharge.latitude)

    discharge = FlexiJoins.innerjoin(
        (discharge, rgi_regions),
        by_pred(:geometry, GO.within, :geometry)
    )

    discharge = rename(discharge, "RGI_CODE" => "rgi" )

    # sum discharge for each region
    discharge_rgi = zeros(drgi, derror; name = "frontal ablation [Gt yr⁻¹]")
    for rgi in drgi
        index = discharge.rgi .== rgi
        discharge_rgi[At(rgi), At(false)] = sum(discharge[index, "discharge_gtyr"])

    end

    # add HMA and global
    discharge_rgi[At(98), At(false)] = sum(discharge_rgi[13..15, At(false)])
    
    # add global
    discharge_rgi[At(99), At(false)] = sum(discharge_rgi[1..19, At(false)])

    # add error (simple fraction of total)
    discharge_rgi[:, At(true)] = discharge_rgi[:, At(false)] * fractional_error
    return discharge_rgi
end

"""
    is_volume_varname(varname) -> Bool

Whether a geotile or regional variable is a volume, reported in km³, rather than a mass, reported in
Gt. Volumes are the variables whose name contains `dv` or `fac`: `dv`, `dv_altim` (`dv_synthesis`
once renamed) and `fac`. Every other variable -- `dm`, `dm_altim`, `smb`, `acc`, `melt`, `refreeze`,
`runoff`, `rain`, `ec`, `discharge`, `net_acc` and the external `dm_grace`, `dm_glambie` -- is a mass.
"""
is_volume_varname(varname) = occursin("dv", varname) || occursin("fac", varname)

# Units tag `individual_geotile_synthesis_gembfit_dv` gives every column it does not convert to Gt.
const _KM3_ICE_EQUIVALENT = "km3 [i.e.]"

# Factor taking one `_gembfit_dv` column to its regional unit: km³ for a volume, Gt for a mass. A mass
# held as km³ of ice equivalent crosses `δice/1000` Gt per km³ here, once.
function _regional_unit_factor(geotiles0, varname, path)
    units = colmetadata(geotiles0, varname, "units", nothing)
    if is_volume_varname(varname)
        units == _KM3_ICE_EQUIVALENT ||
            error("$path: volume column $varname has units $(repr(units)), expected $(repr(_KM3_ICE_EQUIVALENT))")
        return 1.0
    elseif units == "Gt"
        return 1.0
    elseif units == _KM3_ICE_EQUIVALENT
        return δice / 1000
    else
        error("$path: mass column $varname has units $(repr(units)), expected \"Gt\" or $(repr(_KM3_ICE_EQUIVALENT))")
    end
end

"""
    runs2rgi(path2runs, vars2sum=["dv_altim", "runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze", "dv", "dm", "dm_altim", "discharge"])

Aggregate glacier model run data by RGI (Randolph Glacier Inventory) regions.

This function loads data from multiple model runs, aggregates variables by RGI regions,
and creates special aggregations for HMA (regions 13-15) and global (regions 1-19).

Volumes ([`is_volume_varname`](@ref): `dv`, `dv_altim`, `fac`) are summed in km³. Every other variable
is a mass and is summed in Gt: the GEMB components and `discharge`, which the `_gembfit_dv` files hold
in km³ of ice equivalent, are converted with `δice/1000` Gt per km³ on the way in, and `dm` and
`dm_altim` are already Gt. Each column's `"units"` metadata is checked, and an unexpected tag throws.
The unit of each output is in its `"units"` metadata (`"km3"` or `"Gt"`).

Geotiles with no altimetry at all (`dv_altim` NaN at every date) are left out of every sum: they have no
calibrated fit and carry NaN `dv`, `dm` and `discharge`, which would otherwise make their whole region
NaN. They are the sub-km² ice fragments GEMB covers but the binning drops; an error is thrown if those
left out of any run hold 1 km² of glacier or more, and they are listed once per call.

# Arguments
- `path2runs`: Vector of paths to model run files
- `physics_member`: the GEMB physics member whose calibrated synthesis is read
  ([`gemb_physics_path`](@ref)); `gemb_reference_physics_member` by default, `nothing` for a GEMB
  run without members. [`runs2rgi_physics`](@ref) reads several.
- `vars2sum`: List of variable names to aggregate (default includes common glacier variables)

# Returns
- Dictionary containing aggregated data for each variable, with dimensions [run, rgi, date], each
  recording its `units` and `gemb_physics_member` in its metadata

# Examples
```julia
julia> regional_sum = runs2rgi(path2runs; vars2sum=["runoff", "dm"])
julia> runoff_rgi = regional_sum["runoff"]
```
"""
function runs2rgi(path2runs; physics_member=gemb_reference_physics_member, vars2sum=["dv_altim", "runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze", "dv", "dm", "dm_altim", "discharge"])

    # read in example file
    drgi = Dim{:rgi}([1:19; 98; 99])
   
    # drun needs full path as some runs are distiguished by folder not just filename
    drun = Dim{:run}(path2runs)

    regional_sum = Dict() # this needs to be a Dict since variables have different date dimensions

    gembfit_dv_file(path) = gemb_physics_path(replace(path, ".jld2" => "_gembfit_dv.jld2"), physics_member)
    example_data = FileIO.load(gembfit_dv_file(path2runs[1]), "geotiles")
    for varname in vars2sum
        ddate = Dim{:date}(colmetadata(example_data, varname, "date"))
        units = is_volume_varname(varname) ? "km3" : "Gt"
        regional_sum[varname] = DD.rebuild(fill(NaN, drun, drgi, ddate);
                                           metadata=Dict("units" => units, "gemb_physics_member" => physics_member))
    end

    unobserved = [String[] for _ in path2runs]

    Threads.@threads for i in eachindex(path2runs)
        binned_synthesized_file = path2runs[i]
        binned_synthesized_dv_file = gembfit_dv_file(binned_synthesized_file)
        geotiles0 = FileIO.load(binned_synthesized_dv_file, "geotiles")

        no_altim = [all(isnan, x) for x in geotiles0.dv_altim]
        if any(no_altim)
            area = sum(sum.(geotiles0.area_km2[no_altim]))
            area < 1 || error("$(binned_synthesized_dv_file): geotiles without altimetry hold $(round(area; digits=2)) km² " *
                              "of glacier, too much to leave out of the regional sums: $(geotiles0.id[no_altim])")
            unobserved[i] = geotiles0.id[no_altim]
            geotiles0 = geotiles0[.!no_altim, :]
        end

        for varname in vars2sum
            factor = _regional_unit_factor(geotiles0, varname, binned_synthesized_dv_file)
            geotiles0[!, varname] = [v .* factor for v in geotiles0[!, varname]]
        end

        # group by rgi and sum
        geotiles_reg = groupby(geotiles0, :rgiid)

        # sum across timeseries
        regions0 = DataFrames.combine(geotiles_reg, vars2sum .=> Ref ∘ sum; renamecols=false)

        # add a HMA region 
        index = (regions0[:, :rgiid] .<= 15) .& (regions0[:, :rgiid] .>= 13)
        foo = DataFrames.combine(regions0[index, :], vcat("rgiid", vars2sum) .=> Ref ∘ sum; renamecols=false)
        foo.rgiid .= 98
        regions0 = append!(regions0, foo)

        # add a Global region 
        index = (regions0[:, :rgiid] .<= 19) .& (regions0[:, :rgiid] .>= 1)
        foo = DataFrames.combine(regions0[index, :], vcat("rgiid", vars2sum) .=> Ref ∘ sum; renamecols=false)
        foo.rgiid .= 99
        regions0 = append!(regions0, foo)

        for varname in vars2sum
            # reduce all region time series vectors to a matrix
            regional_sum[varname][At(binned_synthesized_file), At(regions0.rgiid), :] = reduce(hcat, regions0[!, varname])'
        end
    end

    left_out = unique(reduce(vcat, unobserved))
    isempty(left_out) || printstyled("    -> runs2rgi: $(length(left_out)) geotiles without altimetry left out of the regional sums " *
        "in $(count(!isempty, unobserved)) of $(length(path2runs)) runs: $(join(left_out, ", "))\n"; color=:light_yellow)

    return regional_sum
end

"""
    rgi_trends(regional_sum::AbstractDict, discharge_rgi, daterange)

Calculate trends, accelerations, amplitudes, and phases for regional glacier variables.

This function fits time series models to various glacier variables across different regions
and model runs, computing key parameters that characterize their temporal evolution.

# Arguments
- `regional_sum`: Dictionary containing time series data for different glacier variables
- `discharge_rgi`: Discharge values by RGI region, used for GSI calculations
- `daterange`: Time period over which to calculate trends

# Returns
- A DimensionalArray with dimensions for variable name, model run, RGI region, and parameter
  (trend, acceleration, amplitude, phase)

The function handles special calculations for derived variables:
- `net_acc`: Calculated as accumulation minus evaporation/condensation
- `gsi`: Glacier sustainability index, calculated as (runoff + discharge) / net_acc

`regional_sum` is in the units [`runs2rgi`](@ref) returns (masses in Gt, volumes in km³) and
`discharge_rgi` in Gt yr⁻¹, so every term of `net_acc` and `gsi` is in Gt.

# Examples
```julia
julia> daterange = DateTime(2000, 3, 1)..DateTime(2024, 12, 15)
julia> region_fit = rgi_trends(regional_sum, discharge_rgi, daterange)
julia> runoff_trends = region_fit[At("runoff"), :, :, At("trend")]
```
"""
function rgi_trends(regional_sum::AbstractDict, discharge_rgi, daterange)
    varnames = vcat(collect(keys(regional_sum)), ["net_acc", "gsi"])
    dparameter = Dim{:parameter}(["trend", "acceleration", "amplitude", "phase"])
    dvarname = Dim{:varname}(varnames)
    drun = dims(regional_sum[varnames[1]], :run)
    drgi = dims(regional_sum[varnames[1]], :rgi)
    region_fit = zeros(dvarname, drun, drgi, dparameter)

    Threads.@threads for binned_synthesized_file in drun

        for varname in setdiff(varnames, ["gsi"])
           
            if varname == "net_acc"
                ddate = dims(regional_sum["runoff"][:,:,minimum(daterange)..maximum(daterange)], :date)
            else
                ddate = dims(regional_sum[varname][:,:,minimum(daterange)..maximum(daterange)], :date)
            end

            x = decimalyear.(ddate)
            x = x .- ceil(mean(x))

            
            for rgi in drgi
                if varname == "net_acc"
                    y = regional_sum["acc"][At(binned_synthesized_file), At(rgi), minimum(daterange)..maximum(daterange)] .- regional_sum["ec"][At(binned_synthesized_file), At(rgi), minimum(daterange)..maximum(daterange)]
                else
                    y = regional_sum[varname][At(binned_synthesized_file), At(rgi), minimum(daterange)..maximum(daterange)]
                end

                # first fit linear trend and seasonal [this creates less spread in linear fit as acceleration fixed to zero]
                dm_fit = curve_fit(offset_trend_seasonal2, x, y, p_offset_trend_seasonal)
                region_fit[At(varname), At(binned_synthesized_file), At(rgi), At("trend")] = dm_fit.param[2]

                # then fit linear trend and acceleration and seasonal
                dm_fit = curve_fit(offset_trend_acceleration_seasonal2, x, y, p3)
                region_fit[At(varname), At(binned_synthesized_file), At(rgi), At("acceleration")] = dm_fit.param[3]
                region_fit[At(varname), At(binned_synthesized_file), At(rgi), At("amplitude")] = hypot(dm_fit.param[4], dm_fit.param[5])

                phase0 = 365.25 * seasonal_peak_fraction(dm_fit.param[4], dm_fit.param[5])
                region_fit[At(varname), At(binned_synthesized_file), At(rgi), At("phase")] = phase0
            end
        end
    end

    varname = "gsi"
    for rgi in drgi
        for param in ["trend", "acceleration"]

            if param == "trend"
                D = discharge_rgi[At(rgi), At(false)];
            else
                D = 0
            end

            region_fit[At(varname), :, At(rgi), At(param)] = ((region_fit[At("runoff"), :, At(rgi), At(param)] .+ D) ./ region_fit[At("net_acc"), :, At(rgi), At(param)])
        end
    end

    return region_fit
end

"""
    mass_closure_residual(region_fits; rgi, parameter="trend", discharge="discharge") -> Float64

`acc - runoff - ec - discharge - dm` for region `rgi` of `region_fits` (dimensions varname, rgi,
parameter, error; the reference values at `error=false` are used), in Gt yr⁻¹ for the trend.

GEMB's surface mass balance is `smb = acc - runoff - ec` (`gemb_add_derived_vars!`, with
`runoff = melt - refreeze`), and the modelled mass change is `dm = smb - discharge`, so the residual is
zero to the fitting tolerance when every term is in Gt ([`runs2rgi`](@ref)) and `discharge` names the
discharge that `dm` subtracts, i.e. the summed geotile `discharge` column. The trend and acceleration
fits are linear in the data, so the fit of the sum is the sum of the fits.
"""
function mass_closure_residual(region_fits; rgi, parameter="trend", discharge="discharge")
    fit(varname) = region_fits[varname=At(varname), rgi=At(rgi), parameter=At(parameter), error=At(false)]
    return fit("acc") - fit("runoff") - fit("ec") - fit(discharge) - fit("dm")
end

"""
    rgi_trends(da::AbstractDimArray, daterange)

Calculate trend, acceleration, amplitude, and phase parameters for each RGI region in a time series.

This function fits temporal models to glacier data for each RGI region within the specified date range.
It performs two curve fits: first a linear trend with seasonal component, then a model with
linear trend, acceleration, and seasonal components.

# Arguments
- `da`: A DimensionalData array with dimensions for RGI region and date
- `daterange`: Date range to analyze (can be specified as a range of DateTime objects)

# Returns
- A DimensionalData array with dimensions for RGI region and parameter, where parameters include
  trend, acceleration, amplitude, and phase of the fitted curves.

# Examples
```julia
julia> daterange = DateTime(2000, 3, 1)..DateTime(2024, 12, 15)
julia> region_fit = rgi_trends(da, daterange)
julia> trends = region_fit[:, At("trend")]
```

# Note
The function first centers the decimal years around their mean to improve numerical stability.
"""
function rgi_trends(da::AbstractDimArray, daterange)

    dparameter = Dim{:parameter}(["trend", "acceleration", "amplitude", "phase"])
    drgi = dims(da, :rgi)

    ddate = dims(da[:,minimum(daterange).. maximum(daterange)], :date)

    d = decimalyear.(ddate)
    d = d .- ceil(mean(d))

    region_fit = zeros(drgi, dparameter)

    for rgi in drgi

        # first fit linear trend and seasonal [this creates less spread in linear fit as acceleration fixed to zero]
        dm_fit = curve_fit(offset_trend_seasonal2, d, da[At(rgi), minimum(daterange)..maximum(daterange)], p_offset_trend_seasonal)
        region_fit[At(rgi), At("trend")] = dm_fit.param[2]

        # then fit linear trend and acceleration and seasonal
        dm_fit = curve_fit(offset_trend_acceleration_seasonal2, d, da[At(rgi), minimum(daterange)..maximum(daterange)], p3)
        region_fit[At(rgi), At("acceleration")] = dm_fit.param[3]
        region_fit[At(rgi), At("amplitude")] = hypot(dm_fit.param[4], dm_fit.param[5])

        phase0 = 365.25 * seasonal_peak_fraction(dm_fit.param[4], dm_fit.param[5])
        region_fit[At(rgi), At("phase")] = phase0
    end

    return region_fit
end

"""
    region_fit_ref_and_err(region_fit, path2reference; p = 0.95, discharge = nothing)

Calculate reference values and error estimates for regional fit parameters.

This function extracts reference values from a specified path and computes error estimates
based on the quantile of absolute differences between all values and the reference value.

# Arguments
- `region_fit`: A DimensionalData array with dimensions for variable name, path, RGI region, and parameter
- `path2reference`: The reference path to use for extracting reference values
- `p`: Quantile level for error estimation (default: 0.95)
- `discharge`: Discharge values with error estimates, required when "net_acc" is present in variables

# Returns
- A DimensionalData array with dimensions for variable name, RGI region, parameter, and error,
  where the error dimension contains the reference values (false) and error estimates (true)

# Examples
```julia
julia> regions_fit = region_fit_ref_and_err(region_fit, path2reference; error_quantile=0.95)
julia> ref_values = regions_fit[:, :, :, At(false)]
julia> err_values = regions_fit[:, :, :, At(true)]
```

# Note
When "net_acc" is present, the function also calculates error propagation for the glacier state index (gsi)
using standard error propagation formulas.
"""
function region_fit_ref_and_err(region_fit, path2reference; error_quantile=0.95, error_scaling=1.0, discharge=nothing)
    derror = Dim{:error}([false, true])
    (dvarname, _, drgi, dparameter) = dims(region_fit)
    regions_fit = fill(NaN, dvarname, drgi, dparameter, derror)
    for varname in dvarname
        for rgi in drgi
            for param in dparameter
                v = region_fit[At(varname), At(path2reference), At(rgi), At(param)]
                regions_fit[At(varname), At(rgi), At(param), At(false)] = v

                err0 = quantile(abs.(region_fit[At(varname), :, At(rgi), At(param)] .- v), error_quantile) * error_scaling
                regions_fit[At(varname), At(rgi), At(param), At(true)] = err0
            end
        end
    end

    if "net_acc" in dvarname
        if isnothing(discharge)
            error("gsi is an exisiting variable, so discharge_error must be provided as a kwarg")
        end

        A = regions_fit[At("runoff"), :, At("trend"), At(false)] .+ discharge[:, At(false)]
        σA = sqrt.(regions_fit[At("runoff"), :, At("trend"), At(true)].^ 2 .+ discharge[:, At(true)] .^ 2)

        B = regions_fit[At("net_acc"), :, At("trend"), At(false)]
        σB = regions_fit[At("net_acc"), :, At("trend"), At(true)]

        err = abs.(regions_fit[At("gsi"), :, At("trend"), At(true)]) .* sqrt.((σA ./ A) .^ 2 .+ (σB ./ B) .^ 2)
        regions_fit[At("gsi"), :, At("trend"), At(true)] = err
    end


    return regions_fit
end


"""
    region_min_frac_error!(region_fits; fractional_error_min, varnames, parameters) -> region_fits

Raise the error of each of `varnames` × `parameters` to at least `fractional_error_min` times the
magnitude of its value: the error becomes `max(ensemble error, fractional_error_min × |value|)`. Modifies
and returns `region_fits` (dimensions `varname`, `rgi`, `parameter`, `error`).
"""
function region_min_frac_error!(region_fits; fractional_error_min, varnames, parameters)
    err = region_fits[varname=At(varnames), parameter = At(parameters), error = At(true)]
    err_frac = abs.(region_fits[varname=At(varnames), parameter = At(parameters), error = At(false)] .* fractional_error_min)
    region_fits[varname=At(varnames), parameter = At(parameters), error = At(true)] = max.(err_frac, err)
    return region_fits
end


"""
    runs_center!(runs_rgi::Dict, daterange)

Center glacier model run data by subtracting the temporal mean for each variable.

# Arguments
- `runs_rgi`: Dictionary containing time series data for different glacier variables
- `daterange`: Time period over which to calculate the mean

# Returns
- The modified `runs_rgi` dictionary with centered data

This function modifies the input dictionary in-place, subtracting the temporal mean
of each variable over the specified date range from all values of that variable.

# Examples
```julia
julia> daterange = DateTime(2000, 3, 1)..DateTime(2024, 12, 15)
julia> runs_center!(runs_rgi, daterange)
```
"""
function runs_center!(runs_rgi::Dict, daterange)
 
    for k in keys(runs_rgi)
        mean_rgi = dropdims(mean(runs_rgi[k][:,:,minimum(daterange)..maximum(daterange)], dims = :date), dims = :date)
        runs_rgi[k] .-= mean_rgi
    end

    return runs_rgi
end

"""
    regions_center!(runs_rgi, daterange)

Center regional glacier data by subtracting the temporal mean.

# Arguments
- `runs_rgi`: DimensionalArray containing regional glacier data
- `daterange`: Time period over which to calculate the mean

# Returns
- The modified `runs_rgi` array with centered data

This function modifies the input array in-place, subtracting the temporal mean
over the specified date range from all values in the array.

# Examples
```julia
julia> daterange = DateTime(2000, 3, 1)..DateTime(2024, 12, 15)
julia> regions_center!(runs_rgi, daterange)
```
"""
function regions_center!(runs_rgi, daterange)
    mean_rgi = dropdims(mean(runs_rgi[:, minimum(daterange)..maximum(daterange), At(false)], dims=:date), dims=:date)
    runs_rgi .-= mean_rgi

    return runs_rgi
end


"""
    runs_delta!(runs_rgi::Dict, path2reference) -> Dict

Calculate differences between each run and a reference run for all variables.

# Arguments
- `runs_rgi`: Dictionary containing time series data for different glacier variables
- `path2reference`: Path to the reference run to subtract from all other runs

# Returns
- The modified `runs_rgi` dictionary with delta values

This function modifies the input dictionary in-place, subtracting the reference run
from all runs for each variable, creating anomalies relative to the reference.

# Examples
```julia
julia> runs_delta!(runs_rgi, path2reference)
julia> # runs_rgi now contains anomalies relative to path2reference
```
"""
function runs_delta!(runs_rgi, path2reference)
    for k in keys(runs_rgi)
        runs_rgi[k] = @d runs_rgi[k] .- runs_rgi[k][At(path2reference), :, :]
    end
    return runs_rgi
end
"""
    runs_quantile(runs_rgi, p; on_abs=true) -> Dict

Calculate quantiles of glacier data across multiple model runs.

# Arguments
- `runs_rgi`: Dictionary containing DimensionalArrays with model run data
- `p`: Quantile to calculate (e.g., 0.95 for 95th percentile)
- `on_abs`: Whether to calculate quantiles on absolute values (default: true)

# Returns
- Dictionary with same keys as input, containing quantile values for each region and date

This function computes quantiles across the run dimension for each variable, region, and date.
When `on_abs=true`, quantiles are calculated on absolute values, useful for error estimation.

# Examples
```julia
julia> runs_q = runs_quantile(runs_rgi, 0.95; on_abs=true)
julia> runoff_95 = runs_q["runoff"]
```
"""
function runs_quantile(runs_rgi, p; on_abs=true)
    runs_quantile = Dict()
    on_abs ? (fun=abs) : (fun = identity)

    for k in keys(runs_rgi)
        (drun, drgi, ddate) = dims(runs_rgi[k])
        runs_quantile[k] = fill(NaN, drgi, ddate)
        for rgi in drgi
            for date in ddate
                v0 = @view runs_rgi[k][:, At(rgi), At(date)]
                if any(isnan.(v0))
                    continue
                else
                    runs_quantile[k][At(rgi), At(date)] = quantile(fun.(v0), p)
                end
            end
        end
    end
    return runs_quantile
end

"""
    runs_select(runs_rgi, runs2select) -> Dict

Extract a subset of model runs from a dictionary of glacier data.

# Arguments
- `runs_rgi`: Dictionary containing DimensionalArrays with model run data
- `runs2select`: Vector of run identifiers to select

# Returns
- Dictionary with same keys as input, containing only the selected runs

This function creates a new dictionary with the same structure as the input,
but including only the specified model runs.

# Examples
```julia
julia> runs_subset = runs_select(runs_rgi, [path2ref, path2run2])
julia> # Dictionary contains only the two specified runs
```
"""
function runs_select(runs_rgi, runs2select)
    runs_selected = Dict()
    for k in keys(runs_rgi)
        runs_selected[k] = runs_rgi[k][At(runs2select), :, :]
    end
    return runs_selected
end


"""
    runs_ref_and_err(runs_rgi, path2reference; p = 0.95) -> Dict

Create a dictionary containing reference values and error estimates for glacier model runs.

# Arguments
- `runs_rgi`: Dictionary containing DimensionalArrays with model run data
- `path2reference`: Path identifier for the reference run
- `p`: Quantile level for error estimation (default: 0.95)

# Returns
- Dictionary with same keys as input, containing DimensionalArrays with an additional
  error dimension where false=reference values and true=error estimates

This function calculates deviations from a reference run, computes quantile-based
error estimates, and combines them into a single data structure with an error dimension.

# Examples
```julia
julia> runs_with_err = runs_ref_and_err(runs_rgi, path2reference; error_quantile=0.95)
julia> ref_runoff = runs_with_err["runoff"][:, :, :, At(false)]
```
"""
function runs_ref_and_err(runs_rgi, path2reference; error_quantile=0.95, error_scaling=1.0)
    # remove reference run
    runs_rgi_delta = runs_delta!(deepcopy(runs_rgi), path2reference)

    # get 95% confidence interval
    regions_err = runs_quantile(runs_rgi_delta, error_quantile; on_abs=true)

    # selected run as reference
    regions_ref = runs_select(runs_rgi, path2reference)

    regions = Dict()
    derror = Dim{:error}([false, true])
    for k in keys(regions_ref)
        regions[k] = zeros(dims(regions_ref[k])..., derror)
        regions[k][:,:,At(false)] = regions_ref[k]
        regions[k][:, :, At(true)] = regions_err[k] * error_scaling
    end

    return regions
end



"""
    rgi_endorheic(path2river_flux, glacier_summary_file; dates4trend=nothing)

Calculate glacier mass change and runoff trends for endorheic and non-endorheic basins by RGI region.

# Arguments
- `path2river_flux`: Path to file containing river flux data with terminal status information
- `glacier_summary_file`: Path to NetCDF file containing glacier runoff and mass change data
- `dates4trend`: Optional date range over which to calculate trends

# Returns
- A DimensionalArray with dimensions for variable name (dm, runoff), RGI region, and endorheic correction
  where false=all glaciers and true=endorheic glaciers only

This function identifies endorheic (landlocked) glaciers by matching glacier IDs to terminal river basins,
calculates mass change and runoff trends, and aggregates results by RGI region including special
aggregations for HMA (region 98) and global (region 99).

# Examples
```julia
julia> glacier_dm_gt = rgi_endorheic(path2river_flux, glacier_summary_file; dates4trend=(DateTime(2000,1,1), DateTime(2024,12,31)))
julia> dm_all = glacier_dm_gt[At("dm"), :, At(false)]
```
"""
function rgi_endorheic(path2river_flux, glacier_summary_file; dates4trend=nothing)
    ds = NCDataset(glacier_summary_file)
    glaciers = DataFrame(rgiid = ds["rgiid"][:], runoff = eachrow(ds["runoff"][:,:]), dm = eachrow(ds["dm"][:,:]))
    DataFrames.colmetadata!(glaciers,"runoff", "date", ds["Ti"][:]; style=:note)
    DataFrames.colmetadata!(glaciers,"dm", "date", ds["Ti"][:]; style=:note)

    
    drgi = Dim{:rgi}([collect(1:19)..., 98, 99])

    # load river flux data
    river_flux = GeoDataFrames.read(path2river_flux)

    # select only terminal rivers
    river_flux = river_flux[river_flux[:, :NextDownID].==0, :]

    # match glaciers to river sink type (endorheic or ocean terminating)
    gid1 = reduce(vcat, river_flux[!, "RGIId"])
    endorheic = falses(length(gid1))

    c = 1
    for r in eachrow(river_flux)
        n = length(r["RGIId"])
        endorheic[c:c+n-1] .= !(r.ocean_terminating)
        c += n
    end

    p1 = sortperm(gid1)
    gid1 = gid1[p1]
    endorheic = endorheic[p1]

    indices = [searchsortedfirst(gid1, x) for x in glaciers[:, "rgiid"]]

    # some glaciers are not in the river_flux file (no rivers in the Antarctic)
    valid_indices = indices .<= length(gid1)

    # add endorheic column to glaciers
    glaciers[!, :endorheic] .= false
    glaciers[valid_indices, :endorheic] .= endorheic[indices[valid_indices]]

    # compute mass change for each region for endorheic and glaciers

    # now compute rates for all glaciers
    varnames = ["dm", "runoff"]  # unit of Gt

    # fit trend to all glaciers
    df_tsfit!(glaciers, varnames; progress=true, datelimits=dates4trend)

    
    # add rgi column
    glaciers[!, :rgi] = [parse(Int16, g[7:8]) for g in glaciers.rgiid]

    # (dh - fac)/1000 * area * volume2mass = dm_gt
    glaciers[!, :dm] = glaciers.dm_trend
    glaciers[!, :runoff] = glaciers.runoff_trend

    # derive statistics for each region
    dvarnames = Dim{:varname}(["dm", "runoff"])
    dendorheic_only = Dim{:endorheic_correction}([false, true])
    glacier_dm_gt = zeros(dvarnames, drgi, dendorheic_only; name="all vs. endorheic only", metadata=Dict("units" => "Gt/yr", "daterange" => dates4trend))

    for rgi in setdiff(collect(drgi), [98, 99])
        index = glaciers.rgi .== rgi
        endorheic = glaciers[:, :endorheic] .& index
        for varname in dvarnames
            glacier_dm_gt[At(varname), At(rgi), At(false)] = sum(glaciers[index, varname])
            glacier_dm_gt[At(varname), At(rgi), At(true)] = sum(glaciers[endorheic, varname])
        end
    end

    # add HMA and global
    glacier_dm_gt[:, At(98), :] = sum(glacier_dm_gt[:, 13..15, :], dims=:rgi)

    # add global
    glacier_dm_gt[:, At(99), :] = sum(glacier_dm_gt[:, 1..19, :], dims=:rgi)

    return glacier_dm_gt
end
"""
    error_bar_table(region_fits, varnames, fit_params, rgi_regions; digits=2)

Create a formatted table of regional glacier parameter values with error estimates.

# Arguments
- `region_fits`: DimensionalArray containing fitted parameters with dimensions for variable name, 
                 RGI region, parameter type, and error
- `varnames`: List of variable names to include in the table
- `fit_params`: List of parameters to include (e.g., "trend", "acceleration", "amplitude", "phase")
- `rgi_regions`: List of RGI region IDs to include in the table
- `digits`: Number of decimal places for rounding numeric values (default: 2)

# Returns
- DataFrame with columns for RGI ID, region name, and formatted parameter values with error estimates
  in the format "value ± error" for each variable and parameter combination

The function automatically determines appropriate units for each parameter type and includes
them in the column headers.

# Examples
```julia
julia> tbl = error_bar_table(region_fits; rgi_regions=[1, 2, 99], digits=2)
julia> # DataFrame with columns for each variable and parameter
```
"""
function error_bar_table(region_fits; varnames=nothing, params=["trend", "acceleration", "amplitude", "phase"], rgi_regions=nothing, digits=2)

    rgi_labels = rginum2label.(rgi_regions)

    regional_results = DataFrames.DataFrame(rgi = rgi_regions, region_name = rgi_labels)

    if isnothing(varnames)
        varnames = val(dims(region_fits, :varname))
    end

    if isnothing(rgi_regions)
        rgi_regions = val(dims(region_fits, :rgi))
    end


    if params == nothing
        params = dims(region_fits, :fit_param)
    end

    ismass = falses(dims(region_fits, :varname))
    for varname in varnames
        ismass[varname = At(varname)] = !is_volume_varname(varname)
    end



    for varname in varnames

        ismass[varname=At(varname)] ? base_unit = "Gt" : base_unit = "km^3"

        #println("    $(varname)")
        for param in params
            if param in ["acceleration"]
                unit = "$(base_unit)/yr^2"
            elseif param in ["trend"]
                unit = "$(base_unit)/yr"
            elseif param in ["amplitude"]
                unit = "$(base_unit)"
            elseif param in ["phase"]
                unit = "day_of_max"
            end

            column_label = Symbol("$(varname)_$(param)_[$unit]")
            regional_results[!, column_label] .= ""
        end
    end

    for (i, rgi) in enumerate(rgi_regions)
        #println("$(rginum2label(rgi))")
        for varname in varnames
            ismass[varname=At(varname)] ? base_unit = "Gt" : base_unit = "km^3"
            #println("    $(varname)")
            for param in params
                v = region_fits[At(varname), At(rgi), At(param), At(false)]
                err = region_fits[At(varname), At(rgi), At(param), At(true)]
                

                if param in ["acceleration"]
                    v = v
                    err = err
                    unit = "$(base_unit)/yr^2"
                elseif param in ["trend"]
                    unit = "$(base_unit)/yr"
                elseif param in ["amplitude"]
                    unit = "$(base_unit)"
                elseif param in ["phase"]
                    unit = "day_of_max"
                end
               
                if param in ["phase"]
                    v = round(Int, v)
                    err = round(Int, err)
                else
                    v = round(v; digits)
                    err = round(err; digits)
                end
                
                column_label = Symbol("$(varname)_$(param)_[$unit]")
                regional_results[i, column_label] = "$(v) ± $(err)"
            end
        end
    end
    return regional_results
end


"""
    region_fits_table(region_fits; varnames=nothing, rgi_regions=nothing, no_err=String[], no_acceleration=String[], amplitude=String[], digits=2)

Create a DataFrame table from region_fits DimArray with trend and acceleration columns for specified variables.

# Arguments
- `region_fits`: DimArray with dimensions (varname, rgi, parameter, error)
- `varnames`: Optional vector of variable names to include, in desired order (default: nothing, includes all variables)
- `rgi_regions`: Optional vector of RGI region IDs to include (default: nothing, includes all regions)
- `no_err`: Vector of variable names to display without error bars (default: empty array)
- `no_acceleration`: Vector of variable names to exclude acceleration column for (default: empty array)
- `amplitude`: Vector of variable names to include amplitude column for (without error) (default: empty array)
- `digits`: Number of decimal places for rounding (default: 2)

# Returns
- DataFrame with rows for each RGI region and columns for each variable's trend and acceleration,
  formatted as "value ± error" strings (unless variable is in no_err list). Amplitude is included
  without error for variables specified in the amplitude list.

# Example
```julia
# Include all variables and regions
result_table = region_fits_table(region_fits)

# Specify variables in desired order
result_table = region_fits_table(region_fits, varnames=["fac", "refreeze", "rain", "acc"])

# Suppress errors for certain variables
result_table = region_fits_table(region_fits, varnames=["gsi", "fac"], no_err=["gsi"])

# Exclude acceleration for certain variables
result_table = region_fits_table(region_fits, varnames=["gsi", "discharge", "fac"], no_acceleration=["gsi", "discharge"])

# Include amplitude column for specific variables
result_table = region_fits_table(region_fits, varnames=["dm", "runoff"], amplitude=["dm"])

# Combine options
result_table = region_fits_table(region_fits,
    varnames=["gsi", "fac", "refreeze"],
    no_err=["gsi"],
    no_acceleration=["gsi"],
    amplitude=["fac"])
```
"""
function region_fits_table(region_fits; varnames=nothing, rgi_regions=nothing, no_err=String[], no_acceleration=String[], amplitude=String[], digits=2)

    # Extract variable names from DimArray if not specified
    if isnothing(varnames)
        varnames = val(dims(region_fits, :varname))
    end

    # Extract RGI regions from DimArray if not specified
    if isnothing(rgi_regions)
        rgi_regions = val(dims(region_fits, :rgi))
    end

    # Create base DataFrame with RGI numbers and region names
    rgi_labels = rginum2label.(rgi_regions)
    df = DataFrames.DataFrame(rgi = rgi_regions, region_name = rgi_labels)

    # Initialize columns for each varname (trend first, then acceleration, then amplitude if requested)
    for varname in varnames
        # Determine unit based on variable type
        base_unit = is_volume_varname(varname) ? "km³" : "Gt"

        # Create column labels
        trend_label = Symbol("$(varname)_trend_[$(base_unit)/yr]")
        df[!, trend_label] .= ""

        # Only create acceleration column if not in no_acceleration list
        if !(varname in no_acceleration)
            accel_label = Symbol("$(varname)_acceleration_[$(base_unit)/yr²]")
            df[!, accel_label] .= ""
        end

        # Create amplitude column if in amplitude list
        if varname in amplitude
            amp_label = Symbol("$(varname)_amplitude_[$(base_unit)]")
            df[!, amp_label] .= ""
        end
    end

    # Populate data for each region and variable
    for (i, rgi) in enumerate(rgi_regions)
        for varname in varnames
            # Determine unit based on variable type
            base_unit = is_volume_varname(varname) ? "km³" : "Gt"

            # Extract and format trend values
            trend_val = region_fits[At(varname), At(rgi), At("trend"), At(false)]
            trend_err = region_fits[At(varname), At(rgi), At("trend"), At(true)]
            trend_val = round(trend_val; digits)
            trend_err = round(trend_err; digits)
            trend_label = Symbol("$(varname)_trend_[$(base_unit)/yr]")

            # Format with or without error based on no_err list
            if varname in no_err
                df[i, trend_label] = "$(trend_val)"
            else
                df[i, trend_label] = "$(trend_val) ± $(trend_err)"
            end

            # Extract and format acceleration values (only if not in no_acceleration list)
            if !(varname in no_acceleration)
                accel_val = region_fits[At(varname), At(rgi), At("acceleration"), At(false)]
                accel_err = region_fits[At(varname), At(rgi), At("acceleration"), At(true)]
                accel_val = round(accel_val; digits)
                accel_err = round(accel_err; digits)
                accel_label = Symbol("$(varname)_acceleration_[$(base_unit)/yr²]")

                # Format with or without error based on no_err list
                if varname in no_err
                    df[i, accel_label] = "$(accel_val)"
                else
                    df[i, accel_label] = "$(accel_val) ± $(accel_err)"
                end
            end

            # Extract and format amplitude values (only if in amplitude list)
            if varname in amplitude
                amp_val = region_fits[At(varname), At(rgi), At("amplitude"), At(false)]
                amp_val = round(amp_val; digits)
                amp_label = Symbol("$(varname)_amplitude_[$(base_unit)]")
                df[i, amp_label] = "$(amp_val)"
            end
        end
    end

    return df
end


"""
    geotile2dimarray_kgm2(binned_synthesized_dv_file; vars2extract=["runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze", "dm"], trim=true, reference_period=nothing)

Convert geotile data to a dimensional array with standardized units (kg/m²).

# Arguments
- `binned_synthesized_dv_file`: Path to the file containing geotile data
- `vars2extract`: List of variable names to extract from the file (default includes common glacier variables)
- `trim`: Whether to trim data to valid range (default: true)
- `reference_period`: Optional time period to use as reference for anomaly calculation

# Returns
- A DimensionalArray with dimensions for variable name, geotile ID, and time, with units in kg/m²
  (except for 'dv' which is in meters)

This function loads geotile data, converts variables to consistent units (kg/m²), and optionally
calculates anomalies relative to a reference period. It handles special unit conversions for
volume ('dv') and mass ('dm') variables.

# Examples
```julia
julia> da = geotile2dimarray_kgm2(binned_synthesized_dv_file; vars2extract=["runoff", "dm"])
julia> runoff_gt = da[At("runoff"), :, :]
```
"""
function geotile2dimarray_kgm2(
    binned_synthesized_dv_file; 
    vars2extract=["runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze", "dm"], 
    trim=true, 
    reference_period=nothing
    )

    geotiles0 = load(binned_synthesized_dv_file, "geotiles")

    geotiles0[!, :area_km2] = sum.(geotiles0.area_km2)

    # trim data to valid range and update date metadata
    if trim
        valid, = validrange(.!isnan.(geotiles0[1, vars2extract[1]]))
    else
        valid = 1:length(geotiles0[1, vars2extract[1]])
    end

    ddate = colmetadata(geotiles0, vars2extract[1], "date")[valid]
    dTi = Ti(ddate)

    varunits = Dict()
    for varname in vars2extract
        if varname == "dv"
            varunits[varname] = u"m"
        else
            varunits[varname] = u"kg/m^2"
        end
    end
    geotile_area = Dict()
    for (id, area) in zip(geotiles0.id, geotiles0.area_km2)
        geotile_area[id] = area * u"km^2"
    end

    dvarname = Dim{:varname}(vars2extract)
    dgeotile = Dim{:geotile}(geotiles0.id)

    # TODO: I've hardcoded the units here to be the units of the first variable in vars2extract ... this will error for mixed units... fixing this would require splitting geotile_out into multiple arrays. 
    geotile_out = fill(NaN * varunits[first(vars2extract)], dvarname, dgeotile, dTi; name=joinpath(splitpath(binned_synthesized_dv_file)[end-2:end]), metadata=Dict("area" => geotile_area))

    # remove unneeded columns
    geotiles0 = geotiles0[:, vcat("geometry", "area_km2", vars2extract)]

    for varname in vars2extract

        if varname == "dv"
            # convert from km3 to m
            geotile_out[At(varname), :, :] = reduce(hcat, getindex.(geotiles0[:, varname], (valid,)) ./ geotiles0[:, :area_km2] * 1000 * 1000 * varunits[varname])'

        elseif varname == "dm"
            # convert from Gt to kg/m²
            geotile_out[At(varname), :, :] = reduce(hcat, getindex.(geotiles0[:, varname], (valid,)) ./ geotiles0[:, :area_km2] * 1000 * 1000 * varunits[varname])'
        else
            # convert from km3 of ice to kg/m²
            geotile_out[At(varname), :, :] = reduce(hcat, getindex.(geotiles0[:, varname], (valid,)) ./ geotiles0[:, :area_km2] * δice * 1000 * varunits[varname])'
        end
    end

    if !isnothing(reference_period)
        geotile_out .-= mean(geotile_out[Ti=reference_period[1]..reference_period[2]], dims=Ti)
    end

    return geotile_out
end


"""
    geotiles_mean_error(path2ref, path2files; error_quantile=0.95, error_scaling=1.0, reference_period=nothing)

Compute reference values and error estimates for geotile-level data across multiple model runs.

This function loads geotile data for a reference run and a set of comparison runs, computes the deviation of each comparison run from the reference, and then estimates the error at each point as the specified quantile of the absolute deviations, optionally scaled. The result is a DimensionalArray with an added error dimension, containing both the reference values and the error estimates.

# Arguments
- `path2ref::AbstractString`: Path to the reference model run file.
- `path2files::AbstractVector{<:AbstractString}`: Vector of paths to model run files to compare against the reference.
- `error_quantile::Real=0.95`: Quantile level for error estimation (e.g., 0.95 for the 95th percentile).
- `error_scaling::Real=1.0`: Optional scaling factor to apply to the error estimate (e.g., 2.0 for 2-sigma).
- `reference_period::Union{Nothing,Tuple}`: Optional time period for centering data (subtracting temporal mean).

# Returns
- `DimArray`: A DimensionalArray with dimensions for variable name, geotile, time, and error,
  where the error dimension contains reference values (`false`) and error estimates (`true`).

# Examples
```julia
julia> geotiles_with_err = geotiles_mean_error(path2ref, path2files; error_quantile=0.95)
julia> ref_runoff = geotiles_with_err[At("runoff"), :, :, At(false)]
```

# Notes
- The error estimate is computed as the specified quantile of the absolute deviation of each run from the reference, for each variable, geotile, and time step.
- The function assumes all input files are compatible and aligned in variables, geotiles, and time.
"""
function geotiles_mean_error(path2ref, path2files; error_quantile=0.95, error_scaling=1.0, reference_period=nothing)

    geotiles_ref = geotile2dimarray_kgm2(path2ref; reference_period)

    # place all runs into array with same number of geotiles... if a geotile is missing then set dv = 0
    drun = Dim{:run}(path2files)
    dgeotile = dims(geotiles_ref, :geotile)
    dvarname = dims(geotiles_ref, :varname)
    dTi = dims(geotiles_ref, :Ti)
    geotiles0 = zeros(drun, dvarname, dgeotile, dTi) * 1u"kg/m^2"

    @showprogress desc = "Loading geotile dv data and calculating error..." for run in drun
        gt = geotile2dimarray_kgm2(run; reference_period)

        # there is an issue here where RGI6 and RGI7 have different geotiles... so we need to intersect the geotiles and then set the values
        dgeotile0 = dims(gt, :geotile)
        dgeotile0 = intersect(dgeotile0, dgeotile)
        
        geotiles0[run=At(run), geotile=At(dgeotile0)] = gt[geotile=At(dgeotile0)]
        geotiles0[run=At(run)] = geotiles0[run=At(run)] .- geotiles_ref
    end

    # calculate the quantile of the absolute deviation relative to the reference run
    geotiles_ref_err = copy(geotiles_ref)
    (dvarname, dgeotile, dTi) = dims(geotiles_ref_err)

    Threads.@threads for varname in dvarname
        for geotile in dgeotile
            for Ti in dTi
                geotiles_ref_err[varname=At(varname), geotile=At(geotile), Ti=At(Ti)] = quantile(abs.(geotiles0[varname=At(varname), geotile=At(geotile), Ti=At(Ti)]), error_quantile) * error_scaling
            end
        end
    end

    geotiles_out = cat(geotiles_ref, geotiles_ref_err; dims=Dim{:error}([false, true]))

    return geotiles_out
end


"""
    geotiles_mean_error_glaciers(glaciers, geotiles; show_progress=true, vars2downscale=nothing)

Downscale geotile-level data to individual glaciers, preserving reference values and error estimates.

# Arguments
- `glaciers`: DataFrame containing glacier data with columns for rgiid, geotile, and area_km2
- `geotiles`: DimensionalArray with dimensions for variable name, geotile, time, and error
- `show_progress`: Whether to display a progress bar during computation (default: true)
- `vars2downscale`: Optional subset of variables to downscale (default: all variables in geotiles)

# Returns
- A DimensionalArray with dimensions for variable name, glacier RGI ID, time, and error,
  where values are scaled by glacier area and the error dimension contains reference values (false)
  and error estimates (true)

# Examples
```julia
julia> glacier_out = geotiles_mean_error_glaciers(glaciers, geotiles; show_progress=true)
julia> runoff_glaciers = glacier_out[At("runoff"), :, :, At(false)]
```
"""
function geotiles_mean_error_glaciers(glaciers, geotiles; show_progress=true, vars2downscale=nothing)

    drgiid = Dim{:rgiid}(glaciers.rgiid)
    (dvarname, dgeotile, dTi, derror) = dims(geotiles)
    
    if !isnothing(vars2downscale)
        dvarname = Dim{:varname}(vars2downscale)
    end

    glacier_out = zeros(dvarname, drgiid, dTi, derror) * 1u"Gt" 

    prog = Progress(nrow(glaciers)*length(dvarname); desc="Downscaling from geotiles to glaciers", enabled=show_progress)

    for varname in dvarname
        Threads.@threads for r in eachrow(glaciers)
            if r.geotile in dgeotile
                glacier_out[varname=At(varname), rgiid=At(r.rgiid)] = geotiles[varname=At(varname), geotile=At(r.geotile)] .* r.area_km2
            end
            next!(prog)
        end
    end
    finish!(prog)

    return glacier_out
end


"""
    glacier_summary_file(binned_synthesized_dv_files, binned_synthesized_dv_file_ref; error_quantile=0.95, error_scaling=1.0, reference_period=(DateTime(2000, 4, 1), DateTime(2024, 12, 31)), surface_mask="glacier", force_remake_before=nothing, geotile_width=2)

Build or update the glacier-level summary NetCDF from geotile runs and save to disk.

Computes reference and error estimates from multiple runs, downscales from geotiles to glaciers,
and writes a NetCDF file with time series and error variables per glacier.

# Arguments
- `binned_synthesized_dv_files`: Vector of paths to binned synthesized run files
- `binned_synthesized_dv_file_ref`: Path to the reference run file
- `error_quantile`: Quantile for error estimation (default: 0.95)
- `error_scaling`: Scaling factor for error (default: 1.0)
- `reference_period`: Tuple of (start, end) DateTime for reference (default: (2000-04-01, 2024-12-31))
- `surface_mask`: Surface mask type (default: "glacier")
- `force_remake_before`: If set, skip if output file is newer than this DateTime
- `geotile_width`: Geotile width in degrees (default: 2)

# Returns
- Nothing; writes the glacier summary file to `pathlocal[:glacier_summary]`.

# Examples
```julia
julia> glacier_summary_file(binned_synthesized_dv_files, binned_synthesized_dv_file_ref; reference_period=(DateTime(2000, 4, 1), DateTime(2024, 12, 31)))
```
"""
function glacier_summary_file(
    binned_synthesized_dv_files, 
    binned_synthesized_dv_file_ref; 
    error_quantile=0.95, 
    error_scaling=1.0, 
    reference_period=(DateTime(2000, 4, 1), DateTime(2024, 12, 31)), 
    surface_mask="glacier", 
    force_remake_before=nothing, 
    geotile_width=2,
    )

    # having strange issues with custom unitful units... 

    glacier_summary_file = pathlocal[:glacier_summary]

    if file_is_current(glacier_summary_file, force_remake_before)
        printstyled("    -> Skipping $(glacier_summary_file) because it was created after force_remake_before:$force_remake_before\n"; color=:light_green)
        return
    end
    
    path2ref = binned_synthesized_dv_file_ref
    path2files = setdiff(binned_synthesized_dv_files, [binned_synthesized_dv_file_ref])

    if length(path2files) == length(binned_synthesized_dv_files)
        error("Reference run not found in binned_synthesized_dv_files")
    end

    # return a DimensionalArray with the mean and error for all runs relative to the reference run, for all modeled variables 
    geotiles0 = geotiles_mean_error(path2ref, path2files; error_quantile, error_scaling, reference_period) # returns variables in units of kg/m2 [90s]

    #TODO: it would be good to save this to disk and a gpkg file of trends and amplitudes currently in synthesis_plots_gis.jl... synthesis_plots_gis.jl should be deprecated at some point
    glacier_geotile_hyps_fn = replace(pathlocal[Symbol("$(surface_mask)_individual")], ".gpkg" => "geotile_hyps.jld2")
    glaciers0 = load(glacier_geotile_hyps_fn, surface_mask)
    glaciers0.area_km2[ismissing.(glaciers0.area_km2)] .= Ref([0.0,])
    
    # load glaciers
    glaciers0.area_km2 = sum.(glaciers0.area_km2)

    # glaciers have been split by geotile by model variables have been averaged by geotile groups so this is no longer needed... combine split glaciers
    glaciers = DataFrame()
    rgiids = unique(glaciers0[!, "RGIId"])
    glaciers[!, :rgiid] = rgiids
    glaciers[!, :area_km2] .= 0.0 * u"km^2"
    glaciers[!, :geotile] .= ""
    glaciers[!, :rgi] .= 0.0
    glaciers[!, :latitude] .= 0.0
    glaciers[!, :longitude] .= 0.0

    geotiles = _geotile_load_align(; surface_mask, geotile_width)
    geotiles[!, :rgi] .= 0

    for rgi in 1:19
        geotiles[geotiles[:, "rgi$rgi"].>0, :rgi] .= rgi
    end

    Threads.@threads for r in eachrow(glaciers)
        index = findall(isequal.(Ref(r.rgiid), glaciers0[!, "RGIId"]))
        _, ind = findmax(glaciers0.area_km2[index])
        ind = index[ind]
        r.area_km2 = sum(glaciers0.area_km2[index]) * u"km^2"
        r.geotile = glaciers0[ind, :geotile]
        r.latitude = glaciers0[ind, :CenLat]
        r.longitude = glaciers0[ind, :CenLon]
        r.rgi = geotiles[findfirst(isequal(r.geotile), geotiles.id), :rgi]
    end

    glacier_out = geotiles_mean_error_glaciers(glaciers, geotiles0) #[2 min]

    # From DimArray to DimStack for saveing
    vars = glacier_out[error=At(false)]
    var_error = glacier_out[error=At(true)]
    dvarerror = Dim{:varname}(parent(val(dims(vars, :varname))) .* "_error")
    var_error = set(var_error, :varname => dvarerror)
    da = cat(vars, var_error, dims=(:varname))
    dstack = DimStack(glacier_out[error=At(false)]; layersfrom=:varname)

    # save to netcdf
    NCDataset(glacier_summary_file, "c") do ds

        data_dims = dims(dstack)

        # First all dimensions need to be defined
        for dim in data_dims
            dname = string(DimensionalData.name(dim))
            defDim(ds, dname, length(dim))
        end

        # now add the variables
        for dim in data_dims
            dname = string(DimensionalData.name(dim))
         
            if isa(dim[1], Date)
                #Type Date not yet supported
                 d = defVar(ds, dname, DateTime.(val(val(dim))), (dname,))
            else
                d = defVar(ds, dname, val(val(dim)), (dname,))
            end
                
            if DateTime <: eltype(dim)
                d.attrib["cf_role"] = "timeseries_id"
            end
        end

        for vaname in keys(dstack)
            v = defVar(ds, "$vaname", ustrip.(parent(dstack[vaname])), string.(DimensionalData.name.(data_dims)))
            v.attrib["units"] = string(Unitful.unit(dstack[vaname][1]))
        end

        v = defVar(ds, "glacier_area", ustrip.(glaciers.area_km2), string.(DimensionalData.name.(data_dims[1:1])))
        v.attrib["units"] = "km^2"

        # add latitude and longitude [This is a hack until I can get geometry into a NetCDF]
        defVar(ds, "latitude", glaciers.latitude, string.(DimensionalData.name.(data_dims[1:1])))
        defVar(ds, "longitude", glaciers.longitude, string.(DimensionalData.name.(data_dims[1:1])))

        # add global attributes
        ds.attrib["title"] = "cumulative glacier mass flux outputs from Gardner et al. 2025"
        ds.attrib["version"] = "beta - " * Dates.format(now(), "yyyy-mm-dd")
        ds.attrib["featureType"] = "timeSeries"
        # The glacier inventory the geotile values were spread over; `glacier_routing` checks it.
        ds.attrib["surface_mask"] = String(surface_mask)

        println("glacier output saved to $glacier_summary_file")
    end
end


"""
    gembfit_dv2gpkg(binned_synthesized_dv_file; outfile_prefix="geotiles_rates", datelimits=(DateTime(2000, 3, 1), DateTime(2025, 1, 1)))

Export geotile-level glacier change trends and amplitudes to GIS-compatible files.

This function processes a GEMB fit output file containing geotile-level glacier change data,
fits temporal trends to all relevant variables, and exports the results to Arrow and GeoPackage (GPKG) formats.
It produces both total rates (in km³/yr) and area-averaged rates (in m/yr), suitable for GIS analysis.

# Arguments
- `binned_synthesized_dv_file::AbstractString`: Path to the input file containing geotile-level glacier change data.
- `outfile_prefix::AbstractString`: Prefix for output files (default: "geotiles_rates").
- `datelimits::Tuple{DateTime, DateTime}`: Date range for fitting temporal trends (default: (2000-03-01, 2025-01-01)).

# Outputs
- Arrow file with total rates in km³/yr.
- GeoPackage files with total rates (km³/yr) and area-averaged rates (m/yr), excluding raw time series columns.

# Examples
```julia
julia> gembfit_dv2gpkg(binned_synthesized_dv_file; outfile_prefix="Gardner2025_geotiles_rates", datelimits=(DateTime(2000, 3, 1), DateTime(2025, 1, 1)))
julia> # Writes Arrow and GeoPackage files to the project output directory
```

# Notes
- The function fits trends to all time series variables, then computes area-averaged rates by dividing by geotile area.
- Variable names containing "dv" are converted to "dh" for area-averaged rates.
- Only non-time-series columns and derived trend/amplitude columns are exported to GIS files.
"""
function gembfit_dv2gpkg(binned_synthesized_dv_file; outfile_prefix="Gardner2025_geotiles_rates", datelimits=(DateTime(2000, 3, 1), DateTime(2025, 1, 1)))
    # export trends and amplitudes for reference_run

    # data is in units of km3 of ice equivalent at `δice`
    geotiles0 = FileIO.load(binned_synthesized_dv_file, "geotiles")

    altim_cols = names(geotiles0)[occursin.("_altim", names(geotiles0))]
    rename!(geotiles0, altim_cols .=> replace.(altim_cols, "_altim" => "_synthesis"))

    # snow is being scaled in GEMB runs to account for avalanching and wind redistibutions, rain is also being scaled by the same factor making it unrealistic.
    rain_cols = occursin.("rain", names(geotiles0))
    geotiles0 = geotiles0[:, .!rain_cols];

    vars_no_write = setdiff(names(geotiles0), ["id", "geometry", "group", "pscale", "ΔT", "area_km2", "rgiid"])
    geotiles0[!, :area_km2] = sum.(geotiles0[:, :area_km2])

    # Fit temporal trends to all variables with dates in metadata
    vars_ts = String[]
    for col in names(geotiles0)
        meta = colmetadata(geotiles0, col)
        if haskey(meta, "date")
            push!(vars_ts, col)
        end
    end

    source_crs1 = GFT.EPSG(4326)
    geotiles0 = df_tsfit!(geotiles0, vars_ts; datelimits)

    colname = names(geotiles0)
    coloffset = colname[occursin.(Ref("_offset"), names(geotiles0))]
    vars_no_write = vcat(vars_no_write, coloffset)
    vars_no_write = vcat(vars_no_write, ["discharge_amplitude", "discharge_phase"])
    geotiles0 = geotiles0[:, Not(vcat(vars_ts, vars_no_write))]

    # Export results, excluding raw time series
    isvec = []
    for i = 1:ncol(geotiles0)
        push!(isvec, typeof(geotiles0[1, 12]) >: Vector)
    end

    outfile_km3yr = joinpath(pathlocal[:data_dir], "project_data", "$(outfile_prefix)_km3yr.gpkg")
    GeoDataFrames.write(outfile_km3yr, geotiles0; crs=source_crs1)

    # now do the same but for area averaged rates of change [from volume to meters]
    
    for varname in names(geotiles0)
        if any(occursin.(vars_ts, Ref(varname)))

            if occursin("_phase", varname)
                rename!(geotiles0, varname => replace(varname, "dv" => "dh"))
                continue
            end

            geotiles0[:, varname] ./= (geotiles0[:, :area_km2] / 1e3)
            if occursin("dv", varname)
                rename!(geotiles0, varname => replace(varname, "dv" => "dh"))
            end
        end
    end

    outfile_myr = joinpath(pathlocal[:data_dir], "project_data", "$(outfile_prefix)_myr.gpkg")
    GeoDataFrames.write(outfile_myr, geotiles0; crs=source_crs1)

    return outfile_km3yr, outfile_myr
end

"""
    ensemble_summary(path2runs_synthesized, ensemble_reference_file; dates4trend=[DateTime(2000, 3, 1), DateTime(2024, 12, 15)], error_quantile=0.95, error_scaling=1.5, discharge_fractional_error=0.15)

Compute regional fit parameters with reference values and error estimates for an ensemble of runs.

Loads discharge by RGI region, aggregates runs by RGI, fits trends over the given date range,
then computes reference values and quantile-based errors relative to the ensemble reference file.

# Arguments
- `path2runs_synthesized`: Vector of paths to synthesized run files
- `ensemble_reference_file`: Path to the reference run file
- `dates4trend`: Date range for trend fitting (default: [2000-03-01, 2024-12-15])
- `error_quantile`: Quantile for error estimation (default: 0.95)
- `error_scaling`: Scaling factor for error (default: 1.5)
- `discharge_fractional_error`: Fractional error for discharge (default: 0.15)

# Returns
- DimensionalArray with dimensions for variable name, RGI region, parameter, and error
  (reference values and error estimates)

# Examples
```julia
julia> dates4trend = [DateTime(2000, 3, 1), DateTime(2024, 12, 15)]
julia> region_fits = ensemble_summary(path2runs_synthesized, ensemble_reference_file; dates4trend, error_quantile=0.95)
julia> runoff_trend_ref = region_fits[At("runoff"), :, At("trend"), At(false)]
```
"""
function ensemble_summary(path2runs_synthesized, ensemble_reference_file; physics_member=gemb_reference_physics_member, dates4trend=[DateTime(2000, 3, 1), DateTime(2024, 12, 15)], error_quantile=0.95, error_scaling=1.5, discharge_fractional_error=0.15)

    paths = pathlocal

    # load discharge for each RGI [<1s]
    path2discharge = discharge_global_file(binned_filled_fileparts(ensemble_reference_file).surface_mask)
    path2rgi_regions = paths[:rgi6_regions_shp]

    discharge_rgi0 = discharge_rgi(path2discharge, path2rgi_regions; fractional_error=discharge_fractional_error)
    
    # load results for all runs for each RGI [18s]
    runs_rgi = runs2rgi(path2runs_synthesized; physics_member)
    runs_rgi_fits = rgi_trends(runs_rgi, discharge_rgi0, dates4trend)
    region_fits = region_fit_ref_and_err(runs_rgi_fits, ensemble_reference_file; error_quantile, error_scaling, discharge=discharge_rgi0)

    return region_fits
end

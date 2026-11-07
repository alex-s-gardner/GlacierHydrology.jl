#density of glacier ice
# Ice density [kg m-3]. Matches the `density_ice` the GEMB runs use, which the tile files record as
# their `model_density_ice` attribute and divide by to build the height-change decomposition. The two
# have to agree: a GEMB volume converted to mass with a different density carries that ratio as a bias.
const δice = 917;
const seasonality_weight = 85/100
# Only used by the `:legacy` penalty mode.
const ΔT_to_pscale_weight = 50/100
const ocean_area_km2 = 362.5 * 1E6

"""
    reference_surface_mask() -> String

Glacier inventory (`"glacier"` for RGI 6.0, `"glacier_rgi7"` for RGI 7.0) of the reference ensemble
member, `GVC.reference_member()`. The glacier-level products spread that member's geotile values over
this inventory's glaciers.
"""
reference_surface_mask() = binned_filled_fileparts(GVC.reference_member()).surface_mask

"""
    ensemble_paths(; geotile_width=2, physics_member=gemb_reference_physics_member) -> NamedTuple

File paths of every altimetry ensemble member of the GlacierVolumeChange dataset
(`GVC.ensemble_members`) and of the reference member (`GVC.reference_member`), in each form the
pipeline reads or writes:

- `filled`, `filled_ref`: the members' filled (`*_aligned.jld2`) files, which name the members;
- `synthesized`, `synthesized_ref`: their synthesized files, the input to the GEMB calibration;
- `gembfit_dv`, `gembfit_dv_ref`: their calibrated syntheses for `physics_member`
  ([`gemb_physics_path`](@ref));
- `surface_masks`: the glacier inventories among the members; `surface_mask_ref`: the reference
  member's ([`reference_surface_mask`](@ref)).
"""
function ensemble_paths(; geotile_width=2, physics_member=gemb_reference_physics_member)
    filled = GVC.ensemble_members(GVC.VolumeDatasetConfig(; geotile_width))
    filled_ref = GVC.reference_member()
    filled_ref in filled || error("the reference member $(filled_ref) is not among the ensemble members")
    gembfit_file(f) = gembfit_dv_path(GVC.synthesized_path(f), physics_member)
    return (;
        filled, filled_ref,
        synthesized=GVC.synthesized_path.(filled), synthesized_ref=GVC.synthesized_path(filled_ref),
        gembfit_dv=gembfit_file.(filled), gembfit_dv_ref=gembfit_file(filled_ref),
        surface_masks=unique(getproperty.(binned_filled_fileparts.(filled), :surface_mask)),
        surface_mask_ref=reference_surface_mask(),
    )
end
# The GEMB physics member (GEMB_GlacierSims `physics_member`) the reported fluxes are calibrated
# against: GEMB_GlacierSims' default, `S-kfit-table`. That is GEMB's default snow set (Flanner & Zender
# 2006 grain growth, SNOTEL-calibrated van Kampenhout et al. 2017 fresh snow, Crocus densification with
# destructive metamorphism at 100 kg m-3), the fitted band decoupling factor, the regional background
# impurity table, ERA5-Land wind and no precipitation gradient. Other members are sensitivities to it.
const gemb_reference_physics_member = "S-kfit-table"

# The forcing grid GEMB run 10 (the pilot sweep) is calibrated on: 7 temperature offsets [K] by 9
# precipitation scalings. Pilot tiles run on the dense grid carry these nodes among others and are read
# at these alone, so every pilot tile enters the ensemble on the same nodes.
const gemb_pilot_forcing_grid = (; pscale=[0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0],
                                   ΔT=[-3.0, -1.5, 0.0, 1.5, 3.0, 4.5, 6.0])

# The geotile-mean pscale nodes a precipitation-gradient member of run 10 is calibrated on. Each band is
# read at pscale · m_b, which must stay on the tiles' 0.25 … 6 axis; at 2.5 % per 100 m the pilot tiles'
# bands do so for geotile-mean pscale 0.75 … 3 (the tallest, N28_E086 and N60_W140, bound it).
const gemb_pilot_gradient_forcing_grid = (; pscale=[0.75, 1.0, 1.5, 2.0, 3.0],
                                            ΔT=gemb_pilot_forcing_grid.ΔT)

# The GEMB runs: where each sweep's tiles are, whether it covers every glacierized tile, and the
# calibration prior derived from its own unpenalized fits.
#
# - `sweep`, `suffix`: the tile directory under `$CLIMATE_CACHE` and the tag of its combined files.
# - `partial_coverage`: the sweep covers a subset of tiles. Such a run is binned on the geotiles with a
#   tile of their own, calibrated on the geotile groups wholly inside it, and skips the fitted-pscale
#   check of the reference geotiles.
# - `allows_gradient`: precipitation-gradient members are defined on its tiles.
# - `forcing_prior`, `distance_from_origin_penalty`: the forcing prior `(; pscale, log_pscale_sd, ΔT,
#   ΔT_sd, corr)` and its weight `wd`, re-derived for each run by the temporal cross-validation of
#   notes/methods_2026-10_seasonal_cycle_and_calibration.md §4 from that run's unpenalized fits.
#   `nothing` and 0 until then: the calibration is unpenalized.
const GEMB_RUNS = Dict(
    9 => (; sweep="tile_runs_v2", suffix="v2", partial_coverage=false, allows_gradient=false,
            forcing_prior=nothing, distance_from_origin_penalty=0.0),
    10 => (; sweep="tile_runs_pilot", suffix="pilot", partial_coverage=true, allows_gradient=true,
             forcing_prior=nothing, distance_from_origin_penalty=0.0),
)

"""
    gemb_info(; gemb_run_id=9, physics_member)

Configuration of a GEMB_GlacierSims tile sweep.

# Arguments
- `gemb_run_id::Int`: 9, the production sweep in `\$CLIMATE_CACHE/tile_runs_v2`, or 10, the pilot sweep
  in `\$CLIMATE_CACHE/tile_runs_pilot` on its own forcing grid and tile subset ([`GEMB_RUNS`](@ref)).
- `physics_member`: the GEMB_GlacierSims physics member (e.g. `"S-kfit-table"`, see
  [`gemb_reference_physics_member`](@ref)). Each run holds one tile sweep per member, under
  `<tile_dir>/<tile member>`, with a combined file of its own, so no tile of another member or sweep can
  enter the ensemble. A run-10 member with a precipitation gradient (`"<tile member>-g<g>"`,
  [`gemb_physics_member_parts`](@ref)) reads its tile member's directory, has a combined file of its own
  and the forcing grid [`gemb_pilot_gradient_forcing_grid`](@ref); run 9 refuses one.

# Returns
Named tuple with
- `gemb_run_id`, `physics_member`: the run and the checked member id;
- `tile_dir`: directory of the member's tile NetCDF files;
- `filename_gemb_combined`: path from which the per-mask geotile ensemble files are named;
- `forcing_grid`: the `(; pscale, ΔT)` nodes the ensemble is read on, or `nothing` for the tiles' own grid;
- `partial_coverage`, `forcing_prior`, `distance_from_origin_penalty`: the run's settings
  ([`GEMB_RUNS`](@ref)).

# Examples
```julia
julia> gembinfo = gemb_info(; gemb_run_id=9, physics_member=gemb_reference_physics_member)
```
"""
function gemb_info(; gemb_run_id=9, physics_member)
    haskey(GEMB_RUNS, gemb_run_id) ||
        error("unrecognized gemb_run_id: $gemb_run_id (runs $(join(sort(collect(keys(GEMB_RUNS))), " and ")) are defined)")
    run = GEMB_RUNS[gemb_run_id]
    member = check_gemb_physics_member(physics_member)
    parts = gemb_physics_member_parts(member)
    gradient = parts.precipitation_gradient > 0
    gradient && !run.allows_gradient &&
        error("precipitation-gradient members ($(repr(member))) are defined for the pilot sweep (GEMB " *
              "run 10) only; run $(gemb_run_id) needs a geotile-mean pscale axis its tiles' bands stay on first")
    return (;
        gemb_run_id,
        physics_member=member,
        tile_dir=joinpath(get(ENV, "CLIMATE_CACHE", "/mnt/bylot-r3/data/era5land"), run.sweep, parts.tile_member),
        filename_gemb_combined="/mnt/bylot-r3/data/gemb/gemb_glacier_sims_tiles_1950to2026_$(run.suffix)_$(member).jld2",
        forcing_grid=!run.partial_coverage ? nothing :
                     gradient ? gemb_pilot_gradient_forcing_grid : gemb_pilot_forcing_grid,
        run.partial_coverage, run.forcing_prior, run.distance_from_origin_penalty,
    )
end


"""
    geotile_groups_forced()

Return manual overrides for geotile groupings where large glaciers span multiple tiles.

Used when glaciers that cross tile boundaries should be grouped or treated separately.
Groupings are only applied when geotile grouping is recomputed (e.g. force_remake_before is set).

# Returns
- Vector of vectors of geotile ID strings; each inner vector is one forced group.

# Examples
```julia
julia> forced = geotile_groups_forced()
julia> geotile_grouping!(geotiles0, glaciers, 100; geotile_groups_manual=forced)
```
"""
function geotile_groups_forced() 
    out = [
        ["lat[+62+64]lon[-148-146]"],
        ["lat[+62+64]lon[-146-144]"],
        ["lat[+60+62]lon[-138-136]"],
        ["lat[+58+60]lon[-138-136]"],
        ["lat[+58+60]lon[-136-134]"],
        ["lat[+58+60]lon[-134-132]"],
        ["lat[+56+58]lon[-134-132]"],
        ["lat[+78+80]lon[-084-082]"],
        ["lat[-52-50]lon[-074-072]"],
        ["lat[+56+58]lon[-130-128]"],
        ["lat[-74-72]lon[-080-078]", "lat[-74-72]lon[-078-076]"],
        ["lat[+78+80]lon[+010+012]", "lat[+78+80]lon[+012+014]", "lat[+78+80]lon[+014+016]"],
        ["lat[+76+78]lon[+062+064]", "lat[+76+78]lon[+064+066]", "lat[+76+78]lon[+066+068]", "lat[+76+78]lon[+068+070]", "lat[+74+76]lon[+062+064]", "lat[+74+76]lon[+064+066]", "lat[+74+76]lon[+066+068]", "lat[+74+76]lon[+066+068]"]
    ]
    return out
end

"""
    gemb_altim_cost(x, dv_altim, dv_gemb, kwargs)

Compute cost between altimetry-derived volume change and GEMB-sampled volume change.

Samples GEMB at (pscale, ΔT) from x, subtracts from dv_altim, then evaluates
model_fit_cost_function on the residuals with the given kwargs.

# Arguments
- `x`: Parameter vector [pscale, ΔT]
- `dv_altim`: DimArray of altimetry-derived volume change (geotile, date)
- `dv_gemb`: GEMB volume change array used for sampling
- `kwargs`: Keyword arguments passed to model_fit_cost_function (e.g. seasonality_weight, distance_from_origin_penalty)

# Returns
- Scalar cost value.

# Examples
```julia
julia> cost = gemb_altim_cost([1.0, 0.0], dv_altim, dv_gemb, (; seasonality_weight=0.85, distance_from_origin_penalty=0.2, ΔT_to_pscale_weight=0.5))
```
"""
function gemb_altim_cost(x, dv_altim, dv_gemb, kwargs)
    pscale = x[1]
    ΔT = x[2]

    res = dv_altim .- gemb_dv_sample(pscale, ΔT, dv_gemb)
    res .-= mean(res)

    cost = model_fit_cost_function(res, pscale, ΔT; kwargs...)

    return cost
end

"""
    model_fit_cost_function(res, pscale, ΔT; seasonality_weight, distance_from_origin_penalty, ΔT_to_pscale_weight, calibrate_to=:all, is_scaling_factor, origin_penalty_mode=:prior, forcing_prior=nothing)

Compute a composite cost function for fitting a model to altimetry data, incorporating trend, seasonality, and parameter penalties.

# Arguments
- `res`: Residuals between observed and modeled values. Should be a DimArray or array-like object with a :date dimension.
- `pscale`: Precipitation scaling factor (numeric).
- `ΔT`: Air temperature offset applied to the GEMB forcing (numeric, in K). An additive parameter whose
  no-op is 0, unlike `pscale`; `is_scaling_factor` selects which distance the `:legacy` penalty uses.
- `seasonality_weight`: Weight (0–1) for the seasonal amplitude in the cost function. Higher values emphasize seasonality.
- `distance_from_origin_penalty`: Strength `wd` of the penalty on the forcing corrections.
- `ΔT_to_pscale_weight`: Weight of the temperature-offset penalty against the precipitation-scaling one
  (`:legacy` only).
- `calibrate_to`` = [:all, :annual, :five_year, :trend]
- `origin_penalty_mode`: `:prior` (default) multiplies the fit cost by `1 + wd * d`, `d` the Mahalanobis
  distance of (log pscale, ΔT) from `forcing_prior`, i.e. the distance in prior standard deviations
  allowing for the correlation between the two. `:legacy` multiplies it by `1 + wd * d`, `d` the weighted
  distance from the no-op point (1, 0).
- `design`: a [`seasonal_design`](@ref) of `res`'s dates, which a caller evaluating many residuals on
  one date axis builds once; `nothing` builds it here.
- `forcing_prior`: `(; pscale, log_pscale_sd, ΔT, ΔT_sd, corr)`, the centre, spreads and correlation used
  by `origin_penalty_mode = :prior`, or `nothing`, which only a zero `distance_from_origin_penalty`
  accepts.

# Returns
A tuple `cost` where:
- `cost`: Composite cost metric, combining RMSE, seasonal amplitude, and penalties for parameter deviation.

# Details
- The function fits a seasonal model to the residuals using `ts_seasonal_model`.
- If `calibrate_to_annual_change_only` is `true`, the cost is computed using only the annual change (residuals at the seasonal minimum).
- If `calibrate_to_trend_only` is `true`, the cost is based only on the absolute value of the linear trend.
- Otherwise, the cost is a weighted sum of RMSE and the amplitude of the seasonal cycle.
- The penalty measures deviation of (`pscale`, `ΔT`) from `forcing_prior`, or from the no-op point (1, 0)
  when `origin_penalty_mode = :legacy`.

# Examples
```julia
julia> cost = model_fit_cost_function(res, 1.0, 0.0; seasonality_weight=0.85, distance_from_origin_penalty=0.35, ΔT_to_pscale_weight=0.5)
```
"""
function model_fit_cost_function(res, pscale, ΔT; seasonality_weight, distance_from_origin_penalty, ΔT_to_pscale_weight, calibrate_to=:all, is_scaling_factor=Dict("pscale" => true, "ΔT" => false), origin_penalty_mode=:prior, forcing_prior=nothing, design=nothing)

    origin_penalty_mode in (:legacy, :prior) || error("origin_penalty_mode must be :legacy or :prior, got $origin_penalty_mode")

    # remove linear trend to emphasize seasonality
    if calibrate_to != :five_year
        fit = ts_seasonal_model(res; design);
    else
        res0 = groupby(res, :date => Bins(year, 5))
        res0 = mean.(res0)
        rmse_cost = sqrt(mean(res0 .^ 2))
    end

    # calibrate to annual change only
    if calibrate_to == :annual
        seasonal_min =  mod1(fit.phase_peak_month + 6, 12)
        res = res[date=Near(DateTime(minimum(year.(dims(res, :date))), seasonal_min, 15):Year(1):DateTime(maximum(year.(dims(res, :date))), seasonal_min, 15))]
        rmse_cost = sqrt(mean(res .^ 2))
    end

    if calibrate_to == :all
        rmse = sqrt(mean(res .^ 2))
        cost = ((1 - seasonality_weight) * rmse) + (seasonality_weight * fit.amplitude)
    elseif calibrate_to == :trend
        cost = (1 - seasonality_weight) * abs(fit.trend)
    else
        cost = (1 - seasonality_weight) * rmse_cost
    end

    if origin_penalty_mode == :prior
        iszero(distance_from_origin_penalty) && return cost
        (is_scaling_factor["pscale"] && !is_scaling_factor["ΔT"]) ||
            error("origin_penalty_mode = :prior assumes a multiplicative pscale and an additive ΔT")
        isnothing(forcing_prior) && error("a forcing-prior penalty needs a forcing_prior")
        return cost * (1 + distance_from_origin_penalty * _forcing_prior_distance(pscale, ΔT, forcing_prior))
    else
        # Distance of a scaling factor from its no-op value of 1, measured symmetrically so that
        # halving and doubling are penalized equally.
        dp = _distance_from_origin(pscale, is_scaling_factor["pscale"])
        dT = _distance_from_origin(ΔT, is_scaling_factor["ΔT"])
        return cost * (1 + sqrt((dT * ΔT_to_pscale_weight)^2 + (dp * (1 - ΔT_to_pscale_weight))^2) * distance_from_origin_penalty)
    end
end

# Mahalanobis distance of (log pscale, ΔT) from the prior centre, in prior standard deviations.
function _forcing_prior_distance(pscale, ΔT, prior)
    z1 = (log(pscale) - log(prior.pscale)) / prior.log_pscale_sd
    z2 = (ΔT - prior.ΔT) / prior.ΔT_sd
    ρ = prior.corr
    return sqrt(max((z1^2 - 2ρ * z1 * z2 + z2^2) / (1 - ρ^2), 0.0))
end

function _distance_from_origin(scale, is_scaling_factor)
    is_scaling_factor || return scale
    return scale < 1 ? 1 / scale - 1 : scale - 1
end

#density of glacier ice
# Ice density [kg m-3]. Matches the `density_ice` the GEMB runs use, which the tile files record as
# their `model_density_ice` attribute and divide by to build the height-change decomposition. The two
# have to agree: a GEMB volume converted to mass with a different density carries that ratio as a bias.
const δice = 917;
const seasonality_weight = 85/100
# Strength of the forcing-prior penalty (`origin_penalty_mode = :prior`): the cost is multiplied by
# 1 + wd × d, d the Mahalanobis distance from `gemb_forcing_prior`. 0.35 is the largest value within 1% of the
# best held-out skill in temporal cross-validation; see notes/methods_2026-10_seasonal_cycle_and_calibration.md.
# Derived, like `gemb_forcing_prior`, from the GEMB run `gemb_forcing_prior_run_id` fits.
const distance_from_origin_penalty = 35 / 100
# Only used by the `:legacy` penalty mode.
const ΔT_to_pscale_weight = 50/100

# Empirical prior on the GEMB forcing corrections for `origin_penalty_mode = :prior`: the area-weighted
# centre, spread and correlation of the unpenalized fits of the well-constrained geotile groups (sd of
# log pscale < 0.15 and of ΔT < 0.5 K within 5% of the minimum cost; 88% of ice area) for the reference
# ensemble member, fill set 6, GEMB run 8. See notes/methods_2026-10_seasonal_cycle_and_calibration.md.
const gemb_forcing_prior = (pscale=1.60, log_pscale_sd=0.45, ΔT=2.20, ΔT_sd=1.65, corr=0.50)
# The GEMB run `gemb_forcing_prior` and `distance_from_origin_penalty` were derived from. Both describe
# that run's forcing bias, so they hold for no other run: `gemb_calibration` throws when it would apply
# them to a different one. Re-derive both from the new run's unpenalized fits
# (`distance_from_origin_penalty = 0`), then update all three together.
const gemb_forcing_prior_run_id = 8
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
    synthesized_file(f) = replace(f, "_aligned.jld2" => "_synthesized.jld2")
    gembfit_file(f) = gemb_physics_path(replace(synthesized_file(f), ".jld2" => "_gembfit_dv.jld2"), physics_member)
    return (;
        filled, filled_ref,
        synthesized=synthesized_file.(filled), synthesized_ref=synthesized_file(filled_ref),
        gembfit_dv=gembfit_file.(filled), gembfit_dv_ref=gembfit_file(filled_ref),
        surface_masks=unique(getproperty.(binned_filled_fileparts.(filled), :surface_mask)),
        surface_mask_ref=binned_filled_fileparts(filled_ref).surface_mask,
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

"""
    gemb_info(; gemb_run_id=9, physics_member)

Configuration of a GEMB_GlacierSims tile sweep.

# Arguments
- `gemb_run_id::Int`: 9, the production sweep in `$CLIMATE_CACHE/tile_runs_v2`, or 10, the pilot sweep
  in `$CLIMATE_CACHE/tile_runs_pilot` on its own forcing grid and tile subset.
- `physics_member`: the GEMB_GlacierSims physics member (e.g. `"S-kfit-table"`, see
  [`gemb_reference_physics_member`](@ref)). Each run holds one tile sweep per member, under
  `<tile_dir>/<tile member>`, with a combined file of its own, so no tile of another member or sweep can
  enter the ensemble. A run-10 member with a precipitation gradient (`"<tile member>-g<g>"`,
  [`gemb_physics_member_parts`](@ref)) reads its tile member's directory, has a combined file of its own
  and the forcing grid [`gemb_pilot_gradient_forcing_grid`](@ref); run 9 refuses one.

# Returns
Named tuple with
- `tile_dir`: directory of the member's tile NetCDF files;
- `filename_gemb_combined`: path from which the per-mask geotile ensemble files are named;
- `physics_member`: the checked member id;
- `forcing_grid`: the `(; pscale, ΔT)` nodes the ensemble is read on, or `nothing` for the tiles' own grid.

# Examples
```julia
julia> gembinfo = gemb_info(; gemb_run_id=9, physics_member=gemb_reference_physics_member)
```
"""
function gemb_info(; gemb_run_id = 9, physics_member = nothing)
    gemb_run_id in (9, 10) || error("unrecognized gemb_run_id: $gemb_run_id (runs 9 and 10 are defined)")
    isnothing(physics_member) &&
        error("GEMB run $(gemb_run_id) holds one tile sweep per physics member; pass physics_member " *
              "(the reference is gemb_reference_physics_member = $(repr(gemb_reference_physics_member)))")
    member = check_gemb_physics_member(physics_member)
    parts = gemb_physics_member_parts(member)
    gemb_run_id == 9 && parts.precipitation_gradient > 0 &&
        error("precipitation-gradient members ($(repr(member))) are defined for the pilot sweep (GEMB " *
              "run 10) only; run 9 needs a geotile-mean pscale axis its tiles' bands stay on first")
    sweep, suffix = gemb_run_id == 9 ? ("tile_runs_v2", "v2") : ("tile_runs_pilot", "pilot")
    gemb_info = (;
        tile_dir=joinpath(get(ENV, "CLIMATE_CACHE", "/mnt/bylot-r3/data/era5land"), sweep, parts.tile_member),
        filename_gemb_combined="/mnt/bylot-r3/data/gemb/raw/gemb_glacier_sims_tiles_1950to2026_$(suffix)_$(member).jld2",
        physics_member=member,
        forcing_grid=gemb_run_id == 9 ? nothing :
                     parts.precipitation_gradient > 0 ? gemb_pilot_gradient_forcing_grid : gemb_pilot_forcing_grid
    )

    return gemb_info
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
    model_fit_cost_function(res, pscale, ΔT; seasonality_weight, distance_from_origin_penalty, ΔT_to_pscale_weight, calibrate_to=:all, is_scaling_factor, origin_penalty_mode=:prior, forcing_prior=gemb_forcing_prior)

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
- `forcing_prior`: `(; pscale, log_pscale_sd, ΔT, ΔT_sd, corr)`, the centre, spreads and correlation used
  by `origin_penalty_mode = :prior`; defaults to `gemb_forcing_prior`.

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
function model_fit_cost_function(res, pscale, ΔT; seasonality_weight, distance_from_origin_penalty, ΔT_to_pscale_weight, calibrate_to=:all, is_scaling_factor=Dict("pscale" => true, "ΔT" => false), origin_penalty_mode=:prior, forcing_prior=gemb_forcing_prior)

    origin_penalty_mode in (:legacy, :prior) || error("origin_penalty_mode must be :legacy or :prior, got $origin_penalty_mode")

    # remove linear trend to emphasize seasonality
    if calibrate_to != :five_year
        fit = ts_seasonal_model(res; interval=nothing);
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
        (is_scaling_factor["pscale"] && !is_scaling_factor["ΔT"]) ||
            error("origin_penalty_mode = :prior assumes a multiplicative pscale and an additive ΔT")
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

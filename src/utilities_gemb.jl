



"""
    gemb_bestfit_grouped(dv_altim, dv_gemb, geotiles0;
                         distance_from_origin_penalty=0.20,
                         ΔT_to_pscale_weight=1,
                         seasonality_weight=GlobalGlacierAnalysis.seasonality_weight,
                         calibrate_to=:all,
                         is_scaling_factor=Dict("pscale"=>true, "ΔT"=>false))

Calibrate the GEMB (Glacier Energy and Mass Balance) forcing parameters to altimetry-derived volume
changes for groups of geotiles. Glaciers often span multiple geotiles, so calibration is performed on
mutually exclusive geotile groups rather than individual tiles to ensure consistent parameter estimates
across connected glacierized regions.

Uses a staged, deterministic grid search (`_gemb_grid_minimize`) to find the `pscale` (precipitation
scaling) and `ΔT` (air temperature offset, K) that minimize the cost function between modeled GEMB
volume changes and altimetry observations. The cost function balances temporal agreement, seasonal cycle
matching, and distance from the empirical forcing prior `gemb_forcing_prior` (`origin_penalty_mode =
:prior`, the default) or from the no-op values 1 for `pscale` and 0 for `ΔT` (`:legacy`).

# Arguments
- `dv_altim`: DimArray of volume change from altimetry [km³ ice equivalent] with dimensions (:date, :geotile)
- `dv_gemb`: DimArray of GEMB modeled volume change [km³ ice equivalent] with dimensions (:date, :geotile, :pscale, :ΔT)
- `geotiles0`: DataFrame containing geotile metadata with columns :id, :group, :rgi, :extent, :area_km2
  - The :group column identifies mutually exclusive geotile groups for joint calibration

# Keyword Arguments
- `distance_from_origin_penalty=0.20`: Weight (0-1) for penalty on deviation from pscale=1, ΔT=0
- `ΔT_to_pscale_weight=1`: Weight (0-1) balancing ΔT vs pscale in optimization (1=equal weight)
- `seasonality_weight`: Weight (0-1) for matching seasonal cycle vs absolute volume change; defaults to
  the package constant `seasonality_weight` (0.85), which `gemb_calibration` also defaults to
- `calibrate_to=:all`: Which component to calibrate to (:all, :seasonal, :trend, etc.)
- `is_scaling_factor`: Dict indicating whether pscale/ΔT are multiplicative (true) or additive (false)

# Returns
- **Normal mode (multiple groups)**: DataFrame with columns [:id, :extent, :group, :pscale, :ΔT, :grid_bound, :rgi]
  - One row per geotile with optimized scaling factors (shared within each group)
  - `grid_bound` names the bounds of the pscale × ΔT search grid the fit lies on (`""` if none; see
    `_grid_bounds_hit`). One `@warn` reports the number of groups, geotiles and ice area so pinned.
- **Debug mode (single group)**: Tuple of (cost grid, geotile_ids) for examining optimization landscape

# Implementation Details
The function operates differently depending on the number of unique groups:
- **Multiple groups**: minimizes each group's cost with `_gemb_grid_minimize`, threaded across groups
- **Single group**: returns the full cost grid for diagnostic plots and prints the grid minimum and the
  value `_gemb_grid_minimize` returns

The search is deterministic, so results do not depend on thread count or order. Data aggregation is
pre-computed outside the parallel loop to avoid thread-safety issues with DimensionalData operations.

# Examples
```julia
# Typical usage for calibration
julia> geotiles_fit = gemb_bestfit_grouped(dv_altim, dv_gemb, geotiles0;
                                           seasonality_weight=0.85)
julia> pscale = geotiles_fit.pscale  # Extract precipitation scaling factors

# Debug mode for single geotile/group
julia> cost_grid, geotile_ids = gemb_bestfit_grouped(dv_altim, dv_gemb, geotiles_subset)
```

See also: [`gemb_altim_cost`](@ref), [`gemb_dv_sample`](@ref)
"""
function gemb_bestfit_grouped(dv_altim, dv_gemb, geotiles0;
    distance_from_origin_penalty = 20/100,
    ΔT_to_pscale_weight = 1,
    seasonality_weight = GlobalGlacierAnalysis.seasonality_weight,
    calibrate_to=:all,
    is_scaling_factor=Dict("pscale" => true, "ΔT" => false),
    origin_penalty_mode=:prior,
    forcing_prior=gemb_forcing_prior
)

    if ΔT_to_pscale_weight > 1 || ΔT_to_pscale_weight < 0
        error("ΔT_to_pscale_weight must be between 0 and 1")
    end

    if seasonality_weight > 1 || seasonality_weight < 0
        error("seasonality_weight must be between 0 and 1")
    end
    
    # there are issues with calibrating the dv_gemb model to individual geotiles as glaciers
    # can cross multiple geotiles. To tackle this issue one needs to process mutually 
    # exclusive geotile groups

    # loop for each geotile
    geotiles = DataFrame()
    geotiles[!, :id] = geotiles0[:, :id]
    geotiles[!, :group] = geotiles0[:, :group]
    geotiles[!, :rgi] = geotiles0[:, :rgi]
    geotiles[!, :extent] = geotiles0[:, :extent]
    geotiles[!, :area_km2] = geotiles0[:, :area_km2]
    

    geotiles[!, :pscale] .= 1.0
    # 0.0 is the no-op for a temperature offset, as 1.0 is for a scaling factor. Overwritten below for
    # every geotile whose group is fitted; it only survives where the fit is skipped.
    geotiles[!, :ΔT] .= 0.0

    dpscale = dims(dv_gemb, :pscale)
    dΔT = dims(dv_gemb, :ΔT)
 
    kwargs2 = (seasonality_weight, distance_from_origin_penalty, ΔT_to_pscale_weight=ΔT_to_pscale_weight, calibrate_to, is_scaling_factor, origin_penalty_mode, forcing_prior)

    geotile_groups = copy(geotiles.group);
    geotile_ids = copy(geotiles.id)
    groups_unique = unique(geotile_groups)

    # for full grid search
    begin
        step_size = 0.05
        s0 = -1 / minimum(dpscale.val)
        e0 = maximum(dpscale.val)
        pscale0 = s0:step_size:e0
        dpscale_search = vcat(-1 ./ pscale0[pscale0.<-1], pscale0[pscale0.>=1])

        if is_scaling_factor["ΔT"]
            s0 = -1 / minimum(dΔT.val)
            e0 = maximum(dΔT.val)
            ΔT0 = s0:step_size:e0
            dΔT_search = vcat(-1 ./ ΔT0[ΔT0.<-1], ΔT0[ΔT0.>=1])
        else
            # 0.1 K steps: fine enough to resolve the narrow pscale-ΔT valley of the cost surface before
            # `_gemb_grid_minimize` refines around the minimum
            dΔT_search = range(minimum(dΔT.val), maximum(dΔT.val); step=0.1)
        end

        dpscale_search = Dim{:pscale}(dpscale_search)
        dΔT_search = Dim{:ΔT}(dΔT_search)
    end
    # loop for each group
    # NOTE: do not do a groupby on geotiles as you need to keep the original geotiles order

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CAUSE OF DAYS LOST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # This had to be moved out of @threads loop as it is not thread safe - your guess as as good as mine [this block only take 1.5s so it's a non-issue]
    dv_altim0 = [Float32.(dropdims(sum(dv_altim[geotile=At(geotile_ids[groups_unique[i] .== geotile_groups])], dims=:geotile), dims=:geotile)) for i in eachindex(groups_unique)]
    dv_gemb0  = [Float32.(dropdims(sum(dv_gemb[geotile=At(geotile_ids[groups_unique[i]  .== geotile_groups])], dims=:geotile), dims=:geotile)) for i in eachindex(groups_unique)]

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if length(groups_unique) != 1

        pscale0 = zeros(length(groups_unique));
        ΔT0 = zeros(length(groups_unique));
 
        Threads.@threads for i in 1:length(groups_unique)
            f2 = Base.Fix{2}(Base.Fix{3}(Base.Fix{4}(gemb_altim_cost, kwargs2), dv_gemb0[i]), dv_altim0[i])
            (pscale0[i], ΔT0[i]) = _gemb_grid_minimize(f2, val(dpscale_search), val(dΔT_search))
        end

        # modify dataframe outside of loop
        geotiles[!, :grid_bound] .= ""
        for i in eachindex(groups_unique)
            gindex = groups_unique[i] .== geotile_groups;
            geotiles[gindex, :pscale] .= pscale0[i]
            geotiles[gindex, :ΔT] .= ΔT0[i]
            geotiles[gindex, :grid_bound] .= _grid_bounds_hit(pscale0[i], ΔT0[i], val(dpscale_search), val(dΔT_search))
        end

        pinned = .!isempty.(geotiles.grid_bound)
        if any(pinned)
            area_pinned = sum(sum, geotiles.area_km2[pinned])
            area_total = sum(sum, geotiles.area_km2)
            bounds = join(("$k: $(count(==(k), geotiles.grid_bound[pinned]))" for k in unique(geotiles.grid_bound[pinned])), "; ")
            @warn "$(length(unique(geotiles.group[pinned]))) of $(length(groups_unique)) geotile groups " *
                  "($(count(pinned)) geotiles, $(round(area_pinned; digits=1)) of $(round(area_total; digits=1)) km² of ice) " *
                  "are fitted at a bound of the pscale × ΔT grid, so their fit is a limit, not a minimum. " *
                  "Geotiles per bound: $(bounds)"
        end

        return geotiles[:, [:id, :extent, :group, :pscale, :ΔT, :grid_bound, :rgi]]

    else
        i = 1;
        geotiles_in_group = geotile_ids[groups_unique[i] .== geotile_groups]
        f2 = Base.Fix{2}(Base.Fix{3}(Base.Fix{4}(gemb_altim_cost, kwargs2), dv_gemb0[i]), dv_altim0[i])

        time_grid_search = @elapsed begin
            cost = fill(NaN, dpscale_search, dΔT_search)
            for pscale in dpscale_search
                for  ΔT in dΔT_search
                    cost[pscale=At(pscale), ΔT=At(ΔT)] = f2([pscale, ΔT])
                end
            end
        end

        show_plot = false
        if show_plot
            center_on_dates = DateTime(2010, 1, 1)..DateTime(2015, 1, 1)
            f = Figure();
            ax = Makie.Axis(f[1, 1]; title="geotile: $(geotiles_in_group[1])");
        
            cost0 = Inf;
            time_grid_search = @elapsed begin
                cost = fill(NaN, dpscale_search, dΔT_search)
                for pscale in dpscale_search
                    for  ΔT in dΔT_search
                        cost1 = f2([pscale, ΔT])

                        if cost0 > cost1
                            cost0 = cost1
                
                            gemb2X = gemb_dv_sample(pscale, ΔT, dv_gemb0[i])
                            gemb2X = gemb2X .- mean(gemb2X[date=center_on_dates])
                            lines!(ax, gemb2X; label = "cost: $(round(cost0, digits=1))")
                        end
                    end
                end
            end

            altim2X = dv_altim0[i] .- mean(dv_altim0[i][date = center_on_dates])

            lines!(ax, altim2X; linewidth=4, color=:black, label="altim")
            Legend(f[1, 2], ax, position=:rt);
            display(f)
        end

    
        index_minimum = argmin(cost)
        pscale_grid_search = DimPoints(cost)[index_minimum][1]
        ΔT_grid_search = DimPoints(cost)[index_minimum][2]
        println("grid_search: pscale: $(round(pscale_grid_search, digits=2)), ΔT: $(round(ΔT_grid_search, digits=2)), time: $(round(time_grid_search, digits=2))")

        # the minimizer the multi-group path uses, so this path reports the calibrated values
        (pscale_fit, ΔT_fit) = _gemb_grid_minimize(f2, val(dpscale_search), val(dΔT_search))
        println("calibrated : pscale: $(round(pscale_fit, digits=2)), ΔT: $(round(ΔT_fit, digits=2))")

        return cost, geotiles_in_group
    end
end

"""
    _gemb_grid_minimize(f, pscale_grid, ΔT_grid; stride=(4, 5), ncandidates=3, nrefine=11) -> (pscale, ΔT, cost)

Minimize `f([pscale, ΔT])` over `pscale_grid × ΔT_grid` deterministically, in three stages:

1. every `stride`-th point of the grid, covering the whole range;
2. the full grid within one coarse step of each of the `ncandidates` best coarse points;
3. an `nrefine × nrefine` grid spanning the grid neighbours of the best point found so far.

A grid search is used because the cost surface often has a long, narrow valley along which `pscale` and
`ΔT` trade off, and can have a second basin at the edge of the forcing grid; a stochastic search (ECA)
stalls partway along the valley or settles in the wrong basin, and its result then depends on the
random seed. Keeping several coarse candidates guards against two separated basins of similar cost.
"""
function _gemb_grid_minimize(f, pscale_grid, ΔT_grid; stride=(4, 5), ncandidates=3, nrefine=11)
    np, nt = length(pscale_grid), length(ΔT_grid)
    sp, st = stride
    cost(ip, it) = f([pscale_grid[ip], ΔT_grid[it]])

    # 1. coarse pass; include the last index so the grid edges are always searched
    ips = unique(vcat(1:sp:np, np)); its = unique(vcat(1:st:nt, nt))
    coarse = [(cost(ip, it), ip, it) for ip in ips, it in its]
    any(c -> isfinite(c[1]), coarse) || error("cost is not finite anywhere on the coarse pscale × ΔT grid")
    candidates = first(sort(vec(coarse); by=first), min(ncandidates, length(coarse)))

    # 2. full-resolution search around each candidate
    best = (Inf, 0, 0)
    for (_, ip0, it0) in candidates
        for ip in max(ip0 - sp, 1):min(ip0 + sp, np), it in max(it0 - st, 1):min(it0 + st, nt)
            c = cost(ip, it)
            c < best[1] && (best = (c, ip, it))
        end
    end
    isfinite(best[1]) ||
        error("cost is not finite anywhere in the full-resolution search around the $(length(candidates)) best coarse points")

    # 3. refine between the grid neighbours of the best point
    (_, ip, it) = best
    fine = (best[1], pscale_grid[ip], ΔT_grid[it])
    for p in range(pscale_grid[max(ip - 1, 1)], pscale_grid[min(ip + 1, np)]; length=nrefine),
        t in range(ΔT_grid[max(it - 1, 1)], ΔT_grid[min(it + 1, nt)]; length=nrefine)
        c = f([p, t])
        c < fine[1] && (fine = (c, p, t))
    end
    return (fine[2], fine[3], fine[1])
end

"""
    _grid_bounds_hit(pscale, ΔT, pscale_grid, ΔT_grid) -> String

Name the bounds of `pscale_grid × ΔT_grid` that a fitted `(pscale, ΔT)` lies on, joined by `","`:
any of `"pscale_min"`, `"pscale_max"`, `"ΔT_min"`, `"ΔT_max"`, or `""` for an interior fit. A fit on a
bound means the cost may still fall beyond the grid, so the value is a limit, not a minimum.
"""
function _grid_bounds_hit(pscale, ΔT, pscale_grid, ΔT_grid)
    hits = String[]
    for (name, x, grid) in (("pscale", pscale, pscale_grid), ("ΔT", ΔT, ΔT_grid))
        lo, hi = extrema(grid)
        atol = 1e-9 * (hi - lo)
        isapprox(x, lo; atol) && push!(hits, "$(name)_min")
        isapprox(x, hi; atol) && push!(hits, "$(name)_max")
    end
    return join(hits, ",")
end


"""
    gemb_calibration(path2runs; surface_mask_default="glacier", gemb_run_id, geotile_grouping_min_feature_area_km2=100, single_geotile_test=nothing, force_remake_before)

Calibrate GEMB (Glacier Energy and Mass Balance) model to altimetry data by finding optimal fits for grouped geotiles.

# Arguments
- `path2runs`: Array of file paths to binned synthesized data files
- `surface_mask_default`: Default surface mask type (default: "glacier")
- `gemb_run_id`: GEMB run identifier
- `geotile_grouping_min_feature_area_km2`: Minimum feature area in km² for geotile grouping (default: 100)
- `single_geotile_test`: Single geotile ID for testing (default: nothing)
- `force_remake_before`: Date to force file regeneration (default: nothing)

# Description
Loads GEMB data and altimetry volume change data, then finds optimal model fits for grouped geotiles.
Saves calibration results as Arrow files with geometry information, each with the tile-run provenance
of `dv_gemb` beside it ([`write_gembfit_provenance`](@ref)). An output is skipped as current only when
both files exist. Throws if `dv_gemb` records no `tile_provenance` ([`gemb_ensemble_dv`](@ref) of a
tile run).

A fit against a GEMB physics member's ensemble ([`gemb_physics_member`](@ref)) is written to that
member's folder ([`gemb_physics_path`](@ref)), so the fits of different members never share a file.
Throws if the member `dv_gemb` records disagrees with the one its tile provenance names.

The pilot sweep (GEMB run 10) covers a subset of tiles, so two defaults follow the run id `dv_gemb`
records:
- `complete_groups_only` (default `true` for run 10, `false` otherwise): fit only the geotile groups
  whose every geotile is in `dv_gemb`, so no group is fitted on part of its ice;
- `pscale_sanity_check` (default `false` for run 10, `true` otherwise): the fitted-pscale windows of
  two reference geotiles. It cannot be turned off for any other run, and with it on, a fit that does
  not include both geotiles throws.
"""
function gemb_calibration(
    path2runs,
    dv_gemb;
    geotile_width=2,
    geotile_grouping_min_feature_area_km2=100,
    single_geotile_test=nothing,
    seasonality_weight=GlobalGlacierAnalysis.seasonality_weight,
    distance_from_origin_penalty=2 / 100,
    ΔT_to_pscale_weight=1,
    force_remake_before=nothing,
    calibrate_to = :all,
    origin_penalty_mode=:prior,
    forcing_prior=gemb_forcing_prior,
    complete_groups_only=nothing,
    pscale_sanity_check=nothing
)
    pilot = get(DD.metadata(dv_gemb), "gemb_run_id", nothing) == 10
    complete_groups_only = something(complete_groups_only, pilot)
    pscale_sanity_check = something(pscale_sanity_check, !pilot)
    pscale_sanity_check || pilot ||
        error("pscale_sanity_check can be turned off only for the pilot sweep (GEMB run 10)")

    # `dv_gemb` is area-weighted by one surface mask's hypsometry, so every run must use that mask

    if !isnothing(single_geotile_test)
        # NOTE: don't prematurely trim data to `single_geotile_test` as fit optimization is done by geotile groups
        @warn "!!!!!!!!!!!!!! SINGLE GEOTILE TEST [$(single_geotile_test)], OUTPUT WILL NOT BE SAVED TO FILE  !!!!!!!!!!!!!!"
    end

    dgeotile = dims(dv_gemb, :geotile)
    params = binned_filled_fileparts.(path2runs)
    surface_masks = unique(getindex.(params, :surface_mask))
    gemb_mask = gemb_surface_mask(dv_gemb)
    surface_masks == [gemb_mask] ||
        error("dv_gemb is built on the $(gemb_mask) surface mask but path2runs uses $(join(surface_masks, ", ")); " *
              "pass a Dict of surface mask => dv_gemb to calibrate each run against its own mask")

    # The prior describes one GEMB run's forcing bias; centring the penalty on it for another run pulls
    # every fit toward the wrong model's bias.
    if origin_penalty_mode == :prior && distance_from_origin_penalty != 0
        gemb_run_id = get(DD.metadata(dv_gemb), "gemb_run_id", nothing)
        isnothing(gemb_run_id) &&
            error("dv_gemb does not record its gemb_run_id; load it with gemb_ensemble_dv")
        forcing_prior === gemb_forcing_prior && gemb_run_id != gemb_forcing_prior_run_id &&
            error("gemb_forcing_prior and distance_from_origin_penalty were derived from GEMB run " *
                  "$(gemb_forcing_prior_run_id) fits, not run $(gemb_run_id). Re-derive both from run " *
                  "$(gemb_run_id)'s unpenalized fits (distance_from_origin_penalty = 0) and update them " *
                  "with gemb_forcing_prior_run_id, or pass a forcing_prior derived for run $(gemb_run_id).")
    end

    # Each output records the tile sweep it was fitted on (`write_gembfit_provenance`), which
    # GEMB_GlacierSims' native run requires before it applies the corrections.
    tile_provenance = get(DD.metadata(dv_gemb), "tile_provenance", nothing)
    isnothing(tile_provenance) &&
        error("dv_gemb records no tile_provenance, so the calibration could not record the tile sweep " *
              "it was fitted on; load a tile-run ensemble with gemb_ensemble_dv")
    physics_member = gemb_physics_member(dv_gemb)
    _check_gemb_physics_provenance(tile_provenance["tile_run_parameters"], physics_member, "dv_gemb")

    # Windows on the fitted pscale of two reference geotiles, wide enough to pass every member of the
    # ensemble and tight enough to catch a fit that has come out wildly wrong. They are not physical
    # bounds: the precipitation bias of the GEMB forcing sets where pscale lands, so they have to be
    # re-checked against the printed values whenever the forcing changes.
    #
    # `lat[+58+60]lon[-138-136]` needs the looser window because its fitted ΔT saturates at the top of
    # the forcing grid. The cost surface is then degenerate along ΔT, which leaves pscale only weakly
    # constrained, so it scatters between local minima across ensemble members rather than tracking one
    # value. A geotile whose ΔT sits at a grid bound cannot be held to a narrow pscale.
    dgeotile_test = Dim{:geotile}([geotiles_golden_test[1], "lat[+58+60]lon[-138-136]"])
    pscale_range_test = DimStack(DimArray([1.3, 1.0], dgeotile_test; name="min"), DimArray([3.0, 2.9], dgeotile_test; name="max"))

    # Process each surface mask to load aligned geotiles
    geotiles = Dict()

    # `ΔT` is an additive air-temperature offset (no-op 0) and `pscale` is a multiplicative
    # precipitation scaling (no-op 1); neither is inferred from the sign of the data they happen to
    # be run with. A grid of `ΔT` values that were all positive, or a `pscale` grid that included 0,
    # would otherwise flip these from the sign check that used to live here.
    is_scaling_factor = Dict("pscale" => true, "ΔT" => false)

    # Restrict dv_gemb to what the altimetry and GEMB share, in both geotile and date, before anything
    # is built on its axes. GEMB spans every geotile holding ice under its own mask, which includes
    # sub-km2 slivers the binning drops, so the altimetry is the smaller set and decides.
    begin
        dv_altim0 = GVC.read_volume(GVC.volume_product_path(first(path2runs))).dv
        dgeotile_altim = dims(dv_altim0, :geotile)

        geotiles2extract = [id for id in val(dgeotile) if in(id, val(dgeotile_altim))]
        isempty(geotiles2extract) &&
            error("GEMB and the altimetry share no geotiles: $(first(path2runs))")
        if length(geotiles2extract) < length(dgeotile)
            dropped = setdiff(val(dgeotile), geotiles2extract)
            printstyled("    -> $(length(dropped)) of $(length(dgeotile)) GEMB geotiles are absent from the altimetry and take no part in the fit: $(join(first(sort(dropped), 8), ", "))\n"; color=:light_yellow)
            dv_gemb = dv_gemb[geotile=At(geotiles2extract)]
            dgeotile = dims(dv_gemb, :geotile)
        end

        ddate = dims(dv_altim0, :date)
        ddate_gemb = dims(dv_gemb, :date)
        index = .!isnan.(dv_gemb[geotile=1, pscale=1, ΔT=1])
        ex_gemb = extrema(ddate_gemb[index])
        index = .!isnan.(GVC.read_volume_height(GVC.volume_product_path(first(path2runs))).dh[geotile=1, height=1])
        ex_dv = extrema(ddate[vec(index)])
        ex = max(ex_gemb[1], ex_dv[1]), min(ex_gemb[2], ex_dv[2])
        index_dv = (ddate .>= ex[1]) .& (ddate .<= ex[2])
        index_gemb = (ddate_gemb .>= ex[1]) .& (ddate_gemb .<= ex[2])
        ddate[index_dv].val == ddate_gemb[index_gemb].val ||
            error("GEMB and altimetry dates differ over their common span $(ex): the GEMB ensemble must be binned on the altimetry date grid")
        dates2extract = ddate[index_dv].val

        dv_gemb = dv_gemb[date=At(ddate_gemb[index_gemb].val)]
    end

    for surface_mask in surface_masks # takes 7 seconds

        # Load and align geotiles with the synthesized data dimensions
        geotiles[surface_mask] = _geotile_load_align(;
            surface_mask,
            geotile_width,
            only_geotiles_w_area_gt_0=true,
            geotile_order=collect(dgeotile)
        )

        # Add in groups
        geotiles_groups = geotile_grouping(; surface_mask, min_area_km2=geotile_grouping_min_feature_area_km2, geotile_width, force_remake_before)

        # Add groups to geotiles
        _, grp_index = intersectindices(geotiles[surface_mask].id, geotiles_groups.id)
        geotiles[surface_mask][!, :group] = geotiles_groups[grp_index, :group]

        if complete_groups_only
            geotiles[surface_mask] = complete_geotile_groups(geotiles[surface_mask], geotiles_groups)
            nrow(geotiles[surface_mask]) > 0 ||
                error("$(surface_mask): no geotile group lies wholly inside the GEMB ensemble")
        end
    
        if !isnothing(single_geotile_test)
            index_geotile = findfirst(geotiles[surface_mask].id .== single_geotile_test)
            geotiles[surface_mask] = geotiles[surface_mask][(geotiles[surface_mask].group .== geotiles[surface_mask][index_geotile, :group]), :]
        end

        if false
            @warn "testing subset"
            index_test = falses(nrow(geotiles[surface_mask]))
            for geotile_test in dgeotile_test
                index_geotile_test = findfirst(geotiles[surface_mask].id .== geotile_test)
                index_test .|= (geotiles[surface_mask].group .== geotiles[surface_mask].group[index_geotile_test])
            end
        
            geotiles[surface_mask] = geotiles[surface_mask][index_test, :]
        end

        # add a single rgi column
        add_single_rgi_column!(geotiles[surface_mask])
    end

    # I'm getting an "geotile volume change contains all NaNs" when run using Threads.. no clue why
    @showprogress desc = "Calibrating GEMB model to altimetry data" for binned_synthesized_file in path2runs

        synthesized_gemb_fit = gemb_physics_path(replace(binned_synthesized_file, ".jld2" => "_gembfit.arrow"),
                                                 physics_member)

        if file_is_current(synthesized_gemb_fit, force_remake_before) &&
           isfile(gembfit_provenance_file(synthesized_gemb_fit)) && isnothing(single_geotile_test)
            printstyled("    -> Skipping $(synthesized_gemb_fit) because it was created after force_remake_before:$force_remake_before\n"; color=:light_green)
            continue
        else

            run_parameters = binned_filled_fileparts(synthesized_gemb_fit)
            surface_mask = run_parameters.surface_mask

            dv_member = GVC.read_volume(GVC.volume_product_path(binned_synthesized_file)).dv

            # Every run is fitted on the one geotile axis established above, so a run that does not
            # cover it would otherwise fail inside `At` with no indication of which file or which
            # geotile is at fault.
            absent = setdiff(val(dgeotile), val(dims(dv_member, :geotile)))
            isempty(absent) ||
                error("$(binned_synthesized_file) is missing $(length(absent)) of the $(length(dgeotile)) " *
                      "geotiles the fit is built on: $(join(first(sort(absent), 8), ", "))")

            # Volume change over the geotiles and dates the altimetry and GEMB share.
            dv_altim = dv_member[date=At(dates2extract), geotile=At(collect(dgeotile))]

            all_nans = dropdims(all(isnan.(dv_altim), dims=:date), dims=:date)

            if any(all_nans)
                error("geotile volume change contains all NaNs: $(binned_synthesized_file)")
            end

            # Find optimal fit to GEMB data
            # There are issues with calibrating the SMB model to individual geotiles since glaciers
            # can cross multiple geotiles, therefore we calibrate the model for groups of
            # distinct geotiles.
            if isnothing(single_geotile_test)

                gembfit = gemb_bestfit_grouped(dv_altim, dv_gemb, geotiles[surface_mask]; seasonality_weight, distance_from_origin_penalty, ΔT_to_pscale_weight, calibrate_to, is_scaling_factor, origin_penalty_mode, forcing_prior)

                for geotile_test in (pscale_sanity_check ? dgeotile_test : String[])
                    any(==(geotile_test), gembfit.id) || error(
                        "pscale sanity geotile $(geotile_test) is not in the fit; it is checked for " *
                        "every run but the pilot sweep (GEMB run 10)")
                    geotile_sanity_check = geotile_test
                    gt_index = findfirst(gembfit.id .== geotile_sanity_check)
                    geotile_row = gembfit[gt_index, :]
                    pscale = round(geotile_row["pscale"], digits=2)
                    ΔT = round(geotile_row["ΔT"], digits=2)
                    dv_altim_mean = dv_altim[geotile=At(geotile_row.id)]
                    dv_altim_mean = mean(dv_altim_mean[.!isnan.(dv_altim_mean)])
                    dv_gemb_mean = dv_gemb[geotile=At(geotile_row.id)]
                    dv_gemb_mean = mean(dv_gemb_mean[.!isnan.(dv_gemb_mean)])

                    println("\n$(geotile_sanity_check): pscale = $pscale, ΔT = $ΔT, dv_altim_mean = $(round(dv_altim_mean, digits=2)), dv_gemb_mean = $(round(dv_gemb_mean, digits=2)), file = $(binned_synthesized_file)")
                    pscale_min = pscale_range_test[:min][geotile=At(geotile_test)]
                    pscale_max = pscale_range_test[:max][geotile=At(geotile_test)]
                    if !(pscale_min < pscale < pscale_max)
                        error("pscale for $geotile_sanity_check is $pscale, outside the expected " *
                              "$pscale_min-$pscale_max: either the fit is wrong or the GEMB forcing " *
                              "has changed enough to move it, in which case widen pscale_range_test")
                    end
                end

                gembfit[!, :area_km2] = sum.(geotiles[surface_mask].area_km2)
                gembfit.geometry = extent2rectangle.(gembfit.extent)
                gembfit = gembfit[:, Not(:extent)]
                mkpath(dirname(synthesized_gemb_fit))
                GeoDataFrames.write(synthesized_gemb_fit, gembfit)
                write_gembfit_provenance(synthesized_gemb_fit, tile_provenance)
            else
                cost, geotiles_in_group = gemb_bestfit_grouped(dv_altim, dv_gemb, geotiles[surface_mask]; seasonality_weight, distance_from_origin_penalty, ΔT_to_pscale_weight, calibrate_to, is_scaling_factor, origin_penalty_mode, forcing_prior)
                return (cost, dv_altim, geotiles_in_group)
            end
        end
    end
end


"""
    complete_geotile_groups(geotiles, groups) -> DataFrame

The rows of `geotiles` (with `id` and `group` columns) whose group has every one of its geotiles in
`geotiles`, against the full grouping `groups` (with `id` and `group`).
"""
function complete_geotile_groups(geotiles, groups)
    present = Set(geotiles.id)
    complete = Set(g for g in unique(groups.group)
                   if all(in(present), groups.id[groups.group .== g]))
    return geotiles[in.(geotiles.group, Ref(complete)), :]
end

"""
    gemb_calibration(path2runs, dv_gemb::AbstractDict; kwargs...)

Calibrate each run in `path2runs` against the GEMB ensemble of its own surface mask: `dv_gemb` maps
surface mask (`String` or `Symbol`) to the `dv` layer of [`gemb_ensemble_dv`](@ref) for that mask.
Throws if a run's mask has no entry. Keywords are passed to the single-mask method; with
`single_geotile_test` set, returns that method's result for the first mask.
"""
function gemb_calibration(path2runs, dv_gemb::AbstractDict; kwargs...)
    for (_, runs, dv_gemb_mask) in _select_by_surface_mask(dv_gemb, path2runs; name="dv_gemb")
        result = gemb_calibration(runs, dv_gemb_mask; kwargs...)
        isnothing(result) || return result
    end
    return nothing
end



"""
    gemb_geotile_dv_file(gembinfo, surface_mask) -> String

Path of the GEMB geotile volume-change ensemble built on the hypsometry of `surface_mask`, as written by
`gemb_tiles_binning.jl` and read by [`gemb_ensemble_dv`](@ref).
"""
gemb_geotile_dv_file(gembinfo, surface_mask) =
    replace(gembinfo.filename_gemb_combined, ".jld2" => "_geotile_dv_$(surface_mask).jld2")

"""
    gemb_ensemble_built(gembinfo, surface_masks) -> DateTime

UTC modification time of the newest of `gembinfo`'s per-mask GEMB ensemble files
([`gemb_geotile_dv_file`](@ref)), on the clock [`file_is_current`](@ref) compares against. Everything
derived from the ensemble -- discharge, calibration and the calibrated synthesis -- is stale if it is
older than this, so it is the `force_remake_before` for those steps.

Throws if any of the files is missing.
"""
function gemb_ensemble_built(gembinfo, surface_masks)
    paths = [gemb_geotile_dv_file(gembinfo, surface_mask) for surface_mask in surface_masks]
    absent = filter(!isfile, paths)
    isempty(absent) ||
        error("GEMB ensemble file(s) missing: $(join(absent, ", ")); build them with gemb_tiles_binning.jl")
    return maximum(path -> Dates.unix2datetime(mtime(path)), paths)
end

"""
    gemb_surface_mask(A) -> String

The surface mask whose hypsometry the GEMB output `A` (a DimArray or one layer of the DimStack from
[`gemb_ensemble_dv`](@ref)) is area-weighted by, read from its metadata. Throws if `A` does not record one.
"""
function gemb_surface_mask(A)
    surface_mask = get(DD.metadata(A), "surface_mask", nothing)
    isnothing(surface_mask) &&
        error("GEMB output does not record the surface mask it was built on; load it with " *
              "gemb_ensemble_dv(; gemb_run_id, surface_mask)")
    return String(surface_mask)
end

"""
    _check_surface_mask(A, expected, name)

Throw unless the `"surface_mask"` in `A`'s metadata is `expected`; `name` identifies `A` in the message.
"""
function _check_surface_mask(A, expected, name)
    found = get(DD.metadata(A), "surface_mask", nothing)
    found == String(expected) ||
        error("$name is for surface mask $(repr(found)), expected $(repr(String(expected)))")
    return nothing
end

"""
    _select_by_surface_mask(by_mask::AbstractDict, path2runs; name) -> Vector

Group `path2runs` by the surface mask in each file name and pair every group with its entry in
`by_mask`, whose keys may be `String`s or `Symbol`s. Returns `(surface_mask, runs, by_mask[surface_mask])`
tuples, one per mask, in order of first appearance. Throws, naming `name`, if a mask has no entry.
"""
function _select_by_surface_mask(by_mask::AbstractDict, path2runs; name)
    by_mask = Dict(String(k) => v for (k, v) in by_mask)
    masks = [String(binned_filled_fileparts(f).surface_mask) for f in path2runs]
    absent = setdiff(unique(masks), keys(by_mask))
    isempty(absent) ||
        error("$name has no entry for surface mask $(join(absent, ", ")), used by path2runs; " *
              "it holds $(join(sort(collect(keys(by_mask))), ", "))")
    return [(mask, path2runs[masks .== mask], by_mask[mask]) for mask in unique(masks)]
end

"""
    gemb_ensemble_dv(; gemb_run_id=9, surface_mask, physics_member=nothing)

Load the GEMB geotile volume-change ensemble of run `gemb_run_id`, area-weighted by the hypsometry of
`surface_mask`, and add the derived variables (see [`gemb_add_derived_vars!`](@ref)).

# Arguments
- `gemb_run_id`: Identifier for the GEMB run.
- `surface_mask`: Surface mask the ensemble was built on, e.g. `:glacier` or `:glacier_rgi7`. Altimetry
  members are calibrated against the ensemble of their own mask.
- `physics_member`: the GEMB physics member of a run that has them (see [`gemb_info`](@ref)).

# Returns
- DimStack over (geotile, date, pscale, ΔT) in km³ of ice equivalent. Every layer records
  `surface_mask` (see [`gemb_surface_mask`](@ref)), `gemb_run_id` and `gemb_physics_member` (see
  [`gemb_physics_member`](@ref)) in its metadata, and a tile run also the `tile_provenance` the file
  was saved with ([`gemb_tile_provenance`](@ref)). That throws if the file records none, was built
  from a tile directory other than the run's `tile_dir`, or from tiles of another physics member.

# Examples
```julia
julia> dv_gemb = gemb_ensemble_dv(; gemb_run_id=8, surface_mask=:glacier_rgi7)
julia> smb = dv_gemb[:smb]
```
"""
function gemb_ensemble_dv(; gemb_run_id=9, surface_mask, physics_member=nothing)
    gembinfo = gemb_info(; gemb_run_id, physics_member)
    gemb_geotile_filename_dv = gemb_geotile_dv_file(gembinfo, surface_mask)
    file_surface_mask, dv_gemb = FileIO.load(gemb_geotile_filename_dv, "surface_mask", "gemb_dv")
    String(file_surface_mask) == String(surface_mask) ||
        error("$(gemb_geotile_filename_dv) was built on the $(file_surface_mask) surface mask, not $(surface_mask)")

    # `gemb_add_derived_vars!` builds a new object via `merge` rather than mutating in place, so
    # its return value has to be captured. Discarding it left :smb, :runoff and :dv absent, and
    # the `dv_gemb[:smb]` usage shown in this function's own docstring would throw a KeyError.
    dv_gemb = gemb_add_derived_vars!(dv_gemb)

    member = get(gembinfo, :physics_member, nothing)
    metadata = Dict{String,Any}("surface_mask" => String(surface_mask), "gemb_run_id" => gemb_run_id,
                                "gemb_physics_member" => member)
    if haskey(gembinfo, :tile_dir)
        provenance = _gemb_ensemble_tile_provenance(gemb_geotile_filename_dv, gembinfo.tile_dir)
        _check_gemb_physics_provenance(provenance["tile_run_parameters"], member, gemb_geotile_filename_dv)
        metadata["tile_provenance"] = provenance
    end
    return DD.maplayers(A -> DD.rebuild(A; metadata), dv_gemb)
end

# The `gemb_tile_provenance` a tile-run ensemble file was saved with, checked against the tile
# directory the run reads.
function _gemb_ensemble_tile_provenance(path, tile_dir)
    provenance = JLD2.jldopen(path, "r") do f
        haskey(f, "tile_provenance") ? f["tile_provenance"] : nothing
    end
    isnothing(provenance) &&
        error("$(path) records no tile_provenance, so the tile sweep it was built from is unknown; " *
              "rebuild it with gemb_tiles_binning.jl")
    isdir(tile_dir) && provenance["tile_dir"] == realpath(tile_dir) ||
        error("$(path) was built from the tile files in $(provenance["tile_dir"]), not $(tile_dir)")
    return provenance
end

"""
    gembfit_provenance_file(gembfit_file) -> String

Path of the tile-run provenance written beside a `gemb_calibration` output: `.arrow` replaced by
`_provenance.toml`, the file GEMB_GlacierSims' `forcing_corrections_provenance_path` names.
"""
function gembfit_provenance_file(gembfit_file)
    endswith(gembfit_file, ".arrow") || error("not an .arrow calibration output: $(gembfit_file)")
    return replace(gembfit_file, r"\.arrow$" => "_provenance.toml")
end

"""
    write_gembfit_provenance(gembfit_file, tile_provenance) -> String

Write [`gembfit_provenance_file`](@ref)`(gembfit_file)` from `tile_provenance` (a
[`gemb_tile_provenance`](@ref)), in the format GEMB_GlacierSims' `assert_corrections_fitted_on` reads:

```toml
tile_dir = "/mnt/bylot-r3/data/era5land/tile_runs_v2"

[tile_run_parameters]
software_gemb_climateforcing_tree = "…"
software_gemb_glaciersims_tree = "…"
software_gemb_tree = "…"
software_gemb_version = "…"
# … every other software_* parameter the tile files record
```

The file is written beside the output and then moved into place, so a reader never sees it half
written. Returns its path.
"""
function write_gembfit_provenance(gembfit_file, tile_provenance)
    tile_dir = tile_provenance["tile_dir"]
    parameters = tile_provenance["tile_run_parameters"]
    tile_dir isa AbstractString || error("tile_provenance has no tile_dir")
    absent = [k for k in GEMB_TILE_PROVENANCE_REQUIRED if !haskey(parameters, k)]
    isempty(absent) || error("tile_provenance does not record $(join(absent, ", "))")

    path = gembfit_provenance_file(gembfit_file)
    partial = path * ".partial"
    open(partial, "w") do io
        TOML.print(io, Dict{String,Any}("tile_dir" => String(tile_dir),
                                        "tile_run_parameters" => Dict{String,Any}(parameters)); sorted=true)
    end
    mv(partial, path; force=true)
    return path
end

"""
    gemb_add_derived_vars!(dv_gemb)

Return a copy of `dv_gemb` with derived variables (smb, runoff, dv) added.

Computes smb = acc - melt + refreeze - ec, runoff = melt - refreeze, dv = smb + fac.

# Arguments
- `dv_gemb`: Symbol-keyed collection of GEMB variables supporting `merge` -- a NamedTuple or
  DimStack, not a `Dict` -- and containing acc, melt, refreeze, ec, fac

# Returns
- A new collection with :smb, :runoff and :dv added

!!! warning
    Despite the `!`, this does **not** mutate its argument: it is built on `merge`, which returns
    a new object. Always assign the result -- `dv_gemb = gemb_add_derived_vars!(dv_gemb)`.

# Examples
```julia
julia> dv_gemb = gemb_add_derived_vars!(dv_gemb)
julia> dv = dv_gemb[:dv]
```
"""
function gemb_add_derived_vars!(dv_gemb)
    dv_gemb = merge(dv_gemb, (; smb = dv_gemb[:acc] .- dv_gemb[:melt] .+ dv_gemb[:refreeze] .- dv_gemb[:ec], runoff = dv_gemb[:melt] .- dv_gemb[:refreeze]))
    dv_gemb = merge(dv_gemb, (; dv = dv_gemb[:smb] .+ dv_gemb[:fac]))

    return dv_gemb
end

"""
    gemb_dv_sample(pscale, ΔT, dv_gemb)

Sample GEMB volume change at a single (pscale, ΔT) by linear interpolation.

# Arguments
- `pscale`: Precipitation scaling factor
- `ΔT`: Air temperature offset applied to the GEMB forcing (K)
- `dv_gemb`: DimArray with dimensions (date, pscale, ΔT)

# Returns
- DimArray of volume change with dimension :date only.

# Examples
```julia
julia> dv_sampled = gemb_dv_sample(1.0, 0.0, dv_gemb)
```
"""
function gemb_dv_sample(pscale, ΔT, dv_gemb)

    ddate = dims(dv_gemb, :date)
    dpscale = dims(dv_gemb, :pscale)
    dΔT = dims(dv_gemb, :ΔT)

    itp = Interpolations.interpolate((eachindex(ddate), dpscale.val, dΔT.val), dv_gemb.data, (NoInterp(), Gridded(Linear()), Gridded(Linear())))
    gemb_dv_out = DimArray(itp.(eachindex(ddate), Ref(pscale), Ref(ΔT)),ddate)

    return gemb_dv_out  
end







"""
    dv_adjust4discharge!(gemb, discharge)

Adjust GEMB volume change (dv) for discharge in-place by subtracting cumulative discharge over time.

# Arguments
- `gemb`: DimStack from `gemb_ensemble_dv` containing :dv (and :date); :dv is modified in-place
- `discharge`: geotile discharge from `global_discharge_filled`, a rate in km³ ice equivalent per year
  applied over decimal years. Its `"surface_mask"` metadata must match `gemb`'s, or this throws.

# Returns
- The modified gemb (gemb[:dv] updated).

# Examples
```julia
julia> dv_adjust4discharge!(gemb, discharge)
julia> # gemb[:dv] now has discharge removed
```
"""
function dv_adjust4discharge!(gemb, discharge)
    _check_surface_mask(discharge, gemb_surface_mask(gemb[:dv]), "discharge")

    decyear = decimalyear.(dims(gemb, :date)) .- decimalyear.(dims(gemb, :date)[1]);
    gemb[:dv][:] = @d gemb[:dv] .- (decyear .* discharge[:discharge]);
    return gemb
end

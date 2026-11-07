# Methods record: altimetry seasonal cycle and GEMB calibration (October 2026)

Changes and findings that the final paper's Extended Methods must reflect. Each item states what the
pipeline does, the evidence for it, and its status. The analysis scripts and per-tile results cited
under "Evidence" are in `/home/gardnera/overnight_2026-10-06/` (and its `evidence/` subfolder) on the
JPL servers; they are not part of the package.

All numbers are for the reference ensemble member
(`glacier_rgi7_dh_cop30_v2_cc_nmad5_v01`, fill set 2 for "production", fill set 6 for the change)
with GEMB run 8, unless stated otherwise.

## 1. The filled altimetry seasonal cycle was distorted by the residual smoothing

**Finding.** `hyps_model_fill!` replaces every elevation-time bin, observed or not, with the fill model
(`model1`, a single annual sine with one phase per geotile) plus the model residual smoothed by the
median of its nearest neighbours in (time, elevation / `smooth_h2t_length_scale`). For ICESat-2
(`smooth_n = 5`) with `smooth_h2t_length_scale = 400`, those neighbours are the bin and ±2 months at
the same elevation: a 5-month running median of the residual. It removes the part of the seasonal
cycle the annual sine cannot represent, which makes the filled cycle more symmetric and smaller than
observed: the seasonal maximum moves earlier and the minimum later.

Because the other missions' seasonal cycles are replaced by ICESat-2's sine-only seasonal model during
amplitude normalization, the distortion also reaches the synthesized product before 2018.

**Evidence.**
- Raw ICESat-2 bins vs filled, 10 strong-signal tiles, per elevation band: filling moves the seasonal
  maximum earlier in all 10 (median ≈ −15 days) and the minimum later in all 10 (≈ +9 days).
- Known-signal injection: an asymmetric 2 m cycle added to the raw ICESat-2 bins of 100 glacier
  geotiles (only where ICESat-2 observed, so sampling and noise are real) is recovered by the production
  fill at 0.85 of its range (IQR 0.82–0.87), with the maximum 16 days early and the minimum 19 days late.
- The apparent seasonal phase offset between altimetry and GEMB is mostly this artifact: across 298
  geotiles, observed-minus-GEMB seasonal-maximum lag of −18.5 days (median) falls to −7.5 days once an
  equivalent operator is applied to GEMB.
- Timestamps are consistent between the two products: both represent the surface height at the centre
  of each 30-day bin. The only systematic timing offset is the ≤ 9-day drift of the decimal-year bin
  edges documented in `project_decyear_bins`.

## 2. Fix: hold a residual climatology out of the smoothing (fill sets 5–8)

**Method.** For ICESat-2 and ICESat, `hyps_model_fill!` fits, in each elevation bin with at least 3
observations in every calendar month, the first two annual harmonics to the model residual
(`preserve_residual_climatology`, `residual_climatology_harmonics = 2`). This climatology is subtracted
before the nearest-neighbour smoothing and added back afterwards, so the smoothing still removes
month-to-month sampling noise but not the mean seasonal shape. Bins without full monthly coverage are
smoothed as before. `hyps_amplitude_normalize!` (`transfer_residual_climatology`) carries ICESat-2's
climatology, together with its sine, into the other missions. Fill sets 5–8 in
`binned_filling_parameters` are sets 1–4 with these options enabled; sets 1–4 are unchanged.

Alternatives tested and rejected:
- adding a semiannual harmonic to the calibration cost (amplitude of the residual): worse fits in about
  half of the golden geotiles, because GEMB's parameters cannot produce a second annual cycle;
- keeping observed residuals unsmoothed: restores timing but increases month-to-month noise 5×;
- smoothing across elevation instead of time (`smooth_h2t_length_scale` up to 4000): ICESat-2 sampling
  noise is shared across elevations within a month, so noise rises 4×;
- a 12-bin (median) climatology: noise 1.9× production; a 3-harmonic climatology: noise 1.13×.

**Evidence (100 geotiles, stratified by RGI region, all 19 regions).**
- Injection test: set 6 recovers 0.95 of the known range (IQR 0.90–0.96), maximum error −1 day,
  minimum error +5 days; it recovers more than 110 % of the range in 0 of 98 tiles (no inflation).
- Noise and gap filling, paired against production: month-to-month noise ratio 1.00 (IQR 0.99–1.03);
  hold-out gap-fill RMSE ratio 1.00 (IQR 1.00–1.03).
- Stable snow-free land (41 geotiles within ±30°): the climatology adds 0.02 m of spurious seasonal
  range, against glacier seasonal ranges of 1–150 m.
- The climatology is built only where coverage allows: no effect in the lowest-coverage third of
  geotiles (13–69 % of months observed).
- Full chain (fill → normalization → synthesis), production code, all geotiles: set 2 reproduces the
  production files exactly; set 6's synthesized seasonal timing is closer to raw ICESat-2 in 65 % of
  tiles, in both 2018–2025 and 2003–2018, with an unchanged noise level. The synthesized seasonal range
  is 1.2× that of set 2.

**Status.** Adopted. All 192 members of fill sets 5–8 were filled and synthesized on 2026-10-06/07
(`*_p5`–`*_p8_*` files next to the set 1–4 files; synthesis error file
`geotile_synthesis_error_fillsets5to8.jld2`). `run_all.jl` uses sets 5–8 and that error file, and
`reference_ensemble_file` is the set-6 member.

## 3. Calibration optimizer: ECA replaced by a staged grid search

**Finding.** `gemb_bestfit_grouped` minimized each group's cost with ECA (`η_max = 1.0, K = 6`,
3000 calls). On the cost surface's long, narrow `pscale`–`ΔT` valley the population collapses early and
stops partway along it; where it stops depends on the random seed, and more calls do not help. For
`lat[-34-32]lon[-070-068]` the stored fit had 12 % (set 2) and 52 % (set 6) higher cost than the true
minimum. ECA with default settings converges more often but still lands in the wrong basin when the
minimum is at the edge of the forcing grid (1 of 30 groups, +60 % cost).

**Method.** `_gemb_grid_minimize` searches deterministically: a coarse pass over the whole range
(every 4th `pscale` and 5th `ΔT` grid point), the full grid (0.05 in `pscale` on its symmetric axis,
0.1 K in `ΔT`) around the three best coarse points, then an 11 × 11 refinement between the grid
neighbours of the best point. About 1000 cost evaluations per group.

**Evidence (30 groups).** Staged grid: 0.085 s per group, cost equal to or lower than the exhaustive
dense grid in 29 of 30 (worst +0.013 %). ECA (defaults): 0.38 s, one failure. Exhaustive dense grid:
0.84 s. A full member calibrates in 1.1 min on 112 threads (ECA: 1.5 min). Against the stored ECA fits,
10 % of ice area moves by more than 0.1 in `pscale` or 0.25 K in `ΔT`.

**Status.** Adopted. Sets 5–8 recalibrated with it on 2026-10-07 (together with item 4). The set 1–4
`*_gembfit.arrow` and `*_gembfit_dv.jld2` files are still the ECA fits.

## 4. Distance-from-origin penalty replaced by a centred empirical prior

The penalty multiplies the cost by `1 + d × wd`, where `d` is the distance of (`pscale`, `ΔT`) from
(1, 0) and `wd = distance_from_origin_penalty = 0.70` (with `ΔT_to_pscale_weight = 0.5`). Study on the
reference member, set 6, 672 groups, staged grid optimizer:

- **Identifiability without the penalty.** 497 groups (88 % of ice area) have a well-defined minimum
  (within 5 % of the minimum cost, sd of log `pscale` < 0.15 and sd of `ΔT` < 0.5 K).
- **Empirical prior from those groups (area-weighted).** `pscale` centred at 1.60 (×/÷ 1.57), `ΔT`
  centred at 2.20 K (sd 1.65 K), correlation(log `pscale`, `ΔT`) 0.50. The penalty's origin (1, 0) lies
  well outside the bulk of the well-constrained fits.
- **Temporal cross-validation** (fit 2000–2013, score 2014–2025, and the reverse; score = held-out cost
  / best achievable): `wd` = 0, 0.1, 0.35, 0.7, 1.5 give 1.537, 1.529, 1.525, 1.585, 1.972. The optimum is
  broad between 0.1 and 0.35; the current 0.7 predicts held-out years worse than no penalty.
- **Effect on runoff.** Total runoff (this member): 921.6 km³ yr⁻¹ at `wd` = 0, 907.4 at 0.35, 886.7 at
  0.7 (−3.8 % relative to no penalty), 818.2 at 1.5. At 0.7, 59 % of ice area moves noticeably from its
  unpenalized fit.
- 215 groups (8 % of ice area) fit at the edge of the forcing grid without the penalty; these are the
  cases a penalty is needed for.

**Centred prior (`origin_penalty_mode = :prior`).** The penalty becomes `cost × (1 + wd × d)` with `d`
the Mahalanobis distance of (log `pscale`, `ΔT`) from the empirical prior (`gemb_forcing_prior`), so
the relative scaling of the two parameters and their correlation come from the data and
`ΔT_to_pscale_weight` is not used. The multiplicative form is kept so the penalty acts the same on groups
whose cost differs by orders of magnitude.

- The prior is stable between halves of the record: `pscale` 1.69 / 1.55, `ΔT` 2.34 / 2.06 K,
  correlation 0.36 / 0.61 (2000–2013 / 2014–2025).
- Cross-validation, each fold using a prior estimated from its training half only: `wd` = 0, 0.1, 0.2,
  0.35, 0.5, 0.7, 1.0, 1.5, 2.5 give 1.537, 1.531, 1.538, 1.545, 1.558, 1.594, 1.714, 1.965, 2.301. A
  stronger penalty does not improve prediction of held-out years, including for poorly constrained
  groups; where the cost surface is flat the penalty selects a point rather than adding information.
- Full record: total runoff 921.6 (`wd` = 0), 936.5 (0.35), 938.9 (0.7), 915.4 km³ yr⁻¹ (2.5); fits on
  the edge of the forcing grid fall from 215 groups (8 % of area) to 99 (0.35), 62 (0.7) and 7 (2.5).
  About 90–95 % of the poorly constrained area moves for `wd` ≥ 0.1, while most well-constrained fits stay.
  Centring the prior removes the legacy penalty's downward runoff bias (−3.8 % at `wd` = 0.7), and runoff
  varies by < 0.3 % for `wd` between 0.35 and 1.0.

**Adopted.** `origin_penalty_mode = :prior` (now the default) with `wd = distance_from_origin_penalty =
0.35`: the largest value within 1 % of the best held-out skill; it halves the number of grid-edge fits and
sits on the stable part of the runoff curve. Sets 5–8 recalibrated with it on 2026-10-07. For the
reference member, area-weighted `pscale` 1.65 and `ΔT` 2.00 K (ECA with the legacy penalty: 1.48,
1.67 K).

## 5. Cost-function options tested and not kept in the code

- An orthogonal cost (trend-only RMS in place of the full-residual RMS), a composite cost (comparing
  median seasonal cycles), an additive quadratic prior centred on (1, 0), weighting residuals by the
  altimetry error, and evaluating the seasonal term over the ICESat-2 era only. None improved the fits
  enough to keep; the composite term was worse: comparing residual seasonal cycles penalizes any timing
  offset, so the fit shrinks the modelled seasonal amplitude and lowers runoff. The code keeps only the
  legacy cost with the `:prior` (default) and `:legacy` penalties.
- The single-group diagnostic path of `gemb_bestfit_grouped` evaluates `ΔT` at 0.1 K steps
  (was 0.45 K) and reports the same minimizer the production path uses.

## 6. Effect on the results: new ensemble (sets 5–8) vs the set 1–4 ensemble

Both ensembles use GEMB run 8; the set 1–4 fits are the ECA fits with the legacy penalty, so this compares
the combined effect of items 2–4. Trends 2000-03 to 2024-12, reference member ± 95 % ensemble error
(× 1.5), units as stored in the `*_gembfit_dv` files:

| global | sets 1–4 | sets 5–8 | change |
|---|---|---|---|
| mass change `dm` | −299.4 ± 29.4 | −302.3 ± 27.6 | −1.0 % |
| runoff | 819.1 ± 30.3 | 928.6 ± 23.2 | +13.4 % |
| accumulation | 680.6 ± 37.2 | 786.2 ± 22.3 | +15.5 % |
| melt | 1009.5 ± 30.7 | 1134.6 ± 27.0 | +12.4 % |
| refreeze | 190.4 ± 4.4 | 206.0 ± 5.1 | +8.2 % |

- The altimetry trend is unchanged; the larger, better-timed seasonal cycle routes more mass through
  the glaciers (more accumulation and more runoff for the same net loss). Runoff rises in 17 of 19
  regions, most in the Southern Andes (+39 %), South Asia West (+20 %), Russian Arctic (+16 %), High
  Mountain Asia (+14 %) and Alaska (+11 %).
- Ensemble spread of global runoff narrows from ±30.3 to ±23.2: part of the earlier spread was ECA's
  seed-dependent error.
- Seasonal cycle, ICESat-2 era, reference member, model `dv` vs altimetry: global seasonal misfit 0.08 →
  0.03 of the observed range, model/observed range 1.09 → 0.96, timing of the seasonal maximum/minimum
  −23/+9 d → +2/−6 d. The misfit improves in 14 of 21 regions; it worsens in the Southern Andes,
  Antarctic periphery, W Canada and Central Europe.
- GRACE (regional_results.jl's comparison: RGI 1–4, 6–12, 16–18 and HMA, endorheic correction,
  2002-05 to 2024-12): RMSE 2.25 → 2.21 Gt yr⁻¹ for the reference member; ensemble median 2.26 → 2.24.
  GRACE constrains the mass trend, which is unchanged, not the accumulation/runoff partition.

## 7. Fixes in the regional aggregation

- `region_min_frac_error!` assigned exactly 20 % of each component rate as its error instead of
  max(ensemble error, 20 %); it now applies the maximum, as the Methods describe. Component-rate errors
  change wherever the ensemble spread exceeds 20 % of the value.
- `runs2rgi` leaves out geotiles with no altimetry at any date (sub-km² ice fragments GEMB covers but the
  binning drops; five, 0.18 km² in total, in the RGI6-mask members); they previously made RGI 8, 10, 11,
  19 and the global sum NaN in those members. It errors if such geotiles hold ≥ 1 km² of glacier.

## 8. GEMB pipeline physics corrections (not yet in any run)

Status: implemented on local branches (GEMB.jl `fix/physics-audit`, GEMB_ClimateForcing
`fix/forcing-audit`, GEMB_GlacierSims `fix/pipeline-audit`, GGA `fix/gemb-audit`), each with
analytical tests; none merged and no GEMB run uses them. Every GEMB number above (run 8) predates them.
The review behind them, with evidence and magnitudes, is `gemb_pipeline_review_2026-10-13.md`. Once a
run uses them, the Extended Methods must describe the model as follows, with before/after calibration
numbers added here:

- **Surface energy balance.** The surface absorbs `ε·LW↓`; upward longwave is emitted plus reflected
  (`known_issue_gemb_longwave_absorption.md`). Stable and unstable turbulent fluxes come from
  Monin–Obukhov similarity with Beljaars–Holtslag (stable) and Paulson (unstable) functions, solved for
  the Obukhov length by four fixed Newton steps from the bulk Richardson number at consistent
  measurement heights. There is no pressure (Exner) factor on the sensible flux. Vapor is exchanged in
  the phase whose latent heat is charged (liquid at L_v at a wet surface, ice at L_s otherwise). A cell
  holding liquid water stays at 0 °C in the thermal solve until that water has refrozen. The explicit
  sub-step is also bounded by the surface-cell relaxation time. Ice saturation vapor pressure is Murphy
  & Koop (2005).
- **Albedo and shortwave.** Exposed ice takes the band's MODIS bare-ice albedo (`albedo_ice`) exactly;
  firn is linear between it and 0.65; snow thinner than 10 cm blends linearly toward the surface below.
  Solar zenith angle is the SW-weighted mean over each forcing hour from solar geometry (Meeus 1991 as in
  NOAA's calculator). The diffuse fraction comes from the clearness index (Erbs et al. 1982), and cloud
  optical thickness from inverting t = 1/(1 + 0.13τ) (Lacis & Hansen 1974) against the FAO-56 clear-sky
  transmissivity, capped at 100. The Lacis & Hansen coefficient was not checked against the paper. Every
  shortwave branch absorbs `(1−α)(SW−SW_d) + (1−α_d)·SW_d`.
- **Mass.** Rain enters the surface cell's pore water, with its sensible heat added to the matrix, so it
  does not change snow density; rain that freezes is counted as refreeze. Dry surface snow grows by
  Marbouty and firn by Arthern (`:hybrid`).
- **Downscaling.** The fitted lapse rate is the ambient (ice-free) slope. Each reanalysis cell is brought
  to the ice-free state with the decoupling fit's mixing model `E_obs = [1 − glm·(1 − k)]·E_ambient`,
  lapsed, then decoupled with the band's full k, so k is applied once. The decoupling fit uses only cells
  above melting. Decoupling leaves LW↓ unchanged and caps vapor pressure at saturation of the cooled
  air. Humidity uses one saturation curve, the IFS curve over water that defines ERA5 dewpoint, and
  temperature and elevation adjustments hold e/e_s(T), so ΔT = 0 and Δz = 0 leave the forcing exactly
  unchanged.
- **Spin-up.** The ceiling is 500 years. Unconverged columns warn and are counted (number and glacier
  area) in each tile file; they do not fail the run.
- **Calibration.** GEMB is weighted by the same glacier mask (RGI6 or RGI7) as the altimetry it is
  calibrated against. The perturbation grid adds pscale 2, 2.5 and 3. Fits pinned at a grid bound are
  labelled and warned about.
- **Discharge.** The geotile discharge term is converted from mass to volume once, and Antarctic
  SMB-derived discharge is converted to Gt with δice (Antarctic total 110.6 → 101.4 Gt/yr on the current
  discharge file).
- **Limitations kept.** Snowfall below the minimum cell thickness resets the surface cell's grain
  properties. The bare-ice albedo product is the darkest 5 % of retrievals.

# GEMB pipeline review, 2026-10-13

Review of the code that produced the paper's GEMB results: GEMB.jl `0cda823`, GEMB_ClimateForcing
`352d7db` (installed copy `~/.julia/packages/GEMB_ClimateForcing/1KYoZ`), GEMB_GlacierSims working
tree on `feat/native-resolution-run`, and the GEMB-consuming parts of this repository. The known
longwave absorption issue (`known_issue_gemb_longwave_absorption.md`) is excluded.

This file records the code as it was reviewed. The fixes and the analytical tests specified below are
on the branches GEMB.jl `fix/physics-audit`, GEMB_ClimateForcing `fix/forcing-audit`,
GEMB_GlacierSims `fix/pipeline-audit` and GlobalGlacierAnalysis `fix/gemb-audit`; the model as it
stands on those branches is described in `methods_2026-10_seasonal_cycle_and_calibration.md` §8.

"Verified" means the mechanism was confirmed by reading the code a second time, independently of
the reviewer. "Reported" means one reviewer's finding, with the evidence they gave; magnitudes from
synthetic forcings show sensitivity, not production values. Scratch scripts behind the evidence are
in `/tmp/gemb_review/<area>/` (not durable).

## Production configuration

`initialize_parameters` defaults except `albedo_ice` per band: `:GardnerSharp` albedo,
`:Arthern` grain growth and densification, explicit thermal solver, constant c_p, Calonne2019
conductivity, ε = 0.97, no subsurface shortwave, `melt_geometry = :thickness`,
`runoff_method = :instantaneous`, 250 m column. Forcing never sets `solar_zenith_angle`,
`cloud_optical_thickness` or `shortwave_downward_diffuse` (all 0). Hourly ERA5-Land, z_T = 2 m,
z_w = 10 m.

## Findings that bias results

| # | Where | Finding | Status | Size |
|---|---|---|---|---|
| 1 | GEMB `calculate_albedo.jl:111` | Bare-ice albedo is `α_ice + (0.58 − α_ice)·exp(−M/200)` with `M` = previous step's surface melt + `water[1]` (a few kg m⁻² hourly), so ice albedo stays at ≈0.58 and the per-band MODIS `albedo_ice` is nearly inert | Verified | absorbed SW on bare ice ~35–40 % low where MODIS ≈ 0.30–0.35 |
| 2 | ClimateForcing `derive_lapse_rate` + GlacierSims `_accumulate_cell_intervals!` | The fitted lapse rate is the on-glacier slope k·b; it is applied to donor temperature and then decoupled again with k_eff, so the vertical term is k_eff·k·b·Δz instead of k·b·Δz | Verified | +0.7 to +2.4 K at bands 1–2.5 km above donor (synthetic, k = 0.8); grows with height, flattens mass-balance gradient |
| 3 | GEMB `calculate_accumulation.jl:204-221` + `remove_melt!` | Rain joins the solid matrix at fixed `dz` (density rises), then melts out under `:thickness` keeping the raised density: rain permanently compacts surface snow, can push it past 830 kg m⁻³ (ice-albedo and impermeable branches). Inherited from MATLAB | Verified | synthetic maritime: runoff −9 to −21 %, refreeze +26–33 %, SMB +0.4 m w.e. yr⁻¹ when fixed |
| 4 | GEMB `turbulent_heat_flux.jl:116-121` | Stable branch: `zL = Ri/(1−5Ri)` below Ri = 0.2 (→ ∞), `zL = Ri` above: flux collapses then jumps back to near-neutral. Ri uses effective height z_w²/z_T = 50 m, so ~73 % of stable melt-season hours are in the `zL = Ri` branch | Verified (code) | SHF −21 to −38 %, melt −0.5 to −12 % with continuous alternatives |
| 5 | GEMB forcing defaults; no producer upstream | Snow albedo computed for overhead sun, clear sky, no diffuse | Verified (two reviewers) | snow albedo 0.03–0.12 too low; partly offsets #1 in calibration but distorts melt with elevation |
| 6 | GGA `utilities_synthesis.jl:1338` | `discharge` from `global_discharge_filled` is already km³ (`mass2volume=true`, line 1088) and is divided by `volume2mass` again | Verified | geotile `discharge` column ×1000/δice too large |
| 7 | GGA `utilities_synthesis.jl:790, 1056` (`discharge2smb`) | Antarctic SMB-derived discharge is km³ ice yr⁻¹ stored as `discharge_gtyr` | Reported | Antarctic discharge ~9 % high (7–10 Gt/yr); 80–115 Gt/yr gate checks the inflated value |
| 8 | GGA `gemb_tiles_binning.jl:33` vs `gemb_calibration` | GEMB dv built with RGI6 hypsometry, calibrated against RGI7-weighted altimetry for `glacier_rgi7_*` members (incl. reference) | Reported | RGI7/RGI6 area ratio 0.67–1.22 per geotile, 0.93–1.08 area-weighted (5–95 %); 7.3 % of ice area in geotiles off by > 10 %; 27 geotiles have ice only in RGI6 (`scripts/audit/discharge_units_and_surface_masks.jl`) |
| 9 | GEMB `turbulent_heat_flux.jl:160`, `gemb_core.jl:66-75` | At a melting surface latent flux uses L_v but vapor mass is added to/removed from ice: L_f·EC energy created/destroyed; energy check shares the assumption. Inherited from MATLAB | Reported (two reviewers) | synthetic SMB −25 to −78 kg m⁻² yr⁻¹ when fixed |
| 10 | ClimateForcing `climate_adjustment.jl:126-135`, `utils.jl:114-133` | ΔT = 0 is not an identity: after decoupling (T cooled, e kept) and in cold saturated air (dewpoint over water, RH against ice, Bolton vs Tetens), the RH clamp removes vapor and LW | Reported (three reviewers) | up to −316 Pa / −6.9 W m⁻² in humid melt season; −1.6 to −2.9 W m⁻² cold saturated |
| 11 | ClimateForcing `downscaling_parameters.jl:258-327, 571-589` | k and Γ fits biased when the melting level falls inside the region (truncated excess fitted as linear) | Reported | k 0.741 vs 0.80, Γ 4.48 vs 4.80 K/km (synthetic) |
| 12 | Spin-up (GlacierSims data) | Spin-ups unconverged at the 100-year ceiling: 11 % of glacier area at baseline, 41 % at −3 K; 565/801 tiles have ≥1 unconverged column | Reported (data) | ~2.2 cm/yr firn-air drift over unconverged area; 3–4 cm/yr tile means in N Greenland, Ross Sea, Alaska |
| 13 | GlacierSims `gemb_tile_sweep.jl:142-147`, GGA `gemb_dv_sample` | Perturbation grid gaps (pscale 1.5→4, ΔT 1→3) where most fits land; linear interpolation of a concave response; no warning for fits pinned at bounds | Reported | ~0.15–0.25 K equivalent mid-gap; ~2 % of area pinned at ΔT +6, ~2 % at pscale 4 |
| 14 | GEMB `calculate_accumulation.jl:196-198` | Any snowfall below `dzmin` resets the whole top cell to fresh-snow grains (MATLAB mass-weights) | Reported | synthetic runoff +1–7 % with mass weighting; truth is between the two |
| 15 | GEMB `calculate_grain_size.jl:200-207` | Default `:Arthern` firn grain-growth law applied to surface snow: dry snow barely coarsens | Reported (physics choice) | r_e 0.14–0.21 mm after 60 d vs 0.30–0.46 (Marbouty) |
| 16 | GEMB `turbulent_heat_flux.jl:102, 156` | Exner factor multiplies SHF and Ri, not LHF | Reported | +5–19 % SHF before stability; ~3 % in full run |
| 17 | GEMB thermal solver | Liquid water in surface cell not used as a latent buffer; wet cell cools below 0 °C within a step | Reported (structural) | not quantified |
| 18 | ClimateForcing `downscaling_parameters.jl:750`, `glacier_adjustment.jl:162` | Decoupling also lowers LW↓ (σT′⁴), not supported by what the k fit measures | Reported (design) | ~4 W m⁻² per K of decoupling |
| 19 | ClimateForcing `pooled_ice_albedo.jl:27-31` | Bare-ice albedo is the mean of the darkest 5 % of MCD43A3 retrievals | Reported (design) | possibly biased low; moot while #1 stands |
| 20 | GGA `utilities_project.jl:207-211` | `project_decyear_bins` groups altimetry up to 9 days after its label; GEMB sampled at `date_center` | Reported | 9–15 % of seasonal amplitude left as residual |

## Latent or diagnostic-only

- GEMB `calculate_melt.jl:391`: `melt_total = max(0, melt − rain)` and rain frozen on contact is
  never counted as refreeze. Runoff and SMB unaffected; melt/refreeze outputs low (refreeze ≥ 3–11 %).
- GEMB `calculate_melt.jl:486`: `flux_dn` omits `freeze2` (unreachable today).
- GEMB `:BrunLefebre` + subsurface SW breaks the SW budget; `:GardnerSharp` + subsurface ignores
  diffuse albedo; GS2010 two-layer term likely uses kg m⁻² where m w.e. is meant (never activates).
- GEMB explicit solver: sub-steps comparable to the surface relaxation time; turbulent fluxes
  10–40 % off in windy transients; an all-low-density column diverges (not reached in production).
- GEMB energy checks rebuild supplied energy from the solver's own fluxes, and are off in
  production; they cannot detect flux-formula errors.
- GlacierSims `append_glacier_tile_netcdf`: cumulative series restart at zero at the seam;
  forcing means recomputed from the appended segment only. No file appended yet.
- GlacierSims: output files record no package versions or commits, and
  `run_parameter_differences` compares only shared keys, so "same experiment" checks cannot see a
  physics change (a `native_runs_lwfix` file passes against the old-physics reference).
- GlacierSims: a failed forcing load drops a cell with a `@warn`, not counted in
  `unrecovered_area`; native runs mark it complete and never retry. A band failing validation in any
  perturbation is dropped from all (only N26_E086, 3.4 km², in merged files).
- GlacierSims: `unreachable::BitVector` written under `Threads.@threads` (race; production runs `-t 1`).
- GlacierSims: `dh` includes the first interval while `dh_mass`/`dh_firn`/`dh_water` start at 0;
  `dv_rate` overstated by N/(N−1); merged files keep global attributes of the first file only.
- GlacierSims band elevations labelled `height_above_reference_ellipsoid` but are EGM2008
  orthometric (datums match each other; label is wrong). `precipitation_mean` written without ×8760
  (overwritten downstream).
- ClimateForcing: longitudes in [0, 360] select the wrong cell; NaN passes
  `validate_climate_forcing_units`; forcing cache frozen at its first `.zmetadata` read (record ends
  2026-08-14, ERA5T months never replaced); `time_range` silently truncated.
- The native pipeline in GlacierSims (`forcing_corrections.jl`, `glacier_native_run.jl`,
  `scripts/gemb_native_run.jl`) is untracked; all native outputs came from uncommitted code.

## Checked and found correct

ERA5-Land ARCO fluxes are hourly de-accumulated (verified on cached data), so ÷3600 and ×1000 are
right; time axis contiguous; storage order; wind, heights, units; geopotential height (List 1968);
vertical datums consistent; MODIS scale/fill/QA. Lapse sign, hypsometric pressure, Buck and
Konzelmann forms, Δz = 0 identity for sub-saturated input; calibration ΔT applied after decoupling.
GEMB grid merge/split/trim (incl. #82) conserve mass, energy and water; whole-run mass budget closes
to 2×10⁻¹⁰ kg m⁻²; densification units; melt and refreeze formulas; enthalpy-form diffusion;
Ψ functions; GS2010 terms; SSA and grain conversions; Brun wet growth (port fixes a MATLAB precedence
bug); surface SW conserves energy. GlacierSims unit conversions, area weighting, ΔT propagation to e
and LW, tile/native settings parity, time coordinate. GGA reading of GEMB outputs, firn term in the
calibration target, δice = 917.

## Tests that encode current (incorrect) behavior

- GlacierSims "elevation interval forcing" (`runtests.jl:1239-1243`): expected value uses the
  implementation's own lapse-then-decouple formula at Δz = 50 m (finding 2).
- GEMB `test_calculate_melt.jl:162` "Rain accounting" pins `melt − rain`.
- GEMB `test_calculate_temperature.jl` steady-state tests (both solvers) balance with
  LW↓ = εσT⁴, and the implicit "energy budget closes" uses full LW↓ (known longwave issue).
- GEMB "Stability functions satisfy their definition" checks formulas retyped in the test file
  and never calls GEMB code.

## Analytical tests to write

"✗" marks tests that fail on the current code. Expected values are independent of the
implementation: closed forms, hand calculations, conservation laws or limiting cases.

### GEMB.jl: thermodynamics and turbulent fluxes

1. ✗ Neutral sensible coefficient: ΔT = 1e-4 K, p = 60 kPa;
   `shf/ΔT = ρa·c_p·κ²U/(ln(z_w/z0)·ln(z_T/z_T0))`, rtol 1e-3 (fails by the Exner factor).
2. Neutral latent flux: Ta = Ts; `LE = κ²U·L·(e_a − e_s(Ts))/(R_v·Ts·ln·ln)`, L = L_s below 0 °C, rtol 1e-4.
3. ✗ Pressure invariance at neutral: `shf/ρa` equal at 100 and 60 kPa, rtol 1e-6.
4. ✗ Stable continuity: scan ΔT 0–15 K by 0.01 K at U = 2–8 m/s; adjacent |SHF| differ < 5 %.
5. ✗ Monin–Obukhov consistency against an iterated Obukhov-length solution with Beljaars–Holtslag
   Ψ at z_w/L and z_T/L, ζ < 1, rtol 5 %.
6. Semi-infinite solid, constant flux F = 100 W m⁻² into cell 1 (turbulence and LW off, ρ = 917,
   constant K, dz = 0.01 m, 2 m deep):
   `T(z,t) = T0 + (2F/k)√(κt/π)·exp(−z²/4κt) − (Fz/k)·erfc(z/2√(κt))`, 2 % of rise at 6 h, both solvers.
7. Discrete linear steady state: uniform K, dz, gradient F/K, F into cell 1; one step changes T by
   < 1e-8 K and `heat_flux_basal = −F` to 1e-8.
8. ✗ Radiative equilibrium: isothermal column at T_eq, LW↓ = σT_eq⁴, no SW, negligible turbulence;
   |ΔT| < 1e-3 K after one day (known longwave issue).
9. ✗ Independent energy budget: ΔH = (ΣSW + ε·LW↓ − εσTs⁴ + H + LE + G)·dt with H, LE recomputed
   outside the solver, to `energy_tolerance`.
10. Melting-point clamp: single active cell at 273.15 K, net flux Q; `M₁·c·(T₁ − 273.15) = Q·dt`, rtol 1e-9.
11. ✗ Time convergence: windy cold case, safety factor 0.8 vs 0.05 flux difference < 2 %.
12. ✗ Wet-cell latent buffer: wet surface cell under 100 W m⁻² cooling holds 273.15 K until its
    water has refrozen, 1e-6 K.
13. ✗ Vapor phase energy: latent energy = EC·L_phase with L_s for mass exchanged with ice, rtol
    1e-12; at a melting surface, condensing c kg raises column enthalpy by c·(h_ice(0 °C) + L_f).
14. Beer–Lambert SW: uniform B, cell i absorbs `swf_ss·(e^{−Bz_{i−1}} − e^{−Bz_i})` plus 0.36 of net
    at the surface, rtol 1e-12; absorbed + transmitted = total.
15. Saturation vapor pressure: 610.78 Pa at 273.16 K exactly; within 0.5 % of Murphy & Koop (2005),
    233–273 K.
16. Spin-up skin temperature: U at floor, no SW: T_skin = (LW↓/σ)^{1/4}, 1e-3 K.
17. Emissivity melt switch: at melting point, next sub-step upward LW = ε_large·σ·273.15⁴, 1e-12.

### GEMB.jl: albedo, shortwave, grains

18. GS2010 closed forms, deep snow r_e = 0.5 mm, c = 0: SZA 0, τ 0 → α = 1.48 − S^−0.07 = 0.73410;
    SZA 60 adds 0.53·α_s(1−α_s)·0.5^1.2; τ = 5 adds 0.1τ·α_s^1.3/(1+1.5τ)^α_s with
    u = 0.64x + (1−x)cosθ; c = 1 ppmw per GS2010 dα_c; rtol 1e-12.
19. Limits: α(SZA = 0) = α_s exactly; α non-decreasing in SZA and τ.
20. ✗ Bare ice: all-ice column, `albedo_ice = 0.30`, M = 2 kg m⁻² → α = 0.30 ± 0.01; α independent
    of dt for the same state.
21. ✗ Bare-ice integration: 20-day bare-ice runs at `albedo_ice` 0.30 vs 0.48; absorbed-SW ratio
    0.70/0.52 = 1.346 ± 2 %.
22. ✗ Thin snow over ice: α → `albedo_ice` as snow depth → 0; semi-infinite value when optically
    thick; monotonic in depth.
23. ✗ (pending GS2010 units) Two-layer A: r_e = 1 mm, z = 3 mm w.e.: A = 2.1·0.003^(1.35(1−α_s)−0.13).
24. ✗ SW conservation for every albedo method, subsurface on and off, diffuse > 0:
    Σswf = (1−α)(SW − SW_d) + (1−α_d)·SW_d, rtol 1e-10.
25. Brun spectral albedo: α = Σ sF_i·a_i(gsz); Σ sF_i = 1.
26. Density-branch extinction, uniform ρ = 300 (B = 10): cell i absorbs
    0.64(1−α)SW·(e^{−BΣ<i} − e^{−BΣ≤i}), plus 0.36(1−α)SW in cell 1; rtol 1e-12.
27. Grain growth closed forms, rtol 1e-12: dendritic dry low gradient (Δδ = −2e8·e^(−6000/T)·dt_d,
    Δs = +1e9·e^(−6000/T)·dt_d, and the diameter formula); dendritic dry high gradient
    (Δδ = Δs = −2e8·e^(−6000/T)|∇T|^0.4·dt_d); dendritic wet (Δδ = −θ³/16·dt_d, θ = 100W/(ρdz) ≤ 9);
    Brun 1989 wet (d_new = (d³ + (6/π)(1.28e-8 + 4.22e-10θ³)dt_s)^(1/3)); Marbouty at −10 °C,
    ρ = 200, ∇T = 30 K/m (Q = 0.0167 mm d⁻¹).
28. Temperature gradient: linear T(z) on a non-uniform grid; `_grain_gradient` returns g exactly at
    every cell including edges.
29. Grain caps: s = 1 caps at 2 mm diameter, s < 1 at 5 mm, exactly.

### GEMB.jl: mass and column

30. ✗ Rain does not compact snow: 0.05 m, 350 kg m⁻³ cell at 0 °C + 5 kg rain at 275.15 K through
    accumulation and melt: ρ = 350, dz = 0.05 (1e-9); water = 5 + 5·c_w·2/L_f = 5.1252 kg (1e-6).
31. ✗ Rain refreeze: 2 kg rain at 274.15 K on 17.5 kg of −10 °C snow: refreeze =
    [c·M·10 − P·c_w·1]/L_f = 1.0747 kg (1e-9), melt = 0; over any run
    Σ(melt + rain − refreeze − runoff) = ΔW (1e-6 kg m⁻²).
32. ✗ Trace snowfall: 0.01 kg on 17.5 kg with r_e = 1 mm: r_e ≥ 1 − (0.01/17.51)·0.95 mm, or the
    docstring states an optical-depth rule that the test then checks.
33. Melt from excess: cell M at 0 °C + ΔT: melt = c·M·ΔT/L_f (rtol 1e-12); ΔT = L_f/c + δ melts the
    cell fully and the next cell gains c·M·δ.
34. Percolation into cold firn: refreeze = min(W − W_irr, cold content/L_f, (ρ_i − ρ)·dz), runoff 0.
35. Runoff on ice: melt on a 917 kg m⁻³ cell gives runoff = melt exactly.
36. Densification one step: ρ₁ = ρ₀ + c·(ρ_i − ρ₀)·dt/365 with Arthern
    c₀ = 0.07·A·g·exp(−60000/RT + 42400/RT̄), rtol 1e-12.
37. Herron–Langway steady state (1980 eq. 7): ln(ρ/(ρ_i − ρ)) linear in depth with
    k₀ = 11·exp(−10160/RT) (stage 1) and k₁√A (stage 2); `steady_state_density` within 0.5 kg m⁻³.
38. Elevation identity: −Σ`ice_flux` = (SMB − ΔW)/ρ_i + ΔFAC, 1e-9 m.
39. Spin-up: constant-climate dry column reaches the Herron–Langway stage-1 profile within 2 kg m⁻³.
40. `flux_dn` mass: a cell with pore water above irreducible and cold content conserves mass,
    flux_dn = in − F1 − dW − F2.

### GEMB_ClimateForcing: ingest

41. Conversions on a local ARCO-shaped fixture (or factored pure function): ssrd = 3600·[0, 100,
    350.5], strd = 3600·300, tp = [0, 1e-3, 2.5e-3] m, (u, v) = (3, 4), (−5, 12) → SW = [0, 100,
    350.5], LW = 300, P = [0, 1, 2.5] kg m⁻², ws = [5, 13]; rtol 1e-12.
42. Accumulation-convention guard (cached data, 46.5° N 8° E, January): ssrd ≤ 1 J m⁻² 20–04 UTC;
    Σtp over 2023 in 0.5–4 m.
43. Interval assignment: SW centroid − (12 − λ/15) = +0.5 ± 0.25 h over a year.
44. ✗ Longitude: on −179.9:0.1:180, lon = 210 selects −150.0; lon = −179.97 selects 180.0.
45. ✗ Dewpoint against IFS: Td = −40, −20, 0, 20 °C → 18.74, 125.27, 610.77, 2335.84 Pa, rtol 2e-3.
46. ✗ RH self-consistency: unclamped e(Td = T)/e_s(T) = 1 ± 1e-3 for −40 to 30 °C.
47. Geopotential height: at the latitude where List g = g₀ (≈45.54°), Φ = 1e4·g₀ → Z = 10015.8 ± 0.5 m; Φ = 0 → 0.
48. Geoid: EGM96 N(0, 0) = 17.16 ± 0.3 m (second reference point to be checked against GeographicLib first).
49. Invariant seam: `surface_elevation` at lon −0.1 equals `z[X = 359.9]/g` exactly.
50. ✗ NaN validation: one NaN per layer raises `ArgumentError` mentioning "non-finite".
51. Time axis: `diff(time)` all 1 h; 2020 request returns 8784 steps.

### GEMB_ClimateForcing: downscaling

52. ✗ End-to-end on-glacier temperature: 20 cells z = 500–2500 m, glm = clamp((z−500)/2000)² + noise,
    k = 0.8, E_amb = 20 − 0.006z; derive k, Γ; band at 3000 m from any donor gives
    T = 273.15 + k(20 − 0.006·3000) = 274.75 K, 1e-6 K.
53. ✗ Truncated-excess recovery: as 52 with E_amb = 14 − 0.006z: k(z̄) = 0.80 ± 0.01, Γ = 4.80 ± 0.05 K/km.
54. ✗ Identity on saturated input: T = 253.15 K, e = `dewpoint_to_vapor_pressure(253.15)`;
    `climate_adjust_for_elevation(s, 0)`, `temperature_adjust(s, 0)` and
    `_cell_forcing_at_interval(s, 0, k = 1)` return e and LW exactly.
55. ✗ ΔT = 0 after decoupling: T = 285.15 K, RH 0.95, k = 0.6; `temperature_adjust(·, 0)` returns its input exactly.
56. ✗ k_eff weighting: glm = 0.5 donor at Δz = 0 above melting: decoupled excess = k·E_amb, 1e-9.
57. Pressure, descending: Δz = −800 m, Γ = 5 K/km: P′ = P·exp(g·800/(R_d·(T + T′)/2)), T′ = T + 4, rtol 1e-12.
58. Decoupling LW: LW′ = (ε_cs(e, T′) + LW/(σT⁴) − ε_cs(e, T))·σT′⁴, rtol 1e-12.
59. Pooled albedo: 400 retrievals 0.30:0.00075:0.59925 → mean of darkest 20 = 0.307125, 1e-6;
    values below 0.25 or rejected QA excluded from the count.
60. Decoupling-factor hold: E_amb crosses 0.5 K at z*; k(z > z*) = k(z*) exactly; `n_held` counts held steps.

### GEMB_GlacierSims

61. Perturbation chain: T = 263.15 K, RH 0.8, LW = 0.8σT⁴, P = 1; ΔT = +2, pscale 1.5 through the
    driver's chain: T′ = T + 2, e′ = 0.8·e_s(T′), LW′ = (ε_cs(e′, T′) + Δε)σT′⁴, P′ = 1.5, all else
    bit-identical; ✗ RH > 1 at ΔT = 0 returns its input.
62. `_tile_totals` by hand: areas [10, 30] km², runoff [1,2,3], [4,5,6] kg m⁻², dh [−0.1,−0.2,−0.3],
    [0,−0.1,−0.2] m → `total_runoff` = [1.3, 3.0, 5.1]e-4 Gt, `dv` = [−0.001, −0.005, −0.009] km³.
63. ✗ Append continuity: dh [−1,−2] then appended [−1,−2] → [−1,−2,−3,−4]; totals continuous.
64. ✗ Components close: `dh = dh_mass + dh_water + dh_firn + strain`, 1e-12.
65. ✗ Bare-ice albedo reaches GEMB: two bands at `albedo_ice` 0.3 and 0.5 on identical
    SW-dominated warm forcing; melt ratio within 10 % of 0.7/0.5.
66. Hypsometry: fully glaciated 0.1° cell at 1050 m → bin area R²·Δλ·(sin φ₂ − sin φ₁), rtol 1e-3;
    two adjacent cells sum to the 0.2° box; a −10 m pixel is counted or throws.
67. Time integration: P = 1 kg m⁻² h⁻¹, weekly output, three weeks → 168 kg m⁻² per interval,
    timestamps at day 7, 14, 21 23:00, `total_precipitation` = 168·k·area·1e-6 Gt.
68. ✗ Summary rate: dv rising r per week → rate r·365.25/7 per year.
69. ✗ Threaded `unreachable` flags: 64 failing intervals, 8 threads, 100 repetitions, all flagged.

### GlobalGlacierAnalysis

70. ✗ Discharge units: a geotile with known discharge D Gt yr⁻¹ gives `discharge` column
    D·1000/δice km³ yr⁻¹ (not D·(1000/δice)²).
71. ✗ Antarctic `discharge2smb`: a synthetic SMB trend of s m w.e. yr⁻¹ over area A returns
    s·A/1000 Gt yr⁻¹.

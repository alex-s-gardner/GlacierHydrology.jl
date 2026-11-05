# Known issue: GEMB absorbs all downward longwave, not emissivity times it

**Status: fixed on GEMB.jl branch `fix/physics-audit` (commit `0d1fb7d`, local, not merged); no run
uses the fix yet.** Every existing GEMB result (the tile sweep, the `gemb_calibration` corrections, the
native-resolution run) was produced with the behavior described here. On the branch the upward
longwave is emitted plus reflected, `εσTs⁴ + (1−ε)·LW↓`, with ε taken per sub-step in both solvers, so
the surface absorbs `ε·LW↓` and every net-longwave diagnostic (`LW↓ − LW↑`) agrees with the solver.
The fix is one of the changes listed in `gemb_pipeline_review_2026-10-13.md`; all GEMB outputs change
and need the sweep, calibration and native run redone.

## What the code does

In GEMB.jl 2.0.0 (the `main` checkout the project manifests pin), the thermal solve adds the whole
downward longwave flux to the top layer and scales only the emitted flux by the surface emissivity
(`mp.emissivity`, default 0.97, `types.jl`):

    Q_surface = LW_down·dt − ε·σ·T⁴·dt + turbulent fluxes

The standard surface energy balance is `(1−α)·SW + ε·LW_down − ε·σ·T⁴ + H + LE`. A gray surface
absorbs `ε·LW_down` and reflects `(1−ε)·LW_down`, so the model absorbs `(1−ε)·LW_down` too much. With
ε = 0.97 and downward longwave of 250–300 W m⁻², that is 7–9 W m⁻² of extra surface heating.

Where, in `src/` of GEMB.jl:

- `calculate_temperature.jl`, explicit solver (the default): lines 241 (`longwave_downward =
  cfs.longwave_downward * dt`), 288 (`Q_surface`) and 351 (the verbose energy check).
- `calculate_temperature.jl`, implicit solver: lines 617, 693, 762 and 805.
- `gemb_core.jl:199`: the net-longwave term of the verbose column energy check uses the same full
  `LW_down`, so the checks agree with the solver and do not detect the problem.
- `initialize_climate_summary.jl:419` (the spin-up initial surface energy balance) does use
  `ε * lw_in`, so the initializer and the solver disagree. The comment at line 93 of that file writes
  the balance without `ε` on `LW_down`, which matches neither the code beside it nor the solver.

No comment or documentation in the package gives a reason for full absorption. Whether this matches
the original MATLAB GEMB has not been checked.

## Evidence: size of the effect on one tile

Tile S50_E068 (Kerguelen, 24 native cells, 531 km²), native-resolution run at the corrections the
production native run used (pscale 1.25, ΔT −2.0 K), 1950–2026 forcing, rates over 2000-01-01 to
2024-12-31. The only change from the production run is that the absorbed longwave is `ε·LW_down`,
using the emissivity in force at each sub-step (the emissivity melt switch can change it inside the
loop). Production output is the baseline.

| | baseline | absorbed `ε·LW_down` | change |
|---|---|---|---|
| mass change, Gt yr⁻¹ | +0.158 | +0.314 | +0.156 (0.29 m w.e. yr⁻¹) |
| melt, Gt yr⁻¹ | 1.271 | 1.142 | −10 % |
| runoff, Gt yr⁻¹ | 1.565 | 1.419 | −9 % |
| refreeze, Gt yr⁻¹ | 0.378 | 0.391 | +3 % |
| rain, precipitation | unchanged | | |

All 24 cells lose less mass (per-cell change in mass-change rate: median +0.005, range +0.0002 to
+0.025 Gt yr⁻¹). The energy-conservation check passed in every cell.

The tile's reference grid (`tile_runs_merged`) gives d(mass change)/d(ΔT) = −0.69 Gt yr⁻¹ K⁻¹ at
pscale 1.25 between ΔT −3 and −1, so the fix is equivalent to an air-temperature offset of −0.23 K.
Refitting to the same altimetry would need ΔT about 0.23 K higher.

## What it affects, and what was not tested

Calibration pins the mass change to altimetry, so the main consequence is a bias in the fitted
`ΔT` (and, through the valley in the cost surface, `pscale`) of roughly +0.2 K for a tile like this
one, plus a second-order change in how the fitted mass budget splits into melt, runoff and refreeze.
This is small against the spread of fitted ΔT across groups (prior standard deviation 1.65 K).

Not tested:

- Any other climate regime. S50_E068 is maritime with high downward longwave; the effect scales with
  `LW_down` and with how much of the extra energy goes to melt rather than being re-emitted.
- The seasonal cycle, and runoff and melt after a refit at the corrected physics.
- The implicit solver (patched, not run; the explicit solver is the default and what every run used).

## Reproducing the test

Copy the installed GEMB.jl, apply the diff below to `src/calculate_temperature.jl` in the copy, and
make the project under test `Pkg.develop` the copy (use a scratch copy of the project, not the
real one). Run `scripts/gemb_native_run.jl` of GEMB_GlacierSims.jl for tile S50_E068 with
`FORCING_CORRECTIONS` set to the reference member's `_gembfit.arrow` that was in force when the
production native run was made, and compare cell totals with the production run. The test output is
in `/mnt/bylot-r3/data/era5land/native_runs_lwfix/S50_E068/`.

```diff
-        Q_surface = (longwave_downward + longwave_upward) + thf
+        Q_surface = (emissivity * longwave_downward + longwave_upward) + thf
-            E_supplied = sw_total + longwave_downward + longwave_upward + thf + heat_flux_basal
+            E_supplied = sw_total + emissivity * longwave_downward + longwave_upward + thf + heat_flux_basal
-                Q_surface = longwave_downward +
+                Q_surface = emissivity * longwave_downward +
-        Q_surface = (longwave_downward + longwave_upward_step) + thf
+        Q_surface = (emissivity * longwave_downward + longwave_upward_step) + thf
-            E_supplied = sw_total + longwave_downward + longwave_upward_step + thf +
+            E_supplied = sw_total + emissivity * longwave_downward + longwave_upward_step + thf +
```

`gemb_core.jl:199` also needs the absorbed term in the same form for its verbose check (not part of
the test above, which ran with `verbose = false`).

## To fix

1. Apply the change above in GEMB.jl (upstream), and decide whether `initialize_climate_summary.jl`
   keeps `ε·LW_down` (it already does) and its line 93 comment is corrected.
2. Rerun the GEMB tile sweep, `gemb_calibration` and the native-resolution run; nothing downstream of
   them can be corrected after the fact, because the fitted corrections and the melt/runoff split both
   move.
3. Record the change in the methods record
   (`notes/methods_2026-10_seasonal_cycle_and_calibration.md`) with the before/after calibration.

## For the Extended Methods

State that GEMB's surface absorbs the full downward longwave flux (no emissivity factor on absorption),
which biases the energy balance warm by `(1−ε)·LW_down`, and that this is absorbed into the fitted
temperature offset (about +0.2 K in the one tile tested).

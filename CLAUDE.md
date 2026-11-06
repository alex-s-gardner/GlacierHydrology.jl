# CLAUDE.md

Guidance for Claude Code in this repository.

## What this package is

GlacierHydrology.jl turns monthly glacier volume change into glacier mass and energy fluxes and their
contribution to rivers and to the people who depend on them. It is the code of the hydrology paper
(`glacier-hydrology-paper`):

1. It reads the volume change of every ensemble member from the products of
   [GlacierVolumeChange.jl](https://github.com/alex-s-gardner/GlacierVolumeChange.jl) (GVC), built by
   `GVC.create_dataset`.
2. It reads a sweep of GEMB runs over precipitation scaling (`pscale`) and temperature offset (`ΔT`), made
   by GEMB_GlacierSims.jl as tile NetCDF.
3. It calibrates GEMB to the altimetric seasonal and long-term volume change for each group of geotiles,
   which gives accumulation, melt, refreeze, runoff and firn air content.
4. It adds solid-ice discharge, routes glacier and land-surface runoff through MERIT Hydro, and computes
   gmax (the largest glacier share of river flux) and the population along glacier-fed rivers.

The altimetry chain up to volume is in GVC, and this package never recomputes `dv`. The snow, firn and
surface energy-balance physics are in GEMB.jl, and the forcing and the tile sweeps are in
GEMB_ClimateForcing.jl and GEMB_GlacierSims.jl. `PARTITION.md` in the paper repositories says which paper
describes which part.

Scripts load the package with `import GlacierHydrology` and qualify every name
(`GlacierHydrology.gemb_calibration`); there is no short alias. GVC is `import GlacierVolumeChange as GVC`.

## Running

- **Hosts:** the package loads only on bylot, devon and baffin. `setpaths()` here and in GVC runs at
  module load and throws on any other host. That includes the tests.
- **Julia:** 1.12 (`julia = "1.12"`). On the servers use `~/.juliaup/bin/julia +1.12`.
- **Workflow:** `src/run_all.jl` runs steps 11–24 after GVC has built the volume dataset (steps 1–10).
  Each step is an include or a function call, so a session can run any one of them after the setup block.
- **Tests:** `julia --project test/runtests.jl` on a JPL server. The tests use synthetic data and mock
  files and read nothing under `/mnt`.

## Steps of `run_all.jl`

| Step | What | Where |
|---|---|---|
| 11 | GEMB tile sweep → geotile ensemble, per surface mask | `gemb_tiles_binning.jl`, `utilities_gemb_tiles.jl` |
| 12 | Solid-ice discharge (flux gates; Antarctic periphery from GEMB SMB) | `global_discharge_filled` |
| 13 | Calibration of (`pscale`, `ΔT`) per geotile group | `gemb_calibration`, `gemb_bestfit_grouped` |
| 14 | Calibrated synthesis for every ensemble member | `geotile_synthesis_gembfit_dv` |
| 15–16 | Glacier-level summary; geotile rates as GIS files | `glacier_summary_file`, `gembfit_dv2gpkg` |
| 17–18 | Glacier and land-surface runoff routing | `glacier_routing`, `land_surface_model_routing.jl` |
| 19–20 | gmax and its figures | `gmax_global.jl`, `gmax_point_figure.jl` |
| 21–23 | Regional results, Extended Data figures, river-buffer population | the scripts of those names |
| 24 | Gridded NetCDF outputs | `gridded_output_files.jl` |

## Key concepts

- **Geotiles:** 2° cells named `lat[+60+62]lon[-142-140]`. A geotile *group* joins the geotiles a glacier
  larger than 100 km² crosses, so ice flow across a boundary does not bias the fit.
- **Ensemble:** `ensemble_paths()` gives the path of every altimetry member in each form the pipeline
  uses (filled, synthesized, calibrated), from `GVC.ensemble_members` and `GVC.reference_member`. Errors
  are ensemble spreads about the reference member (`error_quantile = 0.95`, `error_scaling = 1.5`).
- **GEMB runs and physics members:** `gemb_info(; gemb_run_id, physics_member)`.
  - Run 9 is the production sweep (`$CLIMATE_CACHE/tile_runs_v2/<member>/`).
  - Run 10 is the pilot (`tile_runs_pilot`), on `gemb_pilot_forcing_grid`.
  - A physics member such as `S-kfit-table` names the grain-growth, decoupling and impurity levels of the
    sweep, with optional wind (`w<factor>`, `ws`) and precipitation-gradient (`g<g>`, applied here)
    tokens (`gemb_physics_member_parts`). Production uses GEMB_GlacierSims' default member,
    `gemb_reference_physics_member = "S-kfit-table"`: GEMB's default snow set, ERA5-Land wind and no
    gradient.
  - Calibrated outputs of a member go to `gemb_physics/<member>/` beside the altimetry files
    (`gemb_physics_path`).
- **Calibration cost:** the RMS of the residual blended with its seasonal amplitude
  (`seasonality_weight = 0.85`), times `1 + wd·d`. Here `d` is the Mahalanobis distance from the empirical
  forcing prior `gemb_forcing_prior` and `wd = distance_from_origin_penalty = 0.35`.
  - The prior and `wd` were derived from GEMB run 8 (`gemb_forcing_prior_run_id`).
  - `gemb_calibration` throws when they would be applied to another run. Re-derive all three from that
    run's unpenalized fits (`distance_from_origin_penalty = 0`) and update them together.
  - The method record is `notes/methods_2026-10_seasonal_cycle_and_calibration.md`.
- **Dates:** GEMB is binned on the 30-day bins of `GVC.project_date_bins`, extended back to 1970, as the
  altimetry is. `gemb_calibration` throws if the two date axes differ over their common span.
- **Freshness:** steps 12–14 rebuild any file older than the GEMB ensemble they descend from
  (`gemb_ensemble_built`). `file_is_current` compares modification times in UTC.
- **Units:** GEMB volumes are km³ of ice equivalent. Masses are Gt, using `δice = 917` kg m⁻³ (the
  density the GEMB runs record), and `Unitful.register(GlacierHydrology.MyUnits)` is needed for `Gt`.

## Layout

| Path | Contents |
|---|---|
| `src/GlacierHydrology.jl` | module; the GVC names it imports |
| `src/local_paths.jl` | data paths used only here, merged with `GVC.setpaths()` |
| `src/utilities_project.jl` | calibration constants, `gemb_info`, `ensemble_paths`, cost function |
| `src/utilities_gemb.jl`, `utilities_gemb_tiles.jl`, `utilities_gemb_physics.jl` | GEMB ensemble, calibration, physics members |
| `src/utilities_synthesis.jl` | geotile groups, discharge, calibrated synthesis |
| `src/utilities_routing.jl` | river routing and gmax |
| `src/utilities_postprocessing.jl` | regional sums, trends, summary and GIS files |
| `src/utilities_wgms.jl` | WGMS seasonal-balance scoring (an independent check, not a calibration target) |
| `src/utilities_readers.jl`, `utilities_plotting.jl`, `utilities.jl` | readers for comparison data, figures, NetCDF I/O |
| `src/utilities_response2reviewers.jl` | WGMS helpers used by the Extended Data script (not part of the module) |
| `scripts/audit/` | discharge-unit and surface-mask audit |
| `docs/data_archive/` | readme and file list of the data archive |
| `notes/` | method records of the calibration and of the GEMB pipeline review |

## Working rules

- Fail fast: an unexpected input throws with a message that names the file or setting at fault. Do not
  add silent fallbacks.
- Keep geotile groups independent of each other in the calibration.
- A change to a GVC name in the `using GlacierVolumeChange:` list, to `setpaths()`, or to the volume
  product format is a change to both packages.
- GEMB_GlacierSims.jl reads this package's calibration output (`…_synthesized_gembfit.arrow` under
  `gemb_physics/<member>/`) for its native-resolution runs. Keep that name and location, or change both
  packages together. Its `gga_reference_ensemble_file` looks for `const reference_ensemble_file` in this
  package's `src/utilities_project.jl`, which is now `GVC.reference_ensemble_file`. Until that lookup
  moves to GVC, native runs need `REFERENCE_ENSEMBLE_FILE` or `FORCING_CORRECTIONS` set.

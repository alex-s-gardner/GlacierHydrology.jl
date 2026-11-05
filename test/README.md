# GlacierHydrology.jl tests

```bash
julia --project test/runtests.jl
```

The tests build synthetic inputs with known answers and mock the external files they read, so they
read nothing under `/mnt`. The package itself loads only on bylot, devon and baffin, so the suite runs
there.

| Directory | Contents |
|---|---|
| `unit/` | array helpers; `gemb_info`, `file_is_current` and the calibration cost and grid search |
| `algorithms/` | GEMB ensemble and calibration (`test_gemb.jl`), tile reading (`test_gemb_tiles.jl`), physics members (`test_gemb_physics.jl`), discharge (`test_discharge.jl`), regional units (`test_regional_units.jl`), post-processing, routing, WGMS scoring |
| `io/` | DimArray, JLD2 and NetCDF round-trips; readers for discharge, GRACE and GlaMBIE |
| `integration/` | regional aggregation and statistical invariants |
| `fixtures/` | generators for synthetic geotiles, time series, river networks and mock data files |

New tests use synthetic data with a known answer. They check error messages with
`@test_throws "message" expr`.

# GlacierHydrology.jl

Glacier mass and energy fluxes from satellite altimetry and a calibrated surface mass-balance model, and
the glacier contribution to rivers and to the populations that depend on them.

The package takes monthly glacier volume change for every ensemble member from
[GlacierVolumeChange.jl](https://github.com/alex-s-gardner/GlacierVolumeChange.jl). It calibrates a sweep
of GEMB runs (made by GEMB_GlacierSims.jl) to the seasonal and long-term volume change of each group of
2° geotiles, which gives accumulation, melt, refreeze, runoff and firn air content. It then adds
solid-ice discharge, routes glacier and land-surface runoff through MERIT Hydro, and computes gmax (the
largest glacier share of river flux) and the population along glacier-fed rivers.

It is the code of Gardner et al., *Glacier contributions to rivers and oceans in the early twenty-first
century* (in review).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/alex-s-gardner/GlacierHydrology.jl")
```

The package is set up for the JPL servers bylot, devon and baffin. It reads data from `/mnt/bylot-r3`
and `/mnt/devon-r2`, and it throws at load on any other host. To run elsewhere, edit `src/local_paths.jl`
here and `src/config/paths.jl` in GlacierVolumeChange.jl.

## Running the workflow

Build the volume dataset first with `GlacierVolumeChange.create_dataset` (steps 1–10). Then:

```bash
julia --project src/run_all.jl
```

| Step | What |
|---|---|
| 11 | Bin the GEMB tile sweep onto geotiles, per glacier inventory |
| 12 | Solid-ice discharge, with the Antarctic periphery from modeled SMB |
| 13 | Calibrate GEMB precipitation scaling and temperature offset for each geotile group |
| 14 | Calibrated synthesis for every ensemble member |
| 15 | Glacier-level summary file |
| 16 | Geotile rates as GIS files |
| 17 | Route glacier runoff through river networks |
| 18 | Route land-surface model runoff through the glacier-fed reaches |
| 19–20 | gmax and its figures |
| 21 | Regional results |
| 22 | Extended Data figures |
| 23 | Population within buffers of glacier-fed rivers |
| 24 | Gridded NetCDF outputs |

Each step skips outputs that are newer than their inputs. Steps 12–14 rebuild anything older than the
GEMB ensemble they descend from.

## Tests

```bash
julia --project test/runtests.jl
```

The tests use synthetic data and mock files. They read nothing under `/mnt`, but they have to run on a
JPL server because the package loads only there.

## Data products

- Regional time series and trends of glacier mass, volume and flux components (NetCDF, CSV)
- Geotile rates (GeoPackage)
- River reaches with gmax and glacier runoff (GeoPackage)
- Population near glacier-fed rivers (NetCDF)
- Gridded monthly fluxes at 0.5° and 2° (NetCDF)

`docs/data_archive/` describes the files of the data archive.

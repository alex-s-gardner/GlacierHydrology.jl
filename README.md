# GlobalGlacierAnalysis.jl

A Julia package for comprehensive glacier elevation change and hydrological analysis using multi-mission satellite altimetry and model data.

## Overview

GlobalGlacierAnalysis.jl provides an end-to-end workflow for processing satellite altimetry, glacier elevation change, and model data to quantify glacier mass changes and their impacts on river systems and populations.

## Installation

```julia
using Pkg
Pkg.add("https://github.com/alex-s-gardner/GlobalGlacierAnalysis.jl")
```

## Quick Start

### 1. Data Requirements

**Note:** This package is designed for use on JPL servers (bylot, devon, baffin) and expects data at specific locations. If running elsewhere, adapt `src/local_paths.jl` (paths used only here) and `src/config/paths.jl` in [GlacierVolumeChange.jl](https://github.com/alex-s-gardner/GlacierVolumeChange.jl) (the altimetry chain and the paths both packages use). `local_paths.jl` merges the two tables.

**Required datasets:**
- Satellite altimetry data (GEDI, ICESat, ICESat-2)
- Glacier elevation change stacks (Hugonnet et al.)
- Global DEMs (REMA, COP30, ArcticDEM, NASADEM)
- Glacier outlines (RGI v6/v7)
- Canopy height data (ETH Global Canopy Height)
- River network and basin data (MERIT Hydro, BasinATLAS, GRDC)
- Validation datasets (GRACE, Zemp2019, etc.)

### 2. Run the Workflow

The glacier volume-change dataset is built first, by [GlacierVolumeChange.jl](https://github.com/alex-s-gardner/GlacierVolumeChange.jl) (`GlacierVolumeChange.create_dataset`): altimetry archives, ancillary data, binning, filling, synthesis, and one volume product (NetCDF) per ensemble member. The workflow here is orchestrated by `src/run_all.jl` and reads those products:

```bash
julia --project src/run_all.jl
```

It executes the following sequence (numbered as in `run_all.jl`; steps 1–10 are the volume dataset):

11. **Create glacier model (GEMB) classes** - Generates glacier energy and mass balance classifications
12. **Calculate global glacier discharge** - Computes ice discharge with filled Antarctic data using modeled SMB
13. **Calibrate GEMB model to altimetry data** - Optimizes glacier model parameters for grouped geotiles
14. **Calibrate GEMB model for each synthesized geotile dataset** - Applies model calibration to individual datasets
15. **Generate glacier-level summary files** - Creates key statistics and metrics for further analysis
16. **Export geotile-level glacier change trends** - Produces GIS-compatible files with trends and amplitudes
17. **Route land surface model runoff** - Processes terrestrial water fluxes through river networks
18. **Route glacier runoff** - Processes glacier meltwater through river networks
19. **Calculate gmax** - Determines maximum glacier contribution to river flux
20. **Generate point-based figures** - Creates gmax visualization plots
21. **Produce regional results** - Generates regional analysis outputs for sharing
22. **Generate extended data figures** - Produces manuscript-quality visualizations
23. **Analyze population affected by glacier-fed river changes** - Quantifies population impact using buffer analysis
24. **Create output files for the manuscript** - Generates submission-ready data and figures

Building the full workflow takes several days, depending on data and hardware.

### 3. Configuration

The workflow can be customized by modifying parameters at the top of `src/run_all.jl` and parameters within scripts called by `src/run_all.jl`


### 4. Checkpoint System

The workflow implements checkpoint logic to avoid redundant processing:
- Set `force_remake = false` for normal operation
- Use `force_remake_before` parameters for selective reprocessing
- Comment/uncomment specific steps in `run_all.jl` for partial runs

## Key Components

### Core Processing Scripts
- `run_all.jl` — Main workflow orchestrator
- `land_surface_model_routing.jl` — Terrestrial runoff routing
- `glacier_routing.jl` — Glacier runoff routing (called from `run_all.jl`)
- `gmax_global.jl` — Maximum glacier contribution analysis
- `gmax_point_figure.jl` — Visualization generation
- `regional_results.jl` — Regional analysis outputs
- `manuscript_extended_data_figures.jl` — Publication figures
- `river_buffer_population.jl` — Population impact analysis
- `create_output_files.jl` — Manuscript submission outputs

### Utility Modules

The altimetry chain (archives, ancillary extraction, binning, filling, synthesis), the height-to-volume conversion, its diagnostic plots and the archive maintenance scripts are in [GlacierVolumeChange.jl](https://github.com/alex-s-gardner/GlacierVolumeChange.jl). This package reads volume change from its products (`GlacierVolumeChange.read_volume(GlacierVolumeChange.volume_product_path(member))`) and imports a short list of its utilities (`using GlacierVolumeChange: …` in `src/GlobalGlacierAnalysis.jl`); scripts call the rest as `GVC.<name>` after `import GlacierVolumeChange as GVC`.

**This package:**
- `local_paths.jl` — Data locations (`setpaths`, `pathlocal`), merged with GlacierVolumeChange's
- `utilities.jl` — General utility functions and helpers
- `utilities_project.jl` — GEMB calibration constants and project configuration
- `utilities_readers.jl` — Readers for comparison datasets (GRACE, GlaMBIE, discharge, runoff)
- `utilities_gemb.jl`, `utilities_gemb_tiles.jl` — Glacier Energy and Mass Balance (GEMB) model utilities
- `utilities_synthesis.jl` — Geotile grouping, discharge, and GEMB fits to the synthesized altimetry
- `utilities_routing.jl` — River routing and hydrological analysis utilities
- `utilities_postprocessing.jl` — Post-processing and data export functions
- `utilities_plotting.jl` — Visualization and plotting utilities
- `utilities_response2reviewers.jl` — Response-to-reviewers analysis and outputs

## Dependencies

Core dependencies include:
- **Geospatial:** Proj, GeoArrays, SpaceAltimetry, Geodesy, Rasters, Shapefile, LibGEOS, GeoInterface
- **Data Processing:** DataFrames, Statistics, FileIO, Arrow, HTTP, MAT, CSV, NCDatasets, JLD2
- **Visualization:** CairoMakie, ColorSchemes
- **Scientific Computing:** NonlinearSolve, LsqFit, Distributions, LinearAlgebra
- **Custom Packages:** GlacierVolumeChange, RangeExtractor, GeoTiles, BinStatistics

Install all dependencies with:
```julia
Pkg.instantiate()
```

## Data Outputs

The workflow produces:
- **Geotile-level data** - Processed elevation change and statistics
- **Global datasets** - Synthesized glacier mass change estimates
- **River routing results** - Glacier contribution to river systems
- **Population analysis** - Affected population statistics
- **GIS outputs** - Geospatial files for further analysis
- **Visualization figures** - Publication-ready plots and maps

## License

MIT License

## Contributing

Contributions are welcome! Please submit issues and pull requests on GitHub.

---

**For questions or help, please open an issue on GitHub.**

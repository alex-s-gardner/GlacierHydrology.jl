"""
Test helper functions and utilities for GlacierHydrology.jl tests
"""

using Dates
using Statistics
using DataFrames
using DimensionalData
import DimensionalData as DD

"""
    synthetic_elevation_timeseries(; dates, trend=-0.5, amplitude=0.2, phase=0.0, noise_sigma=0.1)

Generate synthetic elevation time series with known trend and seasonality.

# Arguments
- `dates`: Vector of DateTime objects
- `trend`: Linear trend in m/yr (default: -0.5)
- `amplitude`: Seasonal amplitude in m (default: 0.2)
- `phase`: Phase offset for seasonal signal (default: 0.0)
- `noise_sigma`: Standard deviation of Gaussian noise (default: 0.1)

# Returns
- Vector of elevation values with trend + seasonality + noise
"""
function synthetic_elevation_timeseries(;
    dates,
    trend=-0.5,
    amplitude=0.2,
    phase=0.0,
    noise_sigma=0.1
)
    t = GlacierHydrology.decimalyear.(dates) .- minimum(GlacierHydrology.decimalyear.(dates))
    signal = trend .* t .+ amplitude .* sin.(2π .* (t .+ phase))
    noise = noise_sigma .* randn(length(t))
    return signal .+ noise
end

"""
    synthetic_geotile_id(lat_min, lat_max, lon_min, lon_max)

Create synthetic geotile ID string for testing.

# Arguments
- `lat_min`, `lat_max`: Latitude bounds [-90, 90]
- `lon_min`, `lon_max`: Longitude bounds [-180, 180]

# Returns
- Geotile ID string in format "lat[+YY+YY]lon[+XXX+XXX]"
"""
function synthetic_geotile_id(lat_min, lat_max, lon_min, lon_max)
    # Sign first, then zero-padded magnitude: latitudes get 2 digits, longitudes 3, and the whole
    # thing is bracketed -- e.g. "lat[+00+02]lon[+028+030]". This has to match the real IDs
    # exactly, because `geotile_extent` slices fixed character positions out of the string.
    # (`lpad("+5", 3, '0')` would give "0+5", padding ahead of the sign, so pad the magnitude.)
    signed(v, width) = (v < 0 ? "-" : "+") * lpad(abs(v), width, '0')
    return string(
        "lat[", signed(lat_min, 2), signed(lat_max, 2), "]",
        "lon[", signed(lon_min, 3), signed(lon_max, 3), "]",
    )
end

"""
    synthetic_river_network(n_nodes=5; topology=:linear)

Generate synthetic river network for routing tests.

# Arguments
- `n_nodes`: Number of river nodes (default: 5)
- `topology`: Network structure (:linear, :branching) (default: :linear)

# Returns
- DataFrame with columns: COMID, NextDownID, lengthkm, lon, lat
"""
function synthetic_river_network(n_nodes=5; topology=:linear)
    ids = collect(1:n_nodes)

    if topology == :linear
        # Linear chain to ocean: 1 → 2 → 3 → 4 → 5 → 0
        next_ids = [i < n_nodes ? i+1 : 0 for i in ids]
        lons = range(-120, -110, length=n_nodes)
        lats = fill(45.0, n_nodes)
    elseif topology == :branching
        # Branching: 1 → 3, 2 → 3 → 4 → 5 → 0
        next_ids = [3, 3, 4, 5, 0]
        lons = [-120, -119, -118, -116, -114]
        lats = [45.5, 44.5, 45.0, 45.0, 45.0]
    else
        error("Unknown topology: $topology")
    end

    lengths_km = fill(10.0, n_nodes)

    return DataFrame(
        COMID=ids,
        NextDownID=next_ids,
        lengthkm=lengths_km,
        lon=lons,
        lat=lats
    )
end

"""
    mock_geotile_dimarray(n_geotiles=3, n_dates=12, n_heights=10)

Create mock DimArray with geotile dimensions for testing.

# Arguments
- `n_geotiles`: Number of geotiles (default: 3)
- `n_dates`: Number of time steps (default: 12)
- `n_heights`: Number of elevation bins (default: 10)

# Returns
- DimArray with dimensions (geotile, date, height)
"""
function mock_geotile_dimarray(n_geotiles=3, n_dates=12, n_heights=10)
    geotiles = [synthetic_geotile_id(45+2i, 47+2i, -123-2i, -121-2i) for i in 0:n_geotiles-1]
    dates = range(DateTime(2018,1,1), step=Month(1), length=n_dates)
    heights = range(1000, 3000, length=n_heights)

    data = randn(n_geotiles, n_dates, n_heights)

    return DimArray(
        data,
        (
            DD.Dim{:geotile}(geotiles),
            DD.Dim{:date}(dates),
            DD.Dim{:height}(heights)
        )
    )
end

"""
    approx_equal(a, b; rtol=1e-6, atol=1e-8)

Check approximate equality with relative and absolute tolerance.

# Arguments
- `a`, `b`: Values to compare
- `rtol`: Relative tolerance (default: 1e-6)
- `atol`: Absolute tolerance (default: 1e-8)

# Returns
- Boolean indicating approximate equality
"""
function approx_equal(a, b; rtol=1e-6, atol=1e-8)
    return isapprox(a, b; rtol=rtol, atol=atol)
end

"""
    write_mock_volume_product(member, dv, area) -> String

Write the variables and attributes that `GlacierVolumeChange.read_volume` reads to
`GlacierVolumeChange.volume_product_path(member)` and return that path.

- `dv`: `DimArray(geotile, date)` in km³; its dates must be a run of `GVC.project_date_bins` centers.
- `area`: `DimArray(geotile, height)` in km²; its heights must be a run of `GVC.project_height_bins` centers.
"""
function write_mock_volume_product(member, dv, area)
    path = GlacierHydrology.GVC.volume_product_path(member)
    attrib = ["surface_mask" => "glacier", "reference_dem" => "best", "curvature_correct" => "true",
              "binning_method" => "nmad5", "fill_param" => 1, "amplitude_correct" => "true",
              "binned_folder" => "mock", "project_id" => "v01"]
    GlacierHydrology.NCDatasets.NCDataset(path, "c"; attrib) do ds
        GlacierHydrology.NCDatasets.defDim(ds, "geotile", size(dv, 1))
        GlacierHydrology.NCDatasets.defDim(ds, "time", size(dv, 2))
        GlacierHydrology.NCDatasets.defDim(ds, "height", size(area, 2))
        GlacierHydrology.NCDatasets.defVar(ds, "geotile", collect(dims(dv, :geotile)), ("geotile",))
        GlacierHydrology.NCDatasets.defVar(ds, "time", DateTime.(collect(dims(dv, :date))), ("time",);
                                           attrib=["units" => "days since 1990-01-01"])
        GlacierHydrology.NCDatasets.defVar(ds, "height", collect(dims(area, :height)), ("height",))
        GlacierHydrology.NCDatasets.defVar(ds, "dv", parent(dv), ("geotile", "time"))
        GlacierHydrology.NCDatasets.defVar(ds, "area", parent(area), ("geotile", "height"))
    end
    return path
end

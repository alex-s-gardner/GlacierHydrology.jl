# Read GEMB_GlacierSims tile NetCDF output into the geotile volume-change ensemble the synthesis
# consumes.
#
# The tile files carry band-resolved output on the same 100 m grid as `project_height_bins`, over a
# (ΔT × precipitation scaling) forcing matrix, for most glacierized 2° tiles. That makes the
# hypsometric gap-filling in `process_gemb_geotiles` unnecessary: there is no simulated elevation
# class, and the only interpolation across the forcing grid is the one `gemb_dv_sample` already
# performs downstream.
#
# The sweep does not reach every geotile this package finds ice in, so a geotile with no tile file of
# its own borrows the nearest one ([`gemb_tile_donors`](@ref)) and is rebinned onto its own hypsometry.
# The borrowed series is a worse estimate the further the donor is and the less its band range overlaps
# the target's ice, which is why `process_gemb_tiles` reports both per substitution.
#
# Bands are weighted by this package's own `_geotile_area_km2`, not the tile file's `band_area`. The two
# are different inventories on the same grid and disagree by a few percent per tile; the altimetry these
# series are calibrated against is integrated over `_geotile_area_km2`, so using it for both is what
# keeps the comparison free of an area bias.

# The ice density a tile was run with, kg m-3, read from the file rather than assumed. Mass fluxes are
# stored in kg m-2 and the height layers in metres, so a mass layer crosses this once on its way to a
# volume -- and it must be the same value GEMB divided by to build `dh_mass`, or the budget will not
# close. Read rather than taken from `δice` so a tile run at a non-default density still converts
# correctly; the two agree at GEMB's default of 917.
function _gemb_density_ice(ds)
    haskey(ds.attrib, "model_density_ice") ||
        error("tile file has no model_density_ice attribute; cannot convert mass fluxes to volume " *
              "consistently with its height layers")
    return Float64(ds.attrib["model_density_ice"])
end

# km³ per (km² × m).
const _GEMB_KM3_PER_KM2_M = 1e-3

# Epoch the tile time axes are reduced to before resampling, so a file written against a different
# `days since` origin still lands on the same axis.
const _GEMB_TIME_EPOCH = DateTime(1950, 1, 1)

# Layers `process_gemb_tiles` builds directly. `smb`, `runoff` and `dv` are derived from these by
# `gemb_add_derived_vars!`, which the driver applies before saving.
const GEMB_TILE_LAYERS = (:fac, :acc, :refreeze, :melt, :rain, :ec)

# Band layers read from each tile file. The `dh` layers are `-cumsum(ice_flux)` and already cumulative;
# the mass fluxes are per-interval sums in kg m-2.
const _GEMB_BAND_HEIGHT_LAYERS = ("dh", "dh_mass")
const _GEMB_BAND_MASS_LAYERS = ("melt", "runoff", "refreeze", "rain", "precipitation",
                                "evaporation_condensation")

"""
    gemb_tile_coverage(tile_dir) -> Dict{String,@NamedTuple{path::String, first::DateTime, last::DateTime, n::Int}}

Geotile id, path and time span of every `.nc` in `tile_dir`.

Ids come from each file's `geotile_id` attribute, which `GEMB_GlacierSims` already writes in this
package's own form (`"lat[+62+64]lon[-152-150]"`), so no filename parsing is involved. Only the time
axis and the attributes are read.

The directory accumulates runs over different windows and output intervals -- a tile from a short
development sweep sits beside a full-record one -- so a caller must select on the span rather than
assume every file is usable. [`gemb_tile_paths`](@ref) does that selection.

Throws if a file carries no `geotile_id`, or if two files claim the same geotile.
"""
function gemb_tile_coverage(tile_dir)
    isdir(tile_dir) || error("no GEMB tile directory at $tile_dir")

    coverage = Dict{String,@NamedTuple{path::String, first::DateTime, last::DateTime, n::Int}}()
    for filename in readdir(tile_dir)
        endswith(filename, ".nc") || continue
        path = joinpath(tile_dir, filename)

        entry = NCDataset(path, "r") do ds
            geotile_id = get(ds.attrib, "geotile_id", nothing)
            isnothing(geotile_id) && error("no geotile_id attribute in $path")
            time = DateTime.(ds["time"][:])
            (geotile_id, (; path, first=first(time), last=last(time), n=length(time)))
        end

        geotile_id, span = entry
        if haskey(coverage, geotile_id)
            error("two tile files claim geotile $geotile_id: $(coverage[geotile_id].path) and $path")
        end
        coverage[geotile_id] = span
    end

    isempty(coverage) && error("no .nc tile files found in $tile_dir")
    return coverage
end

"""
    gemb_tile_paths(tile_dir; covering=nothing) -> Dict{String,String}

Map geotile id to tile file path, keeping only tiles whose time axis spans `covering`.

`covering` is a pair or tuple of the first and last date the ensemble needs. Tiles that do not reach
both ends are dropped: the ensemble is built on one date axis, so a tile covering part of it could only
contribute NaN over the rest. Pass `nothing` to keep every tile.

See [`gemb_tile_coverage`](@ref) for why this filter is not optional in practice.
"""
function gemb_tile_paths(tile_dir; covering=nothing)
    coverage = gemb_tile_coverage(tile_dir)
    isnothing(covering) && return Dict(id => span.path for (id, span) in coverage)

    lo, hi = DateTime(first(covering)), DateTime(last(covering))
    return Dict(id => span.path for (id, span) in coverage
                if span.first <= lo && span.last >= hi)
end

"""
    gemb_tile_own_geotiles(geotile_ids, tile_dir; covering) -> Vector{String}

The ids in `geotile_ids` that have a tile file of their own in `tile_dir` spanning `covering` (see
[`gemb_tile_paths`](@ref)), in their given order. A sweep over a subset of tiles (the pilot, GEMB run
10) is binned on these alone, so no geotile borrows a tile of another.
"""
function gemb_tile_own_geotiles(geotile_ids, tile_dir; covering)
    own = keys(gemb_tile_paths(tile_dir; covering))
    return [id for id in geotile_ids if id in own]
end

"""
    gemb_gradient_supported_geotiles(geotile_ids, area_km2, tile_dir, precipitation_gradient, pscale_nodes;
                                     covering) -> (; supported, excluded)

Split `geotile_ids`, each read from its own tile in `tile_dir` ([`gemb_tile_own_geotiles`](@ref)), by
whether a precipitation gradient of `precipitation_gradient` [% per 100 m] can be represented on its
tile at every geotile-mean pscale node in `pscale_nodes`: every band holding ice under the geotile's
hypsometry (`area_km2`, `(geotile, height)`) stays on the tile's pscale axis
([`gemb_gradient_pscale_range`](@ref)). `excluded` is a `DataFrame` of the others with the range each
supports (`lo > hi` when none). A gradient member is calibrated and scored on `supported` alone.
"""
function gemb_gradient_supported_geotiles(geotile_ids, area_km2, tile_dir, precipitation_gradient,
                                          pscale_nodes; covering)
    paths = gemb_tile_paths(tile_dir; covering)
    supported = String[]
    excluded = DataFrame(id=String[], pscale_lo=Float64[], pscale_hi=Float64[])
    for id in geotile_ids
        haskey(paths, id) || error("$(id) has no tile of its own in $(tile_dir)")
        centers, axis = NCDataset(paths[id], "r") do ds
            Float64.(ds["band_center"][:]), Float64.(ds["precipitation_scaling"][:])
        end
        w, _ = _gemb_band_weights(centers, area_km2[geotile=At(id)])
        lo, hi = gemb_gradient_pscale_range(gemb_band_pscale_multipliers(centers, w, precipitation_gradient),
                                            w, axis)
        if lo <= minimum(pscale_nodes) && maximum(pscale_nodes) <= hi
            push!(supported, String(id))
        else
            push!(excluded, (String(id), lo, hi))
        end
    end
    return (; supported, excluded)
end

# The run parameters a tile file records: GEMB_GlacierSims names them in the `run_parameters` global
# attribute and stores each as a global attribute of its own.
function _gemb_tile_run_parameters(ds)
    haskey(ds.attrib, "run_parameters") || return Dict{String,Any}()
    return Dict{String,Any}(k => ds.attrib[k] for k in split(ds.attrib["run_parameters"])
                            if haskey(ds.attrib, k))
end

# The source trees of the code a tile sweep was run with. Every tile of one sweep records the same
# three, and GEMB_GlacierSims requires them in a calibration's provenance.
const GEMB_TILE_PROVENANCE_REQUIRED = ("software_gemb_tree", "software_gemb_climateforcing_tree",
                                       "software_gemb_glaciersims_tree")

"""
    gemb_tile_provenance(tile_dir, paths; physics_member=nothing) -> Dict{String,Any}

The tile sweep the tile files `paths` (in `tile_dir`) come from, as GEMB_GlacierSims'
`assert_corrections_fitted_on` reads it from a calibration's provenance file:

  * `"tile_dir"`: `realpath(tile_dir)`;
  * `"tile_run_parameters"`: every `software_*` and `physics_*` run parameter the files record, read
    from the files.

Only the `software_*` parameters (package versions and source trees) and the `physics_*` parameters
(the GEMB_GlacierSims physics member and its levels) are recorded. They are the same for every tile of
a sweep, while the `model_*` and downscaling parameters are resolved tile by tile, and
GEMB_GlacierSims compares each recorded key with every reference tile it runs against.

Throws if a file lies outside `tile_dir`, lacks one of [`GEMB_TILE_PROVENANCE_REQUIRED`](@ref),
records a `software_*` or `physics_*` parameter differently from another file (including recording
one another lacks), or was run as a physics member other than `physics_member` (`nothing`: as none).
"""
function gemb_tile_provenance(tile_dir, paths; physics_member=nothing)
    isempty(paths) && error("no GEMB tile files to record the provenance of")
    dir = realpath(tile_dir)
    outside = [p for p in paths if dirname(realpath(p)) != dir]
    isempty(outside) || error("tile files outside $(dir): $(join(outside, ", "))")

    recorded = nothing
    reference = ""
    for path in sort(unique(paths))
        parameters = NCDataset(path, "r") do ds
            _gemb_tile_run_parameters(ds)
        end
        software = Dict{String,Any}(k => v for (k, v) in parameters
                                    if startswith(k, "software_") || startswith(k, "physics_"))
        absent = [k for k in GEMB_TILE_PROVENANCE_REQUIRED if !haskey(software, k)]
        isempty(absent) || error("$(path) records no $(join(absent, ", ")), so the code its run came from is unknown")
        _check_gemb_physics_provenance(software, physics_member, path)

        if isnothing(recorded)
            recorded, reference = software, path
            continue
        end
        for key in sort(collect(Base.union(keys(recorded), keys(software))))
            a, b = get(recorded, key, nothing), get(software, key, nothing)
            a == b || error("GEMB tile files disagree on $(key): $(repr(a)) in $(reference) but $(repr(b)) " *
                            "in $(path); the ensemble would mix tile sweeps")
        end
    end
    return Dict{String,Any}("tile_dir" => dir, "tile_run_parameters" => recorded)
end

# Geotile center as (longitude, latitude), the order `Distances.Haversine` takes.
function _geotile_center(geotile_id)
    extent = geotile_extent(geotile_id)
    return ((extent.X[1] + extent.X[2]) / 2, (extent.Y[1] + extent.Y[2]) / 2)
end

"""
    gemb_tile_donors(geotile_ids, paths; max_distance_km=2000) -> Dict{String,@NamedTuple{donor::String, distance_km::Float64}}

For each of `geotile_ids`, the geotile whose tile file supplies its GEMB series.

A geotile present in `paths` is its own donor at zero distance. One that is not borrows the tile whose
geotile center is nearest by great-circle distance; the borrowed series is then rebinned onto the
target's own hypsometry by [`_gemb_band_weights`](@ref), which interpolates across the donor's bands in
height and holds the outermost band constant beyond them.

A donor is only defensible where the two geotiles share a climate and an elevation range, and neither
condition is checked here -- distance is the whole criterion. `max_distance_km` is the backstop: it is
loose enough for the isolated sub-polar islands the sweep omits, whose nearest neighbor is over a
thousand kilometres away, and tight enough that a sweep which dropped a whole region throws rather than
smearing one surviving tile across it.

Throws if `paths` is empty or if any geotile's nearest donor is further than `max_distance_km`.
"""
function gemb_tile_donors(geotile_ids, paths; max_distance_km=2000)
    isempty(paths) && error("no GEMB tile files to draw donors from")

    donor_ids = collect(keys(paths))
    donor_centers = _geotile_center.(donor_ids)
    # Kilometres, so `max_distance_km` and the reported distances share units.
    haversine = Haversine(6371.0)

    donors = Dict{String,@NamedTuple{donor::String, distance_km::Float64}}()
    too_far = Tuple{String,String,Float64}[]

    for geotile_id in geotile_ids
        if haskey(paths, geotile_id)
            donors[geotile_id] = (; donor=geotile_id, distance_km=0.0)
            continue
        end

        center = _geotile_center(geotile_id)
        distance_km, j = findmin(donor_center -> haversine(center, donor_center), donor_centers)

        if distance_km > max_distance_km
            push!(too_far, (geotile_id, donor_ids[j], distance_km))
        else
            donors[geotile_id] = (; donor=donor_ids[j], distance_km)
        end
    end

    if !isempty(too_far)
        report = join(("$(id) (nearest $(donor) at $(round(Int, d)) km)" for (id, donor, d) in too_far),
                      ", ")
        error("$(length(too_far)) geotiles have no GEMB tile within $(max_distance_km) km: $report")
    end

    return donors
end

"""
    gemb_tile_forcing_grid(path) -> (; pscale, ΔT)

The precipitation-scaling and temperature-offset axes of a tile file, as ascending vectors.

`gemb_dv_sample` interpolates over both with `Gridded(Linear())`, which requires ascending knots, so
the ensemble is built on sorted axes and the tile data is permuted into that order on read.
"""
function gemb_tile_forcing_grid(path)
    return NCDataset(path, "r") do ds
        (; pscale = sort(Float64.(ds["precipitation_scaling"][:])),
           ΔT = sort(Float64.(ds["delta_temperature"][:])))
    end
end

"""
    _gemb_band_weights(band_centers, area_km2) -> (w, extrapolated_fraction)

Collapse "fill the ice-bearing height bins GEMB has no band for, weight by glacier area, sum over
height" into one weight per band, so an area-weighted total becomes a single `series * w`.

`area_km2` is one geotile's hypsometry over `project_height_bins` centers. A bin holding ice that no
band covers is filled by linear interpolation in height across the bands that do, held constant beyond
the outermost band. That fill is linear in the band values
and the area weighting is linear in the filled values, so the composition is a fixed linear functional
of the band series and `w` is that functional.

`w` is in km³ per metre of ice equivalent. `extrapolated_fraction` is the share of the geotile's ice
area lying outside the band range, which a caller should report rather than silently accept.
"""
function _gemb_band_weights(band_centers, area_km2)
    issorted(band_centers) || error("band centers are not ascending: $(band_centers)")

    height_centers = collect(val(dims(area_km2, :height)))
    areas = Float64.(collect(area_km2))
    ice = findall(>(0), areas)
    isempty(ice) && error("geotile hypsometry holds no ice")

    n_band = length(band_centers)
    w = zeros(Float64, n_band)
    extrapolated_area = 0.0

    for i in ice
        height = height_centers[i]
        weight = areas[i] * _GEMB_KM3_PER_KM2_M

        j = searchsortedfirst(band_centers, height)
        if j <= n_band && band_centers[j] == height
            # A band covers this bin exactly, which holds for all but a handful of bins.
            w[j] += weight
        elseif j == 1 || j > n_band
            # Outside the band range: hold the nearest band constant.
            w[clamp(j, 1, n_band)] += weight
            extrapolated_area += areas[i]
        else
            lo, hi = band_centers[j-1], band_centers[j]
            f = (height - lo) / (hi - lo)
            w[j-1] += weight * (1 - f)
            w[j] += weight * f
        end
    end

    return w, extrapolated_area / sum(areas[ice])
end

# Cumulative levels resampled onto `target_days`. Linear interpolation of the level is the right
# reduction: the increment over a target interval is then the integral of the flux across it.
function _gemb_resample(series, source_days, target_days)
    return DataInterpolations.LinearInterpolation(series, source_days).(target_days)
end

"""
    read_gemb_tile(path, area_km2, target_days; n_pscale, n_ΔT, grid=nothing, subset=false,
                   precipitation_gradient=0.0)

One tile's volume-change ensemble as `(; layers, extrapolated_fraction, closure_km3, melt_offset_km3)`.

`layers` is a `Dict{Symbol,Array{Float64,3}}` over `(date, pscale, ΔT)` in km³ of ice equivalent,
cumulative and referenced to the first target date. `area_km2` is the geotile's row of
`_geotile_area_km2`; `target_days` are the output date centers as days since `_GEMB_TIME_EPOCH`.
The file's forcing grid must be `n_ΔT` × `n_pscale`, and, when `grid` (a
[`gemb_tile_forcing_grid`](@ref)) is given, must have exactly its axes. With `subset = true` the file
may carry more nodes: only those of `grid` are read, matched by value, and a node of `grid` the file
lacks is an error.

# Precipitation gradient

With `precipitation_gradient = g` [% per 100 m] each band is read at its own precipitation scaling:
at an ensemble node `pscale`, band `b` takes `pscale · m_b` ([`gemb_band_pscale_multipliers`](@ref),
mean-centred on the geotile's band weights), its series interpolated linearly between the file's two
bracketing pscale nodes at the node's ΔT. Every band's series is a full GEMB column at its own forcing
correction, so this is the gradient GEMB would give with total precipitation scaled by `m_b`, up to the
pscale interpolation. A band holding ice whose `pscale · m_b` falls outside the file's pscale axis is an
error that names the geotile-mean range the tile supports ([`gemb_gradient_pscale_range`](@ref)). The
bands are summed only after that, so at `g = 0` this is the band sum read node by node.

# Layer construction

`dv` and `dv_mass` are the area-weighted band `dh` and `dh_mass`. Every mass flux is split by where it
ended up before being accumulated, over each output interval and in each band `b`, and only then
area-weighted over bands (`Σ_b` below, with the weights of [`_gemb_band_weights`](@ref)):

    rain_retained = Σ_b max(rain_b - runoff_b, 0)   -- rain that never left, so it refroze in the ice
    rain_shed     = Σ_b rain_b - rain_retained      -- rain that ran off
    melt_runoff   = Σ_b max(runoff_b - rain_b, 0)   -- the part of runoff that is meltwater

Each band is a separate GEMB column, so the split holds band by band. Splitting the area-weighted
tile totals instead would let rain retained in one band offset meltwater runoff in another, and
understate both `melt_runoff` and `acc` by the same amount.

The stored layers are

    acc      = cumulative (precipitation - rain + rain_retained)   -- snowfall plus refrozen rain
    rain     = cumulative rain_shed
    refreeze = cumulative refreeze                                 -- GEMB's, unaltered
    melt     = cumulative (melt_runoff + refreeze)
    fac      = dv - dv_mass
    ec       = -cumulative evaporation_condensation

This matches the convention the altimetry products are on: **accumulation is snowfall plus the rain that
refroze**, and **runoff excludes rain**, so `gemb_add_derived_vars!`'s `runoff = melt - refreeze` comes out
as `melt_runoff` alone. Neither can be read straight off the file. GEMB's `runoff` carries rain as well
as meltwater despite its `long_name`, which the water balance settles:
`melt + rain - runoff - refreeze - Δstored` closes while the rain-free version misses by exactly the
rain. GEMB's `melt` is separately net of rain (`melt_total = max(0, melt_sum - rain)` in
`calculate_melt.jl`); `melt_offset_km3` reports how far the layer above sits from it, which is a percent
or so.

The split loses nothing. `max(0, a-b) - max(0, b-a) == a - b` for any band and interval, and the band
sum is linear, so the retained and shed halves cancel in `smb`, which reduces to `precipitation - runoff + evaporation_condensation` -- GEMB's
own surface mass balance. Hence

  * `smb == dv_mass`, and `ec` is the real evaporation and condensation rather than a closure residual,
  * `acc`, `rain`, `refreeze`, `melt` and the derived `runoff` are all monotonic by construction, with
    nothing clamped away,
  * `dv == Σ dh·area`, and `dm = (dv - fac)·δice` is exactly GEMB's mass term.

`closure_km3` is the largest `|smb - dv_mass|` over the forcing grid and should sit at the rounding
floor; a nonzero value means one of the assumptions above no longer holds for that tile.

!!! note "Ice density"
    Mass fluxes are divided by the tile's own `model_density_ice`, which is what GEMB divided by to build
    `dh_mass`. Using any other value leaves a closure error of order a percent of the accumulated fluxes
    -- tens of km³ on a large tile. `δice` is set to GEMB's default of 917 so that the volumes here and
    the `volume2mass = δice/1000` applied downstream describe the same ice.
"""
function read_gemb_tile(path, area_km2, target_days; n_pscale, n_ΔT, grid=nothing, subset=false,
                        precipitation_gradient=0.0)
    return NCDataset(path, "r"; maskingvalue = NaN) do ds
        band_centers = Float64.(ds["band_center"][:])
        w, extrapolated_fraction = _gemb_band_weights(band_centers, area_km2)
        w_mass = w ./ _gemb_density_ice(ds)

        source_days = Dates.value.(DateTime.(ds["time"][:]) .- _GEMB_TIME_EPOCH) ./ 86_400_000
        if first(source_days) > first(target_days) || last(source_days) < last(target_days)
            error("$(basename(path)) spans $(first(source_days))..$(last(source_days)) days since " *
                  "$(_GEMB_TIME_EPOCH) but $(first(target_days))..$(last(target_days)) was requested")
        end

        # The ensemble axes are ascending; the file's are not necessarily, so read through a permutation.
        ΔT_stored = Float64.(ds["delta_temperature"][:])
        pscale_stored = Float64.(ds["precipitation_scaling"][:])
        if subset
            isnothing(grid) && throw(ArgumentError("subset = true needs the grid to read"))
            function node(axis, v, name)
                i = findfirst(==(v), axis)
                isnothing(i) && error("$(basename(path)) has no $name node $v; its axis is $(sort(axis))")
                return i
            end
            iΔT = [node(ΔT_stored, v, "ΔT") for v in grid.ΔT]
            ipscale = [node(pscale_stored, v, "pscale") for v in grid.pscale]
        else
            iΔT = sortperm(ΔT_stored)
            ipscale = sortperm(pscale_stored)
        end
        length(iΔT) == n_ΔT && length(ipscale) == n_pscale ||
            error("$(basename(path)) forcing grid is $(length(iΔT))x$(length(ipscale)), expected " *
                  "$(n_ΔT)x$(n_pscale)")
        # Same size is not same grid: the ensemble's axes are the first tile's, and a tile on other
        # nodes would be read onto them silently.
        if !isnothing(grid) && !subset
            ΔT_file = sort(ΔT_stored)
            pscale_file = sort(pscale_stored)
            ΔT_file == grid.ΔT && pscale_file == grid.pscale ||
                error("$(basename(path)) forcing grid is ΔT $(ΔT_file) × pscale $(pscale_file), but the " *
                      "ensemble's is ΔT $(grid.ΔT) × pscale $(grid.pscale); every tile of a sweep must share it")
        end

        # Read each band layer once, whole. Slicing a deflated variable per forcing node instead costs
        # one chunk decompression per slice and dominates the runtime -- 49 nodes over 7 layers is two
        # orders of magnitude more decompression than the file needs. Kept at the stored `Float32`:
        # about 45 MB per layer, and the arithmetic below promotes against the `Float64` weights.
        band = Dict(name => ds[name][:, :, :, :] for name in
                    (_GEMB_BAND_HEIGHT_LAYERS..., _GEMB_BAND_MASS_LAYERS...))
        m = gemb_band_pscale_multipliers(band_centers, w, precipitation_gradient)
        gradient_nodes = iszero(precipitation_gradient) ? nothing :
            _gemb_gradient_band_nodes(path, band_centers, w, m, pscale_stored, unique(ipscale))
        n_time = length(target_days)
        layers = Dict(k => fill(NaN, n_time, n_pscale, n_ΔT)
                      for k in (GEMB_TILE_LAYERS..., :dv, :dv_mass))
        closure_km3 = 0.0
        melt_offset_km3 = 0.0

        for (jp, p) in enumerate(ipscale), (jt, t) in enumerate(iΔT)
            # `w` folds the height fill and the area weighting into one vector, so each layer's
            # area-weighted total is a single matrix-vector product over bands.
            # With a precipitation gradient each band is read at its own pscale; the band sum follows.
            slab = isnothing(gradient_nodes) ? (name -> view(band[name], :, :, t, p)) :
                   (name -> _gemb_band_pscale_slab(band[name], t, gradient_nodes[p]))
            total(name, weights) = slab(name) * weights

            precipitation = total("precipitation", w_mass)
            rain = total("rain", w_mass)
            refreeze = total("refreeze", w_mass)

            # Split the rain by where it ended up. GEMB's `runoff` carries rain as well as meltwater
            # despite its name, so rain in excess of what left the column over an interval must have
            # stayed in the ice: that part is accumulation, and the remainder ran off. Both halves are
            # per-interval quantities, hence the split before any cumulative sum. The split is per
            # band, before the area weighting: each band is its own column, and splitting the tile
            # totals instead would let rain retained in one band cancel meltwater runoff in another.
            rain_band = slab("rain")
            runoff_band = slab("runoff")
            rain_retained = max.(Float64.(rain_band) .- runoff_band, 0.0) * w_mass
            rain_shed = rain .- rain_retained
            melt_runoff = max.(Float64.(runoff_band) .- rain_band, 0.0) * w_mass

            native = Dict{Symbol,Vector{Float64}}(
                :dv => total("dh", w),
                :dv_mass => total("dh_mass", w),
                # Snowfall plus the rain that refroze in the ice.
                :acc => cumsum(precipitation .- rain .+ rain_retained),
                # Only the rain that left; the retained part is in `acc`, not counted twice.
                :rain => cumsum(rain_shed),
                :refreeze => cumsum(refreeze),
                # So `gemb_add_derived_vars!`'s `melt - refreeze` is runoff net of rain.
                :melt => cumsum(melt_runoff .+ refreeze),
                # Positive `ec` is mass loss here, the reverse of GEMB's flux convention.
                :ec => -cumsum(total("evaporation_condensation", w_mass)),
                # GEMB's own melt layer, kept only to report how far the above sits from it.
                :melt_gemb => cumsum(total("melt", w_mass)),
            )

            series = Dict(k => _gemb_resample(v, source_days, target_days) for (k, v) in native)
            for v in values(series)
                v .-= first(v)
            end

            series[:fac] = series[:dv] .- series[:dv_mass]
            smb = series[:acc] .- series[:melt] .+ series[:refreeze] .- series[:ec]
            closure_km3 = max(closure_km3, maximum(abs, smb .- series[:dv_mass]))
            melt_offset_km3 = max(melt_offset_km3, abs(last(series[:melt]) - last(series[:melt_gemb])))

            for k in (GEMB_TILE_LAYERS..., :dv, :dv_mass)
                layers[k][:, jp, jt] = series[k]
            end
        end

        return (; layers, extrapolated_fraction, closure_km3, melt_offset_km3)
    end
end

# For each ensemble pscale node (a file index into `pscale_stored`), each band's pair of bracketing file
# nodes and the weight on the upper one, at `pscale · m_b`. A band without ice reads the node itself.
function _gemb_gradient_band_nodes(path, band_centers, w, m, pscale_stored, ensemble_nodes)
    order = sortperm(pscale_stored)
    axis = pscale_stored[order]
    lo_ok, hi_ok = gemb_gradient_pscale_range(m, w, axis)
    nodes = Dict{Int,Vector{Tuple{Int,Int,Float64}}}()
    for p in ensemble_nodes
        P = pscale_stored[p]
        nodes[p] = map(eachindex(band_centers)) do b
            w[b] > 0 || return (p, p, 0.0)
            q = P * m[b]
            first(axis) <= q <= last(axis) ||
                error("$(basename(path)): a precipitation gradient puts the $(band_centers[b]) m band at " *
                      "pscale $(round(q; digits=4)) (multiplier $(round(m[b]; digits=4))) at the ensemble node " *
                      "pscale $(P), outside the tile's axis $(first(axis))..$(last(axis)); at this gradient " *
                      "the tile supports geotile-mean pscale $(round(lo_ok; digits=3))..$(round(hi_ok; digits=3)) only")
            j = searchsortedfirst(axis, q)
            axis[j] == q && return (order[j], order[j], 0.0)
            (order[j-1], order[j], (q - axis[j-1]) / (axis[j] - axis[j-1]))
        end
    end
    return nodes
end

# One layer's (time, band) slab at ΔT index `t`, each band interpolated between its own pscale nodes.
function _gemb_band_pscale_slab(A, t, band_nodes)
    out = Matrix{Float64}(undef, size(A, 1), length(band_nodes))
    for (b, (lo, hi, a)) in enumerate(band_nodes)
        if iszero(a)
            out[:, b] .= view(A, :, b, t, lo)
        else
            out[:, b] .= (1 - a) .* view(A, :, b, t, lo) .+ a .* view(A, :, b, t, hi)
        end
    end
    return out
end

"""
    process_gemb_tiles(geotiles, area_km2; tile_dir, date_center, show_stats=true, max_donor_distance_km=2000,
                       physics_member=nothing, forcing_grid=nothing)

Build the GEMB volume-change ensemble from `GEMB_GlacierSims` tile NetCDF.

Returns a `DimStack` of [`GEMB_TILE_LAYERS`](@ref) over `(geotile, date, pscale, ΔT)` in km³ of ice
equivalent. `smb`, `runoff` and `dv` are added later by `gemb_add_derived_vars!`. The stack's metadata `"tile_provenance"` is the [`gemb_tile_provenance`](@ref) of the tile files
read, donors included, which must all have been run as `physics_member` (see
[`gemb_tile_provenance`](@ref)), and every tile must share the first tile's forcing grid.

`forcing_grid = (; pscale, ΔT)` builds the ensemble on those nodes instead, read by value from tiles
that may carry more (a dense-grid tile beside reduced-grid ones); every tile must hold every node.

A `physics_member` with a precipitation gradient (`"<tile member>-g<g>"`,
[`gemb_physics_member_parts`](@ref)) reads the tile member's sweep with each band at its own pscale
([`read_gemb_tile`](@ref)), and its provenance names the gradient member and its
`physics_precipitation_gradient`, so the ensemble and its fits stay apart from the tile member's.

Every band is weighted by `area_km2`, so the ensemble and the altimetry it is calibrated against share
one hypsometry. A geotile with no tile file of its own borrows the nearest one within
`max_donor_distance_km` -- see [`gemb_tile_donors`](@ref) -- and is rebinned onto its own hypsometry, so
every geotile passed in gets a series rather than an all-NaN row that would reach the ensemble's NaN
check with no indication of why.

`show_stats` prints the substitutions worst-area first, the geotiles whose ice extends beyond the band
range they were read against, and the largest `ec` closure residual. All three are silent data-quality
problems otherwise, and the first two compound: a borrowed tile whose bands miss the target's ice
entirely is held constant from its nearest band, which is the weakest estimate this produces.
"""
function process_gemb_tiles(geotiles, area_km2; tile_dir, date_center, show_stats=true,
                            max_donor_distance_km=2000, physics_member=nothing, forcing_grid=nothing)
    paths = gemb_tile_paths(tile_dir; covering=extrema(date_center))
    donors = gemb_tile_donors(geotiles.id, paths; max_distance_km=max_donor_distance_km)
    parts = isnothing(physics_member) ? (; tile_member=nothing, precipitation_gradient=0.0) :
            gemb_physics_member_parts(physics_member)
    # From the files actually read, before the read: a mixed sweep fails here rather than after it.
    tile_provenance = _gemb_gradient_provenance(
        gemb_tile_provenance(tile_dir, unique(paths[d.donor] for d in values(donors));
                             physics_member=parts.tile_member), parts.precipitation_gradient)

    subset = !isnothing(forcing_grid)
    grid = subset ?
        (; pscale=sort(Float64.(collect(forcing_grid.pscale))), ΔT=sort(Float64.(collect(forcing_grid.ΔT)))) :
        gemb_tile_forcing_grid(paths[donors[first(geotiles.id)].donor])
    subset && (allunique(grid.pscale) && allunique(grid.ΔT) ||
               throw(ArgumentError("forcing_grid repeats a node: $forcing_grid")))
    ddate = Dim{:date}(collect(date_center))
    dgeotile = Dim{:geotile}(collect(geotiles.id))
    dpscale = Dim{:pscale}(grid.pscale)
    dΔT = Dim{:ΔT}(grid.ΔT)

    gemb_dv = DimStack([DimArray(fill(NaN, (dgeotile, ddate, dpscale, dΔT)); name=k)
                        for k in GEMB_TILE_LAYERS]...; metadata=Dict("tile_provenance" => tile_provenance))

    target_days = Dates.value.(DateTime.(date_center) .- _GEMB_TIME_EPOCH) ./ 86_400_000
    extrapolated = fill(NaN, dgeotile)
    melt_offset = fill(NaN, dgeotile)
    closure = fill(NaN, dgeotile)

    @showprogress desc = "Reading GEMB tiles" Threads.@threads for row in collect(eachrow(geotiles))
        tile = read_gemb_tile(paths[donors[row.id].donor], area_km2[geotile=At(row.id)], target_days;
                              n_pscale=length(grid.pscale), n_ΔT=length(grid.ΔT), grid, subset,
                              parts.precipitation_gradient)

        for k in GEMB_TILE_LAYERS
            gemb_dv[k][geotile=At(row.id)] = tile.layers[k]
        end

        extrapolated[At(row.id)] = tile.extrapolated_fraction
        melt_offset[At(row.id)] = tile.melt_offset_km3
        # `dv` reconstructed through this package's identities against the one the tile file implies.
        smb = tile.layers[:acc] .- tile.layers[:melt] .+ tile.layers[:refreeze] .- tile.layers[:ec]
        closure[At(row.id)] = maximum(abs, (smb .+ tile.layers[:fac]) .- tile.layers[:dv])
    end

    if show_stats
        printstyled("GEMB tiles read: $(length(dgeotile))\n"; color=:light_green)
        @printf("  max |dv| closure error : %.3g km3\n", maximum(closure))
        @printf("  melt vs GEMB's melt    : up to %.1f km3 (%s) -- the rain GEMB nets off its melt\n",
                maximum(melt_offset), val(dgeotile)[argmax(melt_offset)])

        borrowed = [id for id in val(dgeotile) if donors[id].donor != id]
        if isempty(borrowed)
            println("  every geotile has a GEMB tile of its own")
        else
            ice_km2 = Dict(id => sum(area_km2[geotile=At(id)]) for id in val(dgeotile))
            sort!(borrowed; by=id -> ice_km2[id], rev=true)
            total_km2 = sum(values(ice_km2))
            @printf("  %d geotiles have no tile of their own and borrow the nearest (%.0f km2, %.3f%% of the ice):\n",
                    length(borrowed), sum(ice_km2[id] for id in borrowed),
                    100 * sum(ice_km2[id] for id in borrowed) / total_km2)
            for id in first(borrowed, 10)
                @printf("      %-28s %8.2f km2 <- %-28s %5.0f km, %5.1f%% of ice beyond its bands\n",
                        id, ice_km2[id], donors[id].donor, donors[id].distance_km,
                        100 * extrapolated[At(id)])
            end
            length(borrowed) > 10 && println("      ... $(length(borrowed) - 10) more, all smaller")
        end

        beyond = findall(>(0), collect(extrapolated))
        if isempty(beyond)
            println("  every ice-bearing height bin is covered by a GEMB band")
        else
            order = beyond[sortperm(collect(extrapolated)[beyond]; rev=true)]
            @printf("  %d geotiles have ice beyond the band range (median %.2f%%):\n",
                    length(beyond), 100 * median(collect(extrapolated)[beyond]))
            for i in first(order, 10)
                @printf("      %-28s %6.2f%% of ice area\n",
                        val(dgeotile)[i], 100 * extrapolated[i])
            end
        end
    end

    return gemb_dv
end

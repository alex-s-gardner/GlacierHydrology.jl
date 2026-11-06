# WGMS Fluctuations of Glaciers (FoG) seasonal balances at specific elevations, matched to the GEMB
# band series this package calibrates against altimetry, and scored.
#
# Only balances recorded for an elevation are used: stake and pit balances (`mass_balance_point.csv`)
# and elevation-band balances (`mass_balance_band.csv`). The glacier-wide balances of
# `mass_balance.csv` are never read; that table supplies the survey dates of a glacier-year, which the
# band table does not carry, and the glacier area the records are classified by.
#
# Conventions are those of the FoG 2025-02 data package (its `datapackage.json` and readme): balances in
# m w.e., elevations in m, densities in kg m⁻³, dates as `yyyy-mm-dd` with an uncertainty in days. An
# unknown day is stored mid-month with an uncertainty of 14.5–15 days and an unknown month mid-year with
# 182–182.5 days, so a date's uncertainty, not the date alone, says whether it bounds a season.
#
# The reader validates the schema and throws on anything that is not a data-quality judgement. Quality
# exclusions are returned as records with a reason, so every count can be reported.

const WGMS_BALANCE_CODES = ("annual", "summer", "winter", "index")
const WGMS_TIME_SYSTEMS = ("floating", "fixed", "stratigraphic", "combined", "other")
const WGMS_GTNG_REGIONS = ("01_alaska", "02_western_canada_usa", "03_arctic_canada_north",
    "04_arctic_canada_south", "05_greenland_periphery", "06_iceland", "07_svalbard_jan_mayen",
    "08_scandinavia", "09_russian_arctic", "10_north_asia", "11_central_europe",
    "12_caucasus_middle_east", "13_central_asia", "14_south_asia_west", "15_south_asia_east",
    "16_low_latitudes", "17_southern_andes", "18_new_zealand", "19_subantarctic_antarctic_islands",
    "20_antarctic_mainland")

const _WGMS_GLACIER_COLUMNS = ("id", "short_name", "latitude", "longitude", "gtng_region")
const _WGMS_GLACIER_YEAR_COLUMNS = ("glacier_id", "year", "time_system", "begin_date", "begin_date_unc",
    "midseason_date", "midseason_date_unc", "end_date", "end_date_unc", "area")
const _WGMS_POINT_COLUMNS = ("glacier_id", "year", "id", "begin_date", "begin_date_unc", "end_date",
    "end_date_unc", "latitude", "longitude", "elevation", "balance", "balance_unc", "density",
    "balance_code", "remarks")
const _WGMS_BAND_COLUMNS = ("glacier_id", "year", "lower_elevation", "upper_elevation", "area",
    "winter_balance", "winter_balance_unc", "summer_balance", "summer_balance_unc", "remarks")

# Every column as text, so each field is parsed and checked here with its row number in the message.
function _wgms_read_table(path, columns)
    isfile(path) || error("no WGMS table at $path")
    df = CSV.read(path, DataFrame; types=String, missingstring="", strict=true)
    absent = [c for c in columns if !(c in names(df))]
    isempty(absent) || error("$(basename(path)) lacks column(s) $(join(absent, ", "))")
    return df
end

function _wgms_parse(path, df, column, T; required=false, lo=nothing, hi=nothing, values=nothing)
    out = Vector{Union{Missing,T}}(missing, nrow(df))
    for (i, s) in enumerate(df[!, column])
        if ismissing(s)
            required && error("$(basename(path)) row $(i + 1): $column is empty")
            continue
        end
        v = if T === Date
            d = tryparse(Date, s, dateformat"yyyy-mm-dd")
            isnothing(d) && error("$(basename(path)) row $(i + 1): $column = \"$s\" is not a yyyy-mm-dd date")
            d
        elseif T === String
            s
        else
            x = tryparse(T, s)
            # Integer columns are written as e.g. "910" or, in some exports, "910.0".
            if isnothing(x) && T <: Integer
                f = tryparse(Float64, s)
                x = (!isnothing(f) && isinteger(f)) ? T(f) : nothing
            end
            isnothing(x) && error("$(basename(path)) row $(i + 1): $column = \"$s\" is not a $T")
            x
        end
        if !isnothing(values) && !(v in values)
            error("$(basename(path)) row $(i + 1): $column = \"$v\" is not one of $(join(values, ", "))")
        end
        !isnothing(lo) && v < lo && error("$(basename(path)) row $(i + 1): $column = $v is below $lo")
        !isnothing(hi) && v > hi && error("$(basename(path)) row $(i + 1): $column = $v is above $hi")
        out[i] = v
    end
    return required ? convert(Vector{T}, out) : out
end

function _wgms_unique_key(path, df, columns)
    seen = Dict{Any,Int}()
    for (i, row) in enumerate(eachrow(df[!, collect(columns)]))
        k = Tuple(row)
        haskey(seen, k) && error("$(basename(path)): rows $(seen[k] + 1) and $(i + 1) share the key " *
                                 "$(join(columns, ", ")) = $k")
        seen[k] = i
    end
end

"""
    read_wgms(dir) -> (; glaciers, glacier_years, points, bands)

The four WGMS FoG tables this package's seasonal validation uses, read from the `data/` directory of
an FoG 2025-02 data package and validated.

  * `glaciers` (`glacier.csv`): `id`, `short_name`, `latitude`, `longitude`, `gtng_region` and `region`,
    its number (GTN-G regions are numbered as the RGI first-order regions), and `area_km2`, the median
    over the glacier's years of the area of the glacier-wide balance (`missing` where no year reports
    one).
  * `glacier_years` (`mass_balance.csv`): the survey dates of each glacier-year, its `time_system` and
    its `area` [m²]. The glacier-wide balance columns are not read.
  * `points` (`mass_balance_point.csv`): every stake and pit record, all `balance_code`s.
  * `bands` (`mass_balance_band.csv`): every elevation-band record.

Throws on a missing column, an unparseable field, a value outside the documented enumeration of
`balance_code`, `time_system` or `gtng_region`, a latitude or longitude out of range, a negative
uncertainty, a duplicated primary key, or a point or band whose glacier (or glacier-year, for bands)
is not in the glacier (or glacier-year) table.
"""
function read_wgms(dir)
    gpath, ypath = joinpath(dir, "glacier.csv"), joinpath(dir, "mass_balance.csv")
    ppath, bpath = joinpath(dir, "mass_balance_point.csv"), joinpath(dir, "mass_balance_band.csv")

    g = _wgms_read_table(gpath, _WGMS_GLACIER_COLUMNS)
    glaciers = DataFrame(
        id=_wgms_parse(gpath, g, "id", Int; required=true, lo=0),
        short_name=_wgms_parse(gpath, g, "short_name", String; required=true),
        latitude=_wgms_parse(gpath, g, "latitude", Float64; required=true, lo=-90, hi=90),
        longitude=_wgms_parse(gpath, g, "longitude", Float64; required=true, lo=-180, hi=180),
        gtng_region=_wgms_parse(gpath, g, "gtng_region", String; required=true, values=WGMS_GTNG_REGIONS))
    glaciers.region = parse.(Int, first.(glaciers.gtng_region, 2))
    _wgms_unique_key(gpath, glaciers, (:id,))

    y = _wgms_read_table(ypath, _WGMS_GLACIER_YEAR_COLUMNS)
    glacier_years = DataFrame(
        glacier_id=_wgms_parse(ypath, y, "glacier_id", Int; required=true, lo=0),
        year=_wgms_parse(ypath, y, "year", Int; required=true),
        time_system=_wgms_parse(ypath, y, "time_system", String; values=WGMS_TIME_SYSTEMS),
        area=_wgms_parse(ypath, y, "area", Float64; lo=0))
    for f in ("begin_date", "midseason_date", "end_date")
        glacier_years[!, f] = _wgms_parse(ypath, y, f, Date)
        glacier_years[!, f*"_unc"] = _wgms_parse(ypath, y, f * "_unc", Float64; lo=0)
    end
    _wgms_unique_key(ypath, glacier_years, (:glacier_id, :year))
    areas = Dict{Int,Vector{Float64}}()
    for (id, a) in zip(glacier_years.glacier_id, glacier_years.area)
        ismissing(a) || push!(get!(() -> Float64[], areas, id), a)
    end
    glaciers.area_km2 = Union{Missing,Float64}[haskey(areas, id) ? median(areas[id]) / 1e6 : missing
                                                for id in glaciers.id]

    p = _wgms_read_table(ppath, _WGMS_POINT_COLUMNS)
    points = DataFrame(
        id=_wgms_parse(ppath, p, "id", Int; required=true),
        glacier_id=_wgms_parse(ppath, p, "glacier_id", Int; required=true, lo=0),
        year=_wgms_parse(ppath, p, "year", Int; required=true),
        begin_date=_wgms_parse(ppath, p, "begin_date", Date; required=true),
        begin_date_unc=_wgms_parse(ppath, p, "begin_date_unc", Float64; required=true, lo=0),
        end_date=_wgms_parse(ppath, p, "end_date", Date; required=true),
        end_date_unc=_wgms_parse(ppath, p, "end_date_unc", Float64; required=true, lo=0),
        latitude=_wgms_parse(ppath, p, "latitude", Float64; lo=-90, hi=90),
        longitude=_wgms_parse(ppath, p, "longitude", Float64; lo=-180, hi=180),
        elevation=_wgms_parse(ppath, p, "elevation", Float64),
        balance=_wgms_parse(ppath, p, "balance", Float64),
        balance_unc=_wgms_parse(ppath, p, "balance_unc", Float64; lo=0),
        density=_wgms_parse(ppath, p, "density", Float64; lo=0),
        balance_code=_wgms_parse(ppath, p, "balance_code", String; required=true, values=WGMS_BALANCE_CODES),
        remarks=_wgms_parse(ppath, p, "remarks", String))
    _wgms_unique_key(ppath, points, (:id,))

    b = _wgms_read_table(bpath, _WGMS_BAND_COLUMNS)
    bands = DataFrame(
        glacier_id=_wgms_parse(bpath, b, "glacier_id", Int; required=true, lo=0),
        year=_wgms_parse(bpath, b, "year", Int; required=true),
        lower_elevation=Float64.(_wgms_parse(bpath, b, "lower_elevation", Int; required=true)),
        upper_elevation=Float64.(_wgms_parse(bpath, b, "upper_elevation", Int; required=true)),
        area=_wgms_parse(bpath, b, "area", Float64; lo=0),
        winter_balance=_wgms_parse(bpath, b, "winter_balance", Float64),
        winter_balance_unc=_wgms_parse(bpath, b, "winter_balance_unc", Float64; lo=0),
        summer_balance=_wgms_parse(bpath, b, "summer_balance", Float64),
        summer_balance_unc=_wgms_parse(bpath, b, "summer_balance_unc", Float64; lo=0),
        remarks=_wgms_parse(bpath, b, "remarks", String))
    _wgms_unique_key(bpath, bands, (:glacier_id, :year, :lower_elevation, :upper_elevation))

    known = Set(glaciers.id)
    for (name, t) in (("mass_balance_point.csv", points), ("mass_balance_band.csv", bands),
                      ("mass_balance.csv", glacier_years))
        orphans = setdiff(Set(t.glacier_id), known)
        isempty(orphans) || error("$name refers to glacier id(s) absent from glacier.csv: " *
                                  join(sort(collect(orphans))[1:min(10, end)], ", "))
    end
    years = Set(zip(glacier_years.glacier_id, glacier_years.year))
    undated = setdiff(Set(zip(bands.glacier_id, bands.year)), years)
    isempty(undated) || error("mass_balance_band.csv has glacier-years absent from mass_balance.csv: " *
                              join(sort(collect(undated))[1:min(10, end)], ", "))

    return (; glaciers, glacier_years, points, bands)
end

"""
    wgms_geotile_id(latitude, longitude; width=2) -> String

The id of the `width`° geotile that contains (`latitude`, `longitude`), in this package's form
(`"lat[+46+48]lon[+008+010]"`). A point on a tile edge belongs to the tile north and east of it; the
north pole and the antimeridian at +180° belong to the last tile.
"""
function wgms_geotile_id(latitude, longitude; width=2)
    -90 <= latitude <= 90 || error("latitude $latitude is outside [-90, 90]")
    -180 <= longitude <= 180 || error("longitude $longitude is outside [-180, 180]")
    lat0 = min(floor(Int, latitude / width) * width, 90 - width)
    lon0 = min(floor(Int, longitude / width) * width, 180 - width)
    signed(v, n) = (v < 0 ? "-" : "+") * lpad(abs(v), n, '0')
    return "lat[$(signed(lat0, 2))$(signed(lat0 + width, 2))]lon[$(signed(lon0, 3))$(signed(lon0 + width, 3))]"
end

"""
    WGMS_SEASONAL_QC

Default quality-control settings of [`wgms_seasonal_points`](@ref) and
[`wgms_seasonal_bands`](@ref):

  * `max_date_unc_days = 15`: a date with a larger uncertainty is a placeholder for an unknown month
    (182–182.5 d) and does not bound a season; 15 keeps dates of unknown day (14.5–15 d).
  * `window_days`: the accepted lengths of a winter and a summer window. Longer "winter" windows are
    whole years, shorter ones zero-length or corrupt.
  * `balance_limits` [m w.e.]: plausibility bounds per season.
  * `max_balance_unc = 0.5` [m w.e.]: a larger reported random error is excluded.
  * `winter_density_limits = (100, 700)` [kg m⁻³]: outside this a reported winter density is not a
    seasonal snowpack bulk density.
  * `max_point_distance_km = 50`: a point further than this from its glacier's coordinates is
    mislocated (ice caps put real points 35 km out).
  * `duplicate_tolerance = 0.01` [m w.e.]: duplicates whose balances agree within this are kept once;
    any other duplicate group is excluded whole.
  * `max_band_width = 200` [m].
  * `excluded_remarks`: case-insensitive remark substrings that exclude a record (debris cover, which
    GEMB does not model; reconstructed values, which are not measurements).
  * `assumed_density_remarks`: case-insensitive remark substrings that mark a reported winter density
    as assumed rather than measured (`density_class = "assumed"`), e.g. "estimated density of snow
    550 kg/m3" or "For density assumptions … see Machguth et al. (2016)".
"""
const WGMS_SEASONAL_QC = (;
    max_date_unc_days=15.0,
    window_days=(winter=(30, 300), summer=(30, 240)),
    balance_limits=(winter=(-3.0, 10.0), summer=(-15.0, 5.0)),
    max_balance_unc=0.5,
    winter_density_limits=(100.0, 700.0),
    max_point_distance_km=50.0,
    duplicate_tolerance=0.01,
    max_band_width=200.0,
    excluded_remarks=("debris", "reconstruct"),
    assumed_density_remarks=("assum", "estimat", "standard"),
)

_wgms_remark_matches(r, keys) = !ismissing(r) && any(k -> occursin(k, lowercase(r)), keys)
_wgms_remark_excluded(r, qc) = _wgms_remark_matches(r, qc.excluded_remarks)

_wgms_window_ok(season, t0, t1, qc) =
    (lims = getproperty(qc.window_days, season); lims[1] <= Dates.value(t1 - t0) <= lims[2])

_wgms_balance_ok(season, b, qc) =
    (lims = getproperty(qc.balance_limits, season); lims[1] <= b <= lims[2])

# A seasonal record, the schema both sources share and the model step consumes.
_wgms_record(; source, record_id, glacier_id, year, region, season, latitude, longitude,
             location_source, z_lo, z_hi, begin_date, begin_date_unc, end_date, end_date_unc,
             date_source, time_system, obs, obs_unc, density, density_class) =
    (; source, record_id, glacier_id, year, region, season, geotile=wgms_geotile_id(latitude, longitude),
     latitude, longitude, location_source, z_lo, z_hi, begin_date, begin_date_unc, end_date,
     end_date_unc, date_source, time_system, obs, obs_unc, density, density_class)

"""
    wgms_seasonal_points(wgms; qc=WGMS_SEASONAL_QC) -> (; records, excluded)

The winter and summer stake and pit balances of `wgms` (from [`read_wgms`](@ref)) that pass quality
control, one row per record, and every point not used with the reason (`excluded`: `id`, `glacier_id`,
`region`, `balance_code`, `reason`).

Rules, in order, each with the reason it records:

 1. `balance_code` is winter or summer (`"annual"`, `"index"`: index points are "Not used (index
    measurement)" in the FoG documentation and carry no season);
 2. a balance and an elevation are present (`"no balance"`, `"no elevation"`);
 3. the location is the point's own, else its glacier's (`location_source`), and lies within
    `qc.max_point_distance_km` of the glacier's coordinates (`"mislocated"`);
 4. each window end has an uncertainty ≤ `qc.max_date_unc_days`. A placeholder end is replaced by the
    glacier-year's date for the same season boundary (begin and midseason for winter, midseason and end
    for summer) when that date is certain to `qc.max_date_unc_days` and lies inside the placeholder's
    interval: a point coded winter or summer was used in exactly that glacier-wide seasonal balance.
    `date_source` is `"point"`, `"glacier_year"` or `"mixed"`. Otherwise `"date uncertain"`;
 5. window length within `qc.window_days` (`"window length"`);
 6. remarks (`"remarks"`), balance plausibility (`"implausible balance"`), reported uncertainty
    (`"balance uncertainty"`), winter density within `qc.winter_density_limits` (`"implausible
    density"`);
 7. duplicates (same glacier, season, dates, location and elevation) are kept once when their balances
    agree within `qc.duplicate_tolerance` and excluded whole otherwise (`"conflicting duplicate"`;
    the extra copies of an agreeing group are `"duplicate"`).

`density_class` is `"not_reported"`; `"assumed"` (the remarks say the density was assumed or
estimated, see `qc.assumed_density_remarks`); `"single_value"` (every reported winter density of the
glacier-year is one value: one density applied to all probe depths); or `"measured_varied"`. Only
`"measured_varied"` densities are measurements at the record's own point.
"""
function wgms_seasonal_points(wgms; qc=WGMS_SEASONAL_QC)
    glacier = Dict(r.id => r for r in eachrow(wgms.glaciers))
    gyear = Dict((r.glacier_id, r.year) => r for r in eachrow(wgms.glacier_years))
    haversine = Haversine(6371.0)

    excluded = DataFrame(id=Int[], glacier_id=Int[], region=Int[], balance_code=String[], reason=String[])
    exclude!(p, reason) = push!(excluded, (p.id, p.glacier_id, glacier[p.glacier_id].region, p.balance_code, reason))

    # Seasonal boundary dates of the glacier-year, for substitution.
    function boundary(gy, f)
        (isnothing(gy) || ismissing(gy[f])) && return nothing
        u = gy[f*"_unc"]
        (ismissing(u) || u > qc.max_date_unc_days) && return nothing
        return (gy[f], u)
    end
    function resolve(d, u, gy, f)
        u <= qc.max_date_unc_days && return (d, u, false)
        s = boundary(gy, f)
        (isnothing(s) || abs(Dates.value(s[1] - d)) > u) && return nothing
        return (s[1], s[2], true)
    end

    # Winter densities of each glacier-year, to tell assumed from measured.
    densities = Dict{Tuple{Int,Int},Set{Float64}}()
    for p in eachrow(wgms.points)
        (p.balance_code == "winter" && !ismissing(p.density)) || continue
        push!(get!(densities, (p.glacier_id, p.year), Set{Float64}()), p.density)
    end

    rows = NamedTuple[]
    for p in eachrow(wgms.points)
        p.balance_code in ("winter", "summer") || (exclude!(p, p.balance_code); continue)
        ismissing(p.balance) && (exclude!(p, "no balance"); continue)
        ismissing(p.elevation) && (exclude!(p, "no elevation"); continue)
        g = glacier[p.glacier_id]
        season = Symbol(p.balance_code)

        own = !ismissing(p.latitude) && !ismissing(p.longitude)
        lat, lon = own ? (p.latitude, p.longitude) : (g.latitude, g.longitude)
        if own && haversine((lon, lat), (g.longitude, g.latitude)) > qc.max_point_distance_km
            exclude!(p, "mislocated"); continue
        end

        gy = get(gyear, (p.glacier_id, p.year), nothing)
        f0, f1 = season === :winter ? ("begin_date", "midseason_date") : ("midseason_date", "end_date")
        r0 = resolve(p.begin_date, p.begin_date_unc, gy, f0)
        r1 = resolve(p.end_date, p.end_date_unc, gy, f1)
        (isnothing(r0) || isnothing(r1)) && (exclude!(p, "date uncertain"); continue)
        date_source = (r0[3] && r1[3]) ? "glacier_year" : (r0[3] || r1[3]) ? "mixed" : "point"

        _wgms_window_ok(season, r0[1], r1[1], qc) || (exclude!(p, "window length"); continue)
        _wgms_remark_excluded(p.remarks, qc) && (exclude!(p, "remarks"); continue)
        _wgms_balance_ok(season, p.balance, qc) || (exclude!(p, "implausible balance"); continue)
        (!ismissing(p.balance_unc) && p.balance_unc > qc.max_balance_unc) &&
            (exclude!(p, "balance uncertainty"); continue)

        density_class = "not_reported"
        if season === :winter && !ismissing(p.density)
            lo, hi = qc.winter_density_limits
            lo <= p.density <= hi || (exclude!(p, "implausible density"); continue)
            density_class = _wgms_remark_matches(p.remarks, qc.assumed_density_remarks) ? "assumed" :
                            length(densities[(p.glacier_id, p.year)]) == 1 ? "single_value" : "measured_varied"
        end

        push!(rows, merge((; id=p.id), _wgms_record(; source=:point, record_id=string(p.id),
            glacier_id=p.glacier_id, year=p.year, region=g.region, season, latitude=lat, longitude=lon,
            location_source=own ? "point" : "glacier", z_lo=p.elevation, z_hi=p.elevation,
            begin_date=r0[1], begin_date_unc=r0[2], end_date=r1[1], end_date_unc=r1[2], date_source,
            time_system=isnothing(gy) ? missing : gy.time_system, obs=p.balance, obs_unc=p.balance_unc,
            density=season === :winter ? p.density : missing, density_class)))
    end
    records = DataFrame(rows)
    isempty(rows) && return (; records, excluded)

    # Duplicates are judged on the location the point itself reports; a point located by its glacier
    # duplicates another only if both were.
    key = [:glacier_id, :season, :begin_date, :end_date, :latitude, :longitude, :z_lo]
    keep = trues(nrow(records))
    for grp in groupby(records, key)
        nrow(grp) == 1 && continue
        idx = parentindices(grp)[1]
        # Agreeing duplicates keep their first record; conflicting ones are all excluded.
        agree = maximum(grp.obs) - minimum(grp.obs) <= qc.duplicate_tolerance
        for i in (agree ? idx[2:end] : idx)
            keep[i] = false
            push!(excluded, (records.id[i], records.glacier_id[i], records.region[i],
                             String(records.season[i]), agree ? "duplicate" : "conflicting duplicate"))
        end
    end
    return (; records=records[keep, :], excluded)
end

"""
    wgms_seasonal_bands(wgms; qc=WGMS_SEASONAL_QC) -> (; records, excluded)

The winter and summer elevation-band balances of `wgms` that pass quality control, one row per band
and season, with the window taken from the glacier-year in `mass_balance.csv` (winter: begin to
midseason; summer: midseason to end). `excluded` lists every band-season not used (`glacier_id`,
`year`, `lower_elevation`, `upper_elevation`, `region`, `season`, `reason`).

A band-season needs its balance (`"no balance"`), both window dates with uncertainty ≤
`qc.max_date_unc_days` and begin < midseason < end (`"date uncertain"`), a window length within
`qc.window_days` (`"window length"`), a width 0 < upper − lower ≤ `qc.max_band_width` (`"band
width"`), no excluding remark (`"remarks"`), and a plausible balance (`"implausible balance"`) with a
reported uncertainty ≤ `qc.max_balance_unc` (`"balance uncertainty"`). The band's location is its
glacier's coordinates.
"""
function wgms_seasonal_bands(wgms; qc=WGMS_SEASONAL_QC)
    glacier = Dict(r.id => r for r in eachrow(wgms.glaciers))
    gyear = Dict((r.glacier_id, r.year) => r for r in eachrow(wgms.glacier_years))

    excluded = DataFrame(glacier_id=Int[], year=Int[], lower_elevation=Float64[], upper_elevation=Float64[],
                         region=Int[], season=Symbol[], reason=String[])
    rows = NamedTuple[]
    dated(gy, f) = !ismissing(gy[f]) && !ismissing(gy[f*"_unc"]) && gy[f*"_unc"] <= qc.max_date_unc_days

    for b in eachrow(wgms.bands), season in (:winter, :summer)
        g = glacier[b.glacier_id]
        exclude!(reason) = push!(excluded, (b.glacier_id, b.year, b.lower_elevation, b.upper_elevation,
                                            g.region, season, reason))
        obs = season === :winter ? b.winter_balance : b.summer_balance
        unc = season === :winter ? b.winter_balance_unc : b.summer_balance_unc
        ismissing(obs) && (exclude!("no balance"); continue)

        gy = gyear[(b.glacier_id, b.year)]
        f0, f1 = season === :winter ? ("begin_date", "midseason_date") : ("midseason_date", "end_date")
        ok = dated(gy, "begin_date") && dated(gy, "midseason_date") && dated(gy, "end_date") &&
             gy.begin_date < gy.midseason_date < gy.end_date
        ok || (exclude!("date uncertain"); continue)
        t0, t1 = gy[f0], gy[f1]

        _wgms_window_ok(season, t0, t1, qc) || (exclude!("window length"); continue)
        0 < b.upper_elevation - b.lower_elevation <= qc.max_band_width || (exclude!("band width"); continue)
        _wgms_remark_excluded(b.remarks, qc) && (exclude!("remarks"); continue)
        _wgms_balance_ok(season, obs, qc) || (exclude!("implausible balance"); continue)
        (!ismissing(unc) && unc > qc.max_balance_unc) && (exclude!("balance uncertainty"); continue)

        push!(rows, _wgms_record(; source=:band,
            record_id="$(b.glacier_id)/$(b.year)/$(Int(b.lower_elevation))-$(Int(b.upper_elevation))",
            glacier_id=b.glacier_id, year=b.year, region=g.region, season, latitude=g.latitude,
            longitude=g.longitude, location_source="glacier", z_lo=b.lower_elevation,
            z_hi=b.upper_elevation, begin_date=t0, begin_date_unc=gy[f0*"_unc"], end_date=t1,
            end_date_unc=gy[f1*"_unc"], date_source="glacier_year", time_system=gy.time_system,
            obs, obs_unc=unc, density=missing, density_class="not_reported"))
    end
    return (; records=DataFrame(rows), excluded)
end

"""
    GEMBBandSMB(band_centers, time, smb, refreeze=nothing, height=nothing)

One geotile's GEMB band series in the form the WGMS comparison samples:

  * `band_centers` [m], ascending;
  * `time`: the instants the levels refer to, ascending;
  * `smb`: (time × band) cumulative climatic surface mass balance, precipitation − runoff +
    evaporation/condensation, in m w.e. from any reference;
  * `refreeze`: the same for GEMB's refreeze, or `nothing`;
  * `height`: the same for the surface height in m (the negative cumulative `ice_flux`, settling of
    the existing column included), or `nothing`.

A level at an arbitrary instant is interpolated linearly in time, so a window's balance over an output
period assumes a uniform flux inside that period.
"""
struct GEMBBandSMB{T<:Real}
    band_centers::Vector{Float64}
    time::Vector{DateTime}
    smb::Matrix{T}
    refreeze::Union{Nothing,Matrix{T}}
    height::Union{Nothing,Matrix{T}}
    function GEMBBandSMB{T}(band_centers, time, smb, refreeze=nothing, height=nothing) where {T}
        issorted(band_centers; lt=<=) || error("band centers are not strictly ascending: $band_centers")
        issorted(time; lt=<=) || error("GEMB times are not strictly ascending")
        expected = (length(time), length(band_centers))
        for (name, m) in (("smb", smb), ("refreeze", refreeze), ("height", height))
            isnothing(m) && continue
            size(m) == expected || throw(DimensionMismatch("$name is $(size(m)), expected (time, band) = $expected"))
            all(isfinite, m) || error("$name holds non-finite values")
        end
        return new{T}(collect(Float64, band_centers), collect(DateTime, time), smb, refreeze, height)
    end
end
GEMBBandSMB(band_centers, time, smb::AbstractMatrix{T}, refreeze=nothing, height=nothing) where {T} =
    GEMBBandSMB{T}(band_centers, time, Matrix(smb), _maybe_matrix(refreeze), _maybe_matrix(height))
_maybe_matrix(m) = isnothing(m) ? nothing : Matrix(m)

"""
    wgms_band_weights(band_centers, z_lo, z_hi; max_extrapolation_m=50) -> Vector{Float64} or nothing

Weights `w` over `band_centers` such that `w ⋅ v` is the mean over [`z_lo`, `z_hi`] of the profile
that interpolates the band values `v` linearly between centres and holds the outermost value beyond
them; for `z_lo == z_hi` it is that profile's value at the point. This is the height rule
`_gemb_band_weights` fills the hypsometry with, so a point and the altimetric volume see one profile.
The mean is uniform in height.

Returns `nothing` when the interval reaches further than `max_extrapolation_m` beyond the outermost
band centre.
"""
function wgms_band_weights(band_centers, z_lo, z_hi; max_extrapolation_m=50)
    z_lo <= z_hi || error("z_lo = $z_lo is above z_hi = $z_hi")
    c = band_centers
    n = length(c)
    n >= 1 || error("no band centers")
    (z_lo < first(c) - max_extrapolation_m || z_hi > last(c) + max_extrapolation_m) && return nothing

    w = zeros(Float64, n)
    n == 1 && (w[1] = 1.0; return w)

    # The profile value at z as weights: constant outside, linear between neighbours.
    function point!(w, z, f)
        z <= first(c) && (w[1] += f; return w)
        z >= last(c) && (w[n] += f; return w)
        j = searchsortedfirst(c, z)
        c[j] == z && (w[j] += f; return w)
        a = (z - c[j-1]) / (c[j] - c[j-1])
        w[j-1] += f * (1 - a)
        w[j] += f * a
        return w
    end
    z_lo == z_hi && return point!(w, z_lo, 1.0)

    # Exact mean of a piecewise-linear function: integrate each linear piece by the trapezoid rule, on
    # the breakpoints the interval contains.
    knots = sort(unique(vcat(z_lo, z_hi, filter(x -> z_lo < x < z_hi, c))))
    width = z_hi - z_lo
    for k in 1:length(knots)-1
        h = (knots[k+1] - knots[k]) / width
        point!(w, knots[k], h / 2)
        point!(w, knots[k+1], h / 2)
    end
    return w
end

# Level of each band column at `t` by linear interpolation in time; throws outside the series.
function _gemb_level_at(time, levels, t)
    first(time) <= t <= last(time) ||
        error("$t is outside the GEMB series $(first(time)) .. $(last(time))")
    j = searchsortedfirst(time, t)
    time[j] == t && return levels[j, :]
    a = Dates.value(t - time[j-1]) / Dates.value(time[j] - time[j-1])
    return (1 - a) .* levels[j-1, :] .+ a .* levels[j, :]
end

# Change of the weighted level over [t0, t1].
_gemb_window(time, levels, w, t0, t1) =
    sum(w .* (_gemb_level_at(time, levels, t1) .- _gemb_level_at(time, levels, t0)))

"""
    wgms_model_balances(records, model; max_extrapolation_m=50) -> (; matched, excluded)

The modelled balance of every WGMS seasonal record (from [`wgms_seasonal_points`](@ref) or
[`wgms_seasonal_bands`](@ref)) over its own window.

`model(geotile_id)` returns the geotile's [`GEMBBandSMB`](@ref), or `nothing` where there is no model.
A window runs from 00:00 of `begin_date` to 00:00 of `end_date`; the record's elevation or band is
sampled with [`wgms_band_weights`](@ref).

`matched` is `records` plus:
  * `model` [m w.e.]: the modelled balance;
  * `model_refreeze` [m w.e.]: GEMB's refreeze in the window, an upper bound on the internal
    accumulation the model counts and a glaciological survey does not (missing without refreeze);
  * `model_dh` [m]: the modelled surface height change over the window (missing without height);
  * `model_date_lo`, `model_date_hi` [m w.e.]: the modelled balance with each window end moved by its
    date uncertainty, narrowest and widest window.

`excluded` lists records with no model for their geotile (`"no model"`), outside the band range
(`"outside band range"`) or with a window outside the model series (`"outside model record"`).
"""
function wgms_model_balances(records, model; max_extrapolation_m=50)
    cache = Dict{String,Any}()
    getmodel(id) = get!(() -> model(id), cache, id)

    n = nrow(records)
    out = (model=Vector{Union{Missing,Float64}}(missing, n),
           model_refreeze=Vector{Union{Missing,Float64}}(missing, n),
           model_dh=Vector{Union{Missing,Float64}}(missing, n),
           model_date_lo=Vector{Union{Missing,Float64}}(missing, n),
           model_date_hi=Vector{Union{Missing,Float64}}(missing, n))
    reason = Vector{Union{Nothing,String}}(nothing, n)

    for (i, r) in enumerate(eachrow(records))
        m = getmodel(r.geotile)
        isnothing(m) && (reason[i] = "no model"; continue)
        w = wgms_band_weights(m.band_centers, r.z_lo, r.z_hi; max_extrapolation_m)
        isnothing(w) && (reason[i] = "outside band range"; continue)

        t0, t1 = DateTime(r.begin_date), DateTime(r.end_date)
        u0, u1 = Day(round(Int, r.begin_date_unc)), Day(round(Int, r.end_date_unc))
        if t0 - u0 < first(m.time) || t1 + u1 > last(m.time)
            reason[i] = "outside model record"; continue
        end
        out.model[i] = _gemb_window(m.time, m.smb, w, t0, t1)
        narrow = _gemb_window(m.time, m.smb, w, t0 + u0, t1 - u1)
        wide = _gemb_window(m.time, m.smb, w, t0 - u0, t1 + u1)
        out.model_date_lo[i], out.model_date_hi[i] = minmax(narrow, wide)
        isnothing(m.refreeze) || (out.model_refreeze[i] = _gemb_window(m.time, m.refreeze, w, t0, t1))
        isnothing(m.height) || (out.model_dh[i] = _gemb_window(m.time, m.height, w, t0, t1))
    end

    ok = isnothing.(reason)
    matched = records[ok, :]
    for (k, v) in pairs(out)
        matched[!, k] = v[ok]
    end
    excluded = DataFrame(source=records.source[.!ok], record_id=records.record_id[.!ok],
                         region=records.region[.!ok], season=records.season[.!ok],
                         reason=String.(reason[.!ok]))
    return (; matched, excluded)
end

"""
    gemb_tile_band_smb(path; pscale, ΔT) -> GEMBBandSMB

The band series of a GEMB_GlacierSims tile file at forcing correction (`pscale`, `ΔT`), interpolated
bilinearly across the file's forcing nodes as `gemb_dv_sample` interpolates `dv`.

The per-interval band sums of precipitation, runoff and evaporation/condensation [kg m⁻²] give the
SMB, and refreeze is read alongside; each is accumulated over the file's output periods and assigned
to the period's output time, from zero before the first. Throws if (`pscale`, `ΔT`) lies outside the
forcing grid.

`pscale` may be a vector with one value per band, each band interpolated at its own (a precipitation
gradient, [`gemb_band_pscale_multipliers`](@ref)).
"""
function gemb_tile_band_smb(path; pscale, ΔT)
    return NCDataset(path, "r"; maskingvalue=NaN) do ds
        band_centers = Float64.(ds["band_center"][:])
        time = DateTime.(ds["time"][:])
        dts = Float64.(ds["delta_temperature"][:])
        pss = Float64.(ds["precipitation_scaling"][:])
        it, ip = sortperm(dts), sortperm(pss)
        dts, pss = dts[it], pss[ip]
        (all(p -> first(pss) <= p <= last(pss), pscale) && first(dts) <= ΔT <= last(dts)) ||
            error("($pscale, $ΔT) is outside the forcing grid of $(basename(path)): pscale " *
                  "$(first(pss))..$(last(pss)), ΔT $(first(dts))..$(last(dts))")
        pscale isa AbstractVector && length(pscale) != length(band_centers) &&
            throw(DimensionMismatch("$(length(pscale)) band pscales for $(length(band_centers)) bands"))

        # Bilinear weights over the sorted node axes.
        function axis_weights(axis, x)
            j = searchsortedfirst(axis, x)
            (j <= length(axis) && axis[j] == x) && return [(j, 1.0)]
            a = (x - axis[j-1]) / (axis[j] - axis[j-1])
            return [(j - 1, 1 - a), (j, a)]
        end
        layer = if pscale isa AbstractVector
            function (name)
                A = ds[name][:, :, :, :]
                out = zeros(size(A, 1), length(band_centers))
                for b in eachindex(band_centers), (jt, wt) in axis_weights(dts, ΔT),
                    (jp, wp) in axis_weights(pss, pscale[b])
                    out[:, b] .+= (wt * wp) .* Float64.(view(A, :, b, it[jt], ip[jp]))
                end
                out
            end
        else
            nodes = [(it[jt], ip[jp], wt * wp) for (jt, wt) in axis_weights(dts, ΔT) for (jp, wp) in axis_weights(pss, pscale)]
            name -> sum(w .* Float64.(ds[name][:, :, t, p]) for (t, p, w) in nodes)
        end
        smb = cumsum(layer("precipitation") .- layer("runoff") .+ layer("evaporation_condensation"); dims=1) ./ 1000
        refreeze = cumsum(layer("refreeze"); dims=1) ./ 1000
        GEMBBandSMB(band_centers, time, smb, refreeze)
    end
end

"""
    wgms_tile_model(tile_dir, fits; pscale=:pscale, ΔT=:ΔT, physics_member=nothing, area_km2=nothing)
        -> Function

The `model` argument of [`wgms_model_balances`](@ref) for a GEMB_GlacierSims tile sweep in `tile_dir`
at the forcing corrections of `fits`, a table with an `id` column of geotile ids and the fitted
corrections in columns `pscale` and `ΔT` (a `*_gembfit.arrow` read into a `DataFrame`). Pass numbers
instead of column names to evaluate every geotile at one node, e.g. `pscale=1.0, ΔT=0.0` for the
uncalibrated model.

A geotile with no tile file of its own, or no fit, has no model: the comparison scores a geotile's own
series, not a donor's.

With `physics_member` (a GEMB physics member id) every tile file must have been run as its tile member
([`gemb_tile_provenance`](@ref)), so a score cannot pair one member's tiles with another's fits. A member
with a precipitation gradient ([`gemb_physics_member_parts`](@ref)) reads each band at
`pscale · m_b`, the multipliers the calibration applied ([`gemb_band_pscale_multipliers`](@ref)); it
needs `area_km2`, the `(geotile, height)` hypsometry the calibration weighted bands with
(`_geotile_area_km2`), since `m_b` is centred on it. Bands holding no ice in it are read at `pscale`.
"""
function wgms_tile_model(tile_dir, fits; pscale=:pscale, ΔT=:ΔT, physics_member=nothing, area_km2=nothing)
    paths = gemb_tile_paths(tile_dir)
    parts = isnothing(physics_member) ? (; tile_member=nothing, precipitation_gradient=0.0) :
            gemb_physics_member_parts(physics_member)
    isnothing(physics_member) ||
        gemb_tile_provenance(tile_dir, collect(values(paths)); physics_member=parts.tile_member)
    g = parts.precipitation_gradient
    g > 0 && isnothing(area_km2) &&
        error("physics member $(repr(physics_member)) applies a precipitation gradient centred on each " *
              "geotile's hypsometry; pass area_km2")
    fit = Dict{String,Tuple{Float64,Float64}}()
    for r in eachrow(fits)
        ismissing(r.id) && continue
        p = pscale isa Symbol ? r[pscale] : pscale
        t = ΔT isa Symbol ? r[ΔT] : ΔT
        (ismissing(p) || ismissing(t)) && continue
        fit[r.id] = (Float64(p), Float64(t))
    end
    return function (geotile_id)
        (haskey(paths, geotile_id) && haskey(fit, geotile_id)) || return nothing
        p, t = fit[geotile_id]
        iszero(g) && return gemb_tile_band_smb(paths[geotile_id]; pscale=p, ΔT=t)
        centers = NCDataset(ds -> Float64.(ds["band_center"][:]), paths[geotile_id], "r")
        w, _ = _gemb_band_weights(centers, area_km2[geotile=At(geotile_id)])
        m = gemb_band_pscale_multipliers(centers, w, g)
        return gemb_tile_band_smb(paths[geotile_id]; pscale=[w[b] > 0 ? p * m[b] : p for b in eachindex(m)],
                                  ΔT=t)
    end
end

"""
    top_mass_bulk_density(dz, density, mass) -> Float64

Bulk density [kg m⁻³] of the uppermost `mass` [kg m⁻²] of a column given surface-first layer
thicknesses `dz` [m] and densities `density` [kg m⁻³]: `mass` over the depth that holds it, the
deepest layer counted only in part. Throws if the column holds less than `mass`.
"""
function top_mass_bulk_density(dz, density, mass)
    axes(dz) == axes(density) || throw(DimensionMismatch("dz $(axes(dz)) and density $(axes(density)) differ"))
    mass > 0 || error("mass must be positive, got $mass")
    m, depth = 0.0, 0.0
    for i in eachindex(dz, density)
        layer = dz[i] * density[i]
        if m + layer >= mass
            return mass / (depth + (mass - m) / density[i])
        end
        m += layer
        depth += dz[i]
    end
    error("the column holds $m kg m-2, less than the $mass kg m-2 requested")
end

"""
    young_snow_bulk_density(dz, density, age, max_age_days) -> (; density, mass, depth)

Bulk density [kg m⁻³], mass [kg m⁻²] and depth [m] of the surface layers whose mass-weighted age is at
most `max_age_days`, taken from the surface down to the first older layer: the snow accumulated since
`max_age_days` ago. Layers merge old and new snow, so the boundary is as sharp as the model's layering.
`density` is `NaN` when the surface layer is already older.
"""
function young_snow_bulk_density(dz, density, age, max_age_days)
    axes(dz) == axes(density) == axes(age) ||
        throw(DimensionMismatch("dz, density and age must share indices"))
    m, depth = 0.0, 0.0
    for i in eachindex(dz, density, age)
        age[i] <= max_age_days || break
        m += dz[i] * density[i]
        depth += dz[i]
    end
    return (; density=depth > 0 ? m / depth : NaN, mass=m, depth)
end

# Pearson r, NaN for fewer than three pairs or zero variance.
function _wgms_r(x, y)
    length(x) < 3 && return NaN
    sx, sy = std(x), std(y)
    (sx == 0 || sy == 0) && return NaN
    return cor(x, y)
end

"""
    wgms_seasonal_statistics(matched; by=[:source, :region, :season], elevation_bin_m=nothing) -> DataFrame

Agreement of modelled and observed seasonal balances per group of `by`: `n`, `n_glaciers`,
`n_glacier_years`, `obs_mean`, `model_mean`, `bias` (mean model − obs), `rmse`, `r` (Pearson).
With `elevation_bin_m`, records are also grouped by the bin of their mean elevation, `z_bin` being the
bin's lower edge.
"""
function wgms_seasonal_statistics(matched; by=[:source, :region, :season], elevation_bin_m=nothing)
    df = copy(matched)
    groupcols = copy(by)
    if !isnothing(elevation_bin_m)
        df.z_bin = floor.((df.z_lo .+ df.z_hi) ./ 2 ./ elevation_bin_m) .* elevation_bin_m
        push!(groupcols, :z_bin)
    end
    return DataFrames.combine(groupby(df, groupcols; sort=true),
        nrow => :n,
        :glacier_id => (x -> length(unique(x))) => :n_glaciers,
        [:glacier_id, :year] => ((g, y) -> length(unique(zip(g, y)))) => :n_glacier_years,
        :obs => mean => :obs_mean,
        :model => mean => :model_mean,
        [:model, :obs] => ((m, o) -> mean(m .- o)) => :bias,
        [:model, :obs] => ((m, o) -> sqrt(mean((m .- o) .^ 2))) => :rmse,
        [:model, :obs] => ((m, o) -> _wgms_r(m, o)) => :r)
end

"""
    wgms_vertical_gradients(matched; min_n=3, min_span_m=200) -> DataFrame

The least-squares balance gradient with elevation [m w.e. per 100 m], observed and modelled at the same
elevations, for every source × glacier × year × season with at least `min_n` records spanning at least
`min_span_m`. Columns `source`, `region`, `glacier_id`, `year`, `season`, `n`, `z_min`, `z_max`,
`span_m`, `obs_gradient`, `obs_gradient_se` (the standard error of the observed slope from its
residuals, `NaN` for `n == 2`), `model_gradient`. Elevations are record mid-points.
"""
function wgms_vertical_gradients(matched; min_n=3, min_span_m=200)
    slope(z, b) = (zc = z .- mean(z); sum(zc .* (b .- mean(b))) / sum(zc .^ 2) * 100)
    function slope_se(z, b)
        length(z) <= 2 && return NaN
        zc = z .- mean(z)
        β = sum(zc .* (b .- mean(b))) / sum(zc .^ 2)
        r = b .- mean(b) .- β .* zc
        return sqrt(sum(r .^ 2) / (length(z) - 2) / sum(zc .^ 2)) * 100
    end
    df = transform(matched, [:z_lo, :z_hi] => ((a, b) -> (a .+ b) ./ 2) => :z)
    g = DataFrames.combine(groupby(df, [:source, :region, :glacier_id, :year, :season]),
        nrow => :n, :z => minimum => :z_min, :z => maximum => :z_max,
        :z => (z -> maximum(z) - minimum(z)) => :span_m,
        [:z, :obs] => slope => :obs_gradient, [:z, :obs] => slope_se => :obs_gradient_se,
        [:z, :model] => slope => :model_gradient)
    return g[(g.n .>= min_n) .& (g.span_m .>= min_span_m), :]
end

"""
    wgms_anomaly_correlation(matched; min_years=5, by=[:source, :region, :season]) -> DataFrame

Temporal skill: the Pearson r of interannual anomalies at fixed sites, pooled per group of `by`. A site
is a band (`record_id` without its year), or a point location with its elevation rounded to 10 m; it
needs `min_years` distinct years, and an anomaly is relative to the site's own mean. Columns `n_sites`,
`n`, `r_anomaly`, `obs_anomaly_sd`, `model_anomaly_sd`.
"""
function wgms_anomaly_correlation(matched; min_years=5, by=[:source, :region, :season])
    df = copy(matched)
    df.site = [r.source === :band ? "$(r.glacier_id)/$(Int(r.z_lo))-$(Int(r.z_hi))" :
               "$(r.glacier_id)/$(round(r.latitude; digits=3))/$(round(r.longitude; digits=3))/$(round(Int, r.z_lo / 10))"
               for r in eachrow(df)]
    sites = DataFrames.combine(groupby(df, [:site, :season]), :year => (y -> length(unique(y))) => :n_years)
    good = Set((r.site, r.season) for r in eachrow(sites) if r.n_years >= min_years)
    df = df[[(s, x) in good for (s, x) in zip(df.site, df.season)], :]
    isempty(df) && return DataFrame()
    df = transform(groupby(df, [:site, :season]), :obs => (x -> x .- mean(x)) => :obs_anom,
                   :model => (x -> x .- mean(x)) => :model_anom)
    return DataFrames.combine(groupby(df, by; sort=true),
        :site => (s -> length(unique(s))) => :n_sites, nrow => :n,
        [:model_anom, :obs_anom] => ((m, o) -> _wgms_r(m, o)) => :r_anomaly,
        :obs_anom => std => :obs_anomaly_sd, :model_anom => std => :model_anomaly_sd)
end

"""
    wgms_seasonal_ratio(matched; by=[:source, :region]) -> DataFrame

The winter-to-summer ratio mean(b_w) / mean(−b_s), observed and modelled, over the site-years that
have both seasons (a band, or a point location and elevation, in one glacier-year). Columns `n_pairs`,
`obs_ratio`, `model_ratio`.
"""
function wgms_seasonal_ratio(matched; by=[:source, :region])
    df = copy(matched)
    df.site = [r.source === :band ? r.record_id :
               "$(r.glacier_id)/$(r.year)/$(round(r.latitude; digits=3))/$(round(r.longitude; digits=3))/$(round(Int, r.z_lo / 10))"
               for r in eachrow(df)]
    w = df[df.season .=== :winter, vcat(by, [:site, :obs, :model])]
    s = df[df.season .=== :summer, vcat(by, [:site, :obs, :model])]
    j = innerjoin(w, s; on=vcat(by, [:site]), renamecols="_w" => "_s")
    isempty(j) && return DataFrame()
    return DataFrames.combine(groupby(j, by; sort=true), nrow => :n_pairs,
        [:obs_w, :obs_s] => ((a, b) -> mean(a) / mean(-b)) => :obs_ratio,
        [:model_w, :model_s] => ((a, b) -> mean(a) / mean(-b)) => :model_ratio)
end

# Per source × region × season medians of the glacier-year vertical gradients.
_wgms_gradient_summary(gradients) = isempty(gradients) ? DataFrame() :
    DataFrames.combine(groupby(gradients, [:source, :region, :season]; sort=true), nrow => :n,
            :obs_gradient => median => :obs_median, :model_gradient => median => :model_median,
            [:model_gradient, :obs_gradient] => ((m, o) -> median(m .- o)) => :median_difference,
            [:model_gradient, :obs_gradient] => ((m, o) -> median(abs.(m .- o))) => :median_abs_difference)

# The tables a scorecard reads, recomputed from `matched` (one label), which restricted or paired
# records need: `wgms_validation`'s own tables describe all of its records.
function _wgms_rescore(matched)
    :label in propertynames(matched) || error("matched has no label column")
    nrow(matched) == 0 && return (; matched, summary=DataFrame(), gradient_summary=DataFrame(),
                                  anomaly=DataFrame(), ratio=DataFrame())
    label = only(unique(matched.label))
    labelled(df) = ncol(df) == 0 ? df : insertcols!(df, 1, :label => fill(label, nrow(df)))
    return (; matched,
            summary=labelled(wgms_seasonal_statistics(matched)),
            gradient_summary=labelled(_wgms_gradient_summary(wgms_vertical_gradients(matched))),
            anomaly=labelled(wgms_anomaly_correlation(matched)),
            ratio=labelled(wgms_seasonal_ratio(matched)))
end

"""
    WGMS_PROVENANCE

The values of a matched record's `provenance`, most to least anchored in measurement at the record's own
elevation: `:stake` (a point record), `:band_in_stake_range` (a band whose centre lies within
[`wgms_record_covariates`](@ref)'s tolerance of the elevations of the same glacier-year-season's
stakes), `:band_beyond_stakes` (stakes exist but the band centre is further outside their range) and
`:band_no_stakes` (no stake of that glacier-year-season).
"""
const WGMS_PROVENANCE = (:stake, :band_in_stake_range, :band_beyond_stakes, :band_no_stakes)

"""
    WGMS_ANCHORED_PROVENANCE

The provenance classes whose elevation is measured or bracketed by stakes: `(:stake,
:band_in_stake_range)`.
"""
const WGMS_ANCHORED_PROVENANCE = (:stake, :band_in_stake_range)

"""
    wgms_record_covariates(wgms, records; stakes=:database, tolerance_m=50, qc=WGMS_SEASONAL_QC)
        -> DataFrame

`records` (seasonal point and band records with `source`, `glacier_id`, `year`, `season`, `z_lo` and
`z_hi`) with three covariates that say how far a record is from the glacier it was measured on:

  * `glacier_area_km2`: the median over years of `mass_balance.area` of the record's glacier, in km²
    (`missing` where no year reports one; [`read_wgms`](@ref)'s `area_km2`);
  * `provenance` (see [`WGMS_PROVENANCE`](@ref)): `:stake` for a point. For a band, the band centre
    `(z_lo + z_hi) / 2` against the lowest and highest stake of the same glacier, year and season, each
    side widened by `tolerance_m`;
  * `stake_span_fraction`: for a point, the elevation span of the stakes of its glacier-year-season over
    the elevation range of the glacier's bands (lowest `lower_elevation` to highest `upper_elevation`
    over every year of `mass_balance_band.csv`), `missing` where the glacier has no band range. A
    glacier-year with one stake has fraction 0, and the fraction can exceed 1 where the stakes reach
    beyond the surveyed bands. `missing` for a band.

A stake is a `mass_balance_point` record coded winter or summer with a balance and an elevation.
`stakes=:database` (default) counts every such record of the database, including those
[`wgms_seasonal_points`](@ref) excludes. Quality control judges whether a point is a fair comparison
target (its dates, window, remarks and uncertainty); the glacier's band balance was computed from the
whole stake network that season, so a point excluded from scoring still anchors the band. A point
without an elevation cannot anchor any elevation, so it does not count. `stakes=:qc` counts only the
points that pass `qc`; it gives the narrower, scored-points-only reading.
"""
function wgms_record_covariates(wgms, records; stakes=:database, tolerance_m=50, qc=WGMS_SEASONAL_QC)
    stakes in (:database, :qc) || error("stakes must be :database or :qc, got $(repr(stakes))")
    tolerance_m >= 0 || error("tolerance_m must be non-negative, got $tolerance_m")
    :area_km2 in propertynames(wgms.glaciers) ||
        error("wgms.glaciers has no area_km2 column; read the data with read_wgms")
    ncol(records) == 0 && return records
    absent = setdiff([:source, :glacier_id, :year, :season, :z_lo, :z_hi], propertynames(records))
    isempty(absent) || error("records lack column(s) $(join(absent, ", "))")

    # Lowest and highest stake elevation of each glacier-year-season.
    extremes = Dict{Tuple{Int,Int,Symbol},Tuple{Float64,Float64}}()
    function add_stake!(glacier_id, year, season, z)
        k = (glacier_id, year, season)
        lo, hi = get(extremes, k, (Inf, -Inf))
        extremes[k] = (min(lo, z), max(hi, z))
    end
    if stakes === :database
        for p in eachrow(wgms.points)
            (p.balance_code in ("winter", "summer") && !ismissing(p.balance) && !ismissing(p.elevation)) ||
                continue
            add_stake!(p.glacier_id, p.year, Symbol(p.balance_code), p.elevation)
        end
    else
        for r in eachrow(wgms_seasonal_points(wgms; qc).records)
            add_stake!(r.glacier_id, r.year, r.season, r.z_lo)
        end
    end

    # Elevation range of each glacier's bands over every year.
    band_range = Dict{Int,Tuple{Float64,Float64}}()
    for b in eachrow(wgms.bands)
        lo, hi = get(band_range, b.glacier_id, (Inf, -Inf))
        band_range[b.glacier_id] = (min(lo, b.lower_elevation), max(hi, b.upper_elevation))
    end
    area = Dict(zip(wgms.glaciers.id, wgms.glaciers.area_km2))

    n = nrow(records)
    glacier_area_km2 = Vector{Union{Missing,Float64}}(missing, n)
    provenance = Vector{Symbol}(undef, n)
    stake_span_fraction = Vector{Union{Missing,Float64}}(missing, n)
    for (i, r) in zip(eachindex(provenance), eachrow(records))
        haskey(area, r.glacier_id) || error("record $(r.record_id): glacier $(r.glacier_id) is not in wgms.glaciers")
        glacier_area_km2[i] = area[r.glacier_id]
        e = get(extremes, (r.glacier_id, r.year, r.season), nothing)
        if r.source === :point
            isnothing(e) && error("point record $(r.record_id) has no stake under stakes=$(repr(stakes)); " *
                                  "the records were not drawn from wgms")
            provenance[i] = :stake
            z_bands = get(band_range, r.glacier_id, nothing)
            if !isnothing(z_bands) && z_bands[2] > z_bands[1]
                stake_span_fraction[i] = (e[2] - e[1]) / (z_bands[2] - z_bands[1])
            end
        elseif r.source === :band
            centre = (r.z_lo + r.z_hi) / 2
            provenance[i] = isnothing(e) ? :band_no_stakes :
                            (e[1] - tolerance_m <= centre <= e[2] + tolerance_m) ? :band_in_stake_range :
                            :band_beyond_stakes
        else
            error("record $(r.record_id): source must be :point or :band, got $(repr(r.source))")
        end
    end
    out = copy(records)
    out.glacier_area_km2 = glacier_area_km2
    out.provenance = provenance
    out.stake_span_fraction = stake_span_fraction
    return out
end

"""
    wgms_size_class(area_km2; edges_km2=(1, 10)) -> String

The size class of a glacier of `area_km2`: `"<1"`, `"1-10"` or `">=10"` for the default edges, one class
per interval of the ascending `edges_km2`, `"missing"` where the area is unknown. An area on an edge
belongs to the class above it.
"""
function wgms_size_class(area_km2; edges_km2=(1, 10))
    ismissing(area_km2) && return "missing"
    e = collect(edges_km2)
    return _wgms_size_labels(e)[searchsortedlast(e, area_km2)+1]
end

# Class labels in ascending order of area.
function _wgms_size_labels(edges)
    (!isempty(edges) && issorted(edges; lt=<=)) ||
        error("edges_km2 must be non-empty and strictly ascending, got $(Tuple(edges))")
    return vcat("<$(first(edges))", ["$a-$b" for (a, b) in zip(edges, Iterators.drop(edges, 1))],
                ">=$(last(edges))")
end

"""
    wgms_validation(model; wgms, label="", qc=WGMS_SEASONAL_QC, elevation_bin_m=250,
                    max_extrapolation_m=50, stakes=:database, stake_tolerance_m=50) -> NamedTuple

Score one GEMB configuration against the WGMS seasonal point and band balances: quality control
([`wgms_seasonal_points`](@ref), [`wgms_seasonal_bands`](@ref)), matching ([`wgms_model_balances`](@ref))
and statistics, points and bands kept apart.

`model(geotile_id)` returns a [`GEMBBandSMB`](@ref) or `nothing` (see [`wgms_tile_model`](@ref) for a
tile sweep at fitted corrections); `wgms` is [`read_wgms`](@ref)'s result, read once and shared across
configurations. Every table carries a `label` column, so the results of several configurations
`vcat` into one comparison.

`matched` carries the covariates of [`wgms_record_covariates`](@ref) (`glacier_area_km2`, `provenance`,
`stake_span_fraction`), computed with `stakes` and `stake_tolerance_m`; [`wgms_scorecard_strata`](@ref),
[`wgms_scorecard`](@ref)`(v; provenance)` and [`wgms_member_difference`](@ref) score by them.

Returns `(; matched, excluded_qc, excluded_model, summary, by_elevation, gradients, gradient_summary,
anomaly, ratio)`.
"""
function wgms_validation(model; wgms, label="", qc=WGMS_SEASONAL_QC, elevation_bin_m=250,
                         max_extrapolation_m=50, stakes=:database, stake_tolerance_m=50)
    pts, bds = wgms_seasonal_points(wgms; qc), wgms_seasonal_bands(wgms; qc)
    common = [:source, :record_id, :glacier_id, :year, :region, :season, :geotile, :latitude, :longitude,
              :location_source, :z_lo, :z_hi, :begin_date, :begin_date_unc, :end_date, :end_date_unc,
              :date_source, :time_system, :obs, :obs_unc, :density, :density_class]
    records = vcat(isempty(pts.records) ? DataFrame() : pts.records[:, common],
                   isempty(bds.records) ? DataFrame() : bds.records[:, common]; cols=:union)
    records = wgms_record_covariates(wgms, records; stakes, tolerance_m=stake_tolerance_m, qc)
    m = wgms_model_balances(records, model; max_extrapolation_m)
    matched = m.matched

    counts(df, source) = insertcols!(DataFrames.combine(groupby(df, [:region, :reason]), nrow => :n), 1,
                                     :source => fill(source, length(unique(zip(df.region, df.reason)))))
    excluded_qc = vcat(counts(pts.excluded, :point), counts(bds.excluded, :band))
    excluded_model = DataFrames.combine(groupby(m.excluded, [:source, :region, :season, :reason]), nrow => :n)

    labelled(df) = ncol(df) == 0 ? df : insertcols!(df, 1, :label => fill(label, nrow(df)))
    if isempty(matched)
        return (; matched=labelled(matched), excluded_qc=labelled(excluded_qc),
                excluded_model=labelled(excluded_model), summary=DataFrame(), by_elevation=DataFrame(),
                gradients=DataFrame(), gradient_summary=DataFrame(), anomaly=DataFrame(), ratio=DataFrame())
    end
    gradients = wgms_vertical_gradients(matched)
    gradient_summary = _wgms_gradient_summary(gradients)
    return (; matched=labelled(matched), excluded_qc=labelled(excluded_qc),
            excluded_model=labelled(excluded_model),
            summary=labelled(wgms_seasonal_statistics(matched)),
            by_elevation=labelled(wgms_seasonal_statistics(matched; elevation_bin_m)),
            gradients=labelled(gradients), gradient_summary=labelled(gradient_summary),
            anomaly=labelled(wgms_anomaly_correlation(matched)),
            ratio=labelled(wgms_seasonal_ratio(matched)))
end

"""
    wgms_scorecard(v; min_n=30, provenance=nothing) -> DataFrame

One row per source × season summarizing a [`wgms_validation`](@ref) result `v` with every region
weighted equally, so configurations and ensemble members can be ranked on one table. Only regions
with at least `min_n` matched records count. Columns:

  * `label`, `restriction`, `source`, `season`, `n_regions`, `n`, `n_glaciers`, `n_glacier_years` (the
    counts are of the records in the regions that count);
  * `mean_abs_bias`, `mean_rmse` [m w.e.]: means over regions of |bias| and RMSE;
  * `median_r`: median over regions of the Pearson r of the records;
  * `median_r_anomaly`: median over regions of the interannual-anomaly r at fixed sites;
  * `median_gradient_difference` [m w.e. per 100 m]: median over regions of the median model − obs
    vertical gradient of a glacier-year ([`wgms_vertical_gradients`](@ref));
  * `median_abs_gradient_difference` [m w.e. per 100 m]: the same for |model − obs|, the gradient
    error a configuration is ranked on;
  * `median_gradient_ratio`: median over regions of the median modelled over the median observed
    gradient (1 is a correct profile, 0 a flat one; `NaN` where the observed median is 0);
  * `median_abs_log_ratio_error`: median over regions of |log(model / obs winter-to-summer ratio)|,
    the same for both seasons of a source.

A statistic with no region to draw on is `NaN`.

`provenance`, a collection of [`WGMS_PROVENANCE`](@ref) classes (for example
[`WGMS_ANCHORED_PROVENANCE`](@ref)), restricts the scorecard to the records of those classes, recomputed
from `v.matched`; `restriction` is `"all"` without it and the classes joined by `+` with it.
"""
function wgms_scorecard(v; min_n=30, provenance=nothing)
    restriction = "all"
    if !isnothing(provenance)
        classes = Tuple(provenance)
        unknown = setdiff(classes, WGMS_PROVENANCE)
        isempty(unknown) || error("unknown provenance class(es) $(join(repr.(unknown), ", ")); " *
                                  "expected members of $WGMS_PROVENANCE")
        (:matched in propertynames(v) && :provenance in propertynames(v.matched)) ||
            error("v.matched has no provenance column; score it with wgms_validation")
        v = _wgms_rescore(v.matched[in.(v.matched.provenance, Ref(classes)), :])
        restriction = join(string.(classes), "+")
    end
    isempty(v.summary) && return DataFrame()
    label = first(v.summary.label)
    s = v.summary[v.summary.n .>= min_n, :]
    med(x) = isempty(x) ? NaN : median(x)
    rows = NamedTuple[]
    for (source, season) in sort(unique(zip(s.source, s.season)))
        k = (s.source .=== source) .& (s.season .=== season)
        regions = Set(s.region[k])
        pick(df, col) = (isempty(df) || !(col in propertynames(df))) ? Float64[] :
            Float64[r[col] for r in eachrow(df) if r.source === source && r.season === season &&
                    r.region in regions && isfinite(r[col])]
        gradient_ratio(df) = (isempty(df) || !(:obs_median in propertynames(df))) ? Float64[] :
            Float64[r.model_median / r.obs_median for r in eachrow(df)
                    if r.source === source && r.season === season && r.region in regions && r.obs_median != 0]
        ratio = isempty(v.ratio) ? Float64[] :
            Float64[abs(log(r.model_ratio / r.obs_ratio)) for r in eachrow(v.ratio)
                    if r.source === source && r.region in regions && r.model_ratio > 0 && r.obs_ratio > 0]
        push!(rows, (; label, restriction, source, season, n_regions=length(regions), n=sum(s.n[k]),
                     n_glaciers=sum(s.n_glaciers[k]), n_glacier_years=sum(s.n_glacier_years[k]),
                     mean_abs_bias=mean(abs.(s.bias[k])), mean_rmse=mean(s.rmse[k]),
                     median_r=med(filter(isfinite, s.r[k])),
                     median_r_anomaly=med(pick(v.anomaly, :r_anomaly)),
                     median_gradient_difference=med(pick(v.gradient_summary, :median_difference)),
                     median_abs_gradient_difference=med(pick(v.gradient_summary, :median_abs_difference)),
                     median_gradient_ratio=med(gradient_ratio(v.gradient_summary)),
                     median_abs_log_ratio_error=med(ratio)))
    end
    return DataFrame(rows)
end

"""
    WGMS_RANK_METRICS

The scorecard columns [`wgms_rank`](@ref) ranks on, each with the direction that is better:
the regional |bias| and the vertical-gradient error (lower), and the interannual-anomaly correlation
(higher). Vertical gradients are a primary metric because the calibration constrains only the
area-integrated volume, so a profile error is what it cannot correct.
"""
const WGMS_RANK_METRICS = (mean_abs_bias=:lower, median_abs_gradient_difference=:lower,
                           median_r_anomaly=:higher)

"""
    wgms_rank(cards; metrics=WGMS_RANK_METRICS) -> DataFrame

Rank the configurations in `cards`, the `vcat` of several [`wgms_scorecard`](@ref)s, within each source
× season on each of `metrics` (1 is best, ties share the lower rank, `NaN` ranks last). Returns `cards`
with a `rank_<metric>` column per metric and `mean_rank`, sorted by source, season and `mean_rank`.

Cards of several `restriction`s (scorecards with different `provenance`) are ranked within each
restriction, so one table reports the ranking for every restriction; see [`wgms_rank_restricted`](@ref)
for cards built on common records.
"""
function wgms_rank(cards; metrics=WGMS_RANK_METRICS)
    out = copy(cards)
    groups = :restriction in propertynames(out) ? [:restriction, :source, :season] : [:source, :season]
    for (metric, better) in pairs(metrics)
        better in (:lower, :higher) || error("rank direction for $metric must be :lower or :higher, got $better")
        out[!, "rank_$(metric)"] = zeros(Int, nrow(out))
    end
    for g in groupby(out, groups)
        for (metric, better) in pairs(metrics)
            x = Float64.(g[!, metric])
            key = [isnan(v) ? Inf : (better === :lower ? v : -v) for v in x]
            g[!, "rank_$(metric)"] = [1 + count(<(k), key) for k in key]
        end
    end
    out.mean_rank = vec(mean(Matrix(out[:, ["rank_$(m)" for m in keys(metrics)]]); dims=2))
    return sort!(out, [groups; :mean_rank])
end

"""
    WGMS_PROVENANCE_RESTRICTIONS

The record sets [`wgms_rank_restricted`](@ref) ranks by default: every record (`nothing`) and the
records anchored in stake measurements, [`WGMS_ANCHORED_PROVENANCE`](@ref).
"""
const WGMS_PROVENANCE_RESTRICTIONS = (nothing, WGMS_ANCHORED_PROVENANCE)

# The (source, record_id, season) key of every record of a `matched` table, which must be unique;
# `name` labels the table in the error.
function _wgms_record_keys(matched, name)
    k = collect(zip(matched.source, matched.record_id, matched.season))
    length(unique(k)) == length(k) || error("$name holds duplicate (source, record_id, season) keys")
    return k
end

"""
    wgms_rank_restricted(vs; restrictions=WGMS_PROVENANCE_RESTRICTIONS, min_n=30,
                         metrics=WGMS_RANK_METRICS) -> DataFrame

[`wgms_rank`](@ref) of the [`wgms_validation`](@ref) results `vs` (one per configuration, distinct
labels) under each of `restrictions`, a collection of `nothing` (every record) or provenance classes
([`wgms_scorecard`](@ref)`(v; provenance)`).

Every configuration is scored on the records all of them matched, so a configuration whose sweep omits
a tile does not rank on a different sample. The result is the ranked cards of all restrictions,
`restriction` telling them apart.
"""
function wgms_rank_restricted(vs; restrictions=WGMS_PROVENANCE_RESTRICTIONS, min_n=30,
                              metrics=WGMS_RANK_METRICS)
    vs = collect(vs)
    length(vs) >= 2 || error("ranking needs at least two configurations, got $(length(vs))")
    labels = [only(unique(v.matched.label)) for v in vs]
    length(unique(labels)) == length(labels) || error("configuration labels are not distinct: $labels")
    record_keys = [_wgms_record_keys(v.matched, "matched of $l") for (v, l) in zip(vs, labels)]
    common = Set(reduce(intersect, record_keys))
    isempty(common) && error("the configurations share no matched record")
    cards = DataFrame[]
    for (v, k) in zip(vs, record_keys)
        shared = v.matched[in.(k, Ref(common)), :]
        scored = _wgms_rescore(shared)
        for r in restrictions
            push!(cards, wgms_scorecard(scored; min_n, provenance=r))
        end
    end
    return wgms_rank(reduce(vcat, cards); metrics)
end

# The stratum of every record of `m` under `kind`, and the strata in reporting order.
function _wgms_strata(m, kind, size_edges_km2)
    if kind === :all
        return fill("all", nrow(m)), ["all"]
    elseif kind === :size_class
        return wgms_size_class.(m.glacier_area_km2; edges_km2=size_edges_km2),
               [_wgms_size_labels(collect(size_edges_km2)); "missing"]
    elseif kind === :provenance
        return string.(m.provenance), string.(collect(WGMS_PROVENANCE))
    elseif kind === :region
        return string.(m.region), string.(sort(unique(m.region)))
    end
    error("unknown stratification $(repr(kind)); expected :all, :size_class, :provenance or :region")
end

"""
    wgms_scorecard_strata(v; by=(:all, :size_class, :provenance, :region), min_n=30,
                          size_edges_km2=(1, 10)) -> DataFrame

The metrics of [`wgms_scorecard`](@ref) for subsets of the records of a [`wgms_validation`](@ref)
result `v`, one row per stratum × source × season:

  * `:all`: every record, the scorecard's own rows;
  * `:size_class`: the glacier's [`wgms_size_class`](@ref) (`size_edges_km2`);
  * `:provenance`: the record's class in [`WGMS_PROVENANCE`](@ref);
  * `:region`: the GTN-G region.

Each row carries `stratum_by` and `stratum`, then `n`, `n_glaciers` and `n_glacier_years`, the counts of
all records of the stratum; `n_scored`, the records in regions with at least `min_n` of them, which are
the records the metrics are computed from; `n_regions`; and the scorecard's metric columns (`NaN` where
no region has `min_n` records). Metrics weight regions equally within a stratum, so a stratum's rows
are comparable with the pooled rows and with each other. A stratum with no record is absent.
"""
function wgms_scorecard_strata(v; by=(:all, :size_class, :provenance, :region), min_n=30,
                               size_edges_km2=(1, 10))
    m = v.matched
    needed = [:glacier_area_km2, :provenance, :region, :glacier_id, :year, :source, :season]
    absent = setdiff(needed, propertynames(m))
    isempty(absent) || error("v.matched lacks column(s) $(join(absent, ", ")); score it with wgms_validation")
    blocks = DataFrame[]
    for kind in by
        class, levels = _wgms_strata(m, kind, size_edges_km2)
        for level in levels
            sub = m[class .== level, :]
            nrow(sub) == 0 && continue
            counts = DataFrames.combine(groupby(sub, [:source, :season]), nrow => :n,
                :glacier_id => (x -> length(unique(x))) => :n_glaciers,
                [:glacier_id, :year] => ((g, y) -> length(unique(zip(g, y)))) => :n_glacier_years)
            card = wgms_scorecard(_wgms_rescore(sub); min_n)
            block = isempty(card) ? counts :
                leftjoin(counts, select(card, Not(:label, :restriction, :n, :n_glaciers, :n_glacier_years),
                                        :n => :n_scored); on=[:source, :season])
            insertcols!(block, 1, :label => fill(only(unique(sub.label)), nrow(block)),
                        :stratum_by => fill(kind, nrow(block)), :stratum => fill(level, nrow(block)))
            push!(blocks, sort!(block, [:source, :season]))
        end
    end
    isempty(blocks) && return DataFrame()
    out = reduce(vcat, blocks; cols=:union)
    ids = [:label, :stratum_by, :stratum, :source, :season, :n, :n_glaciers, :n_glacier_years]
    for c in (:n_scored, :n_regions)
        out[!, c] = c in propertynames(out) ? coalesce.(out[!, c], 0) : zeros(Int, nrow(out))
    end
    for c in setdiff(propertynames(out), ids, [:n_scored, :n_regions])
        out[!, c] = coalesce.(out[!, c], NaN)
    end
    return select(out, ids, :n_scored, :n_regions, :)
end

"""
    wgms_glacier_scatter(matched; by=[:source, :region, :season], elevation_bin_m=250) -> DataFrame

The scatter between glaciers of the observed balance at the same elevation: per group of `by`, the
pooled standard deviation `glacier_scatter` [m w.e.] of the glacier means of `obs` within each
`elevation_bin_m` bin, using the bins that hold at least two glaciers (`n_dof`, the pooled degrees of
freedom, is 0 and the scatter `NaN` where none does), and `n_glaciers`. Each glacier counts once per bin,
whatever its number of records or years, so a densely sampled glacier does not set the scatter.
"""
function wgms_glacier_scatter(matched; by=[:source, :region, :season], elevation_bin_m=250)
    elevation_bin_m > 0 || error("elevation_bin_m must be positive, got $elevation_bin_m")
    df = transform(matched, [:z_lo, :z_hi] => ((a, b) -> floor.((a .+ b) ./ 2 ./ elevation_bin_m)) => :z_bin)
    return DataFrames.combine(groupby(df, by; sort=true)) do g
        means = DataFrames.combine(groupby(g, [:z_bin, :glacier_id]), :obs => mean => :x)
        ss, dof = 0.0, 0
        for b in groupby(means, :z_bin)
            nrow(b) >= 2 || continue
            ss += sum(abs2, b.x .- mean(b.x))
            dof += nrow(b) - 1
        end
        (; n_glaciers=length(unique(g.glacier_id)), n_dof=dof, glacier_scatter=dof == 0 ? NaN : sqrt(ss / dof))
    end
end

# Records both members matched, with each member's model − obs error. Throws when the members are not
# scored on one set of observations.
function _wgms_pair_members(a, b)
    need = [:source, :record_id, :season, :glacier_id, :year, :region, :z_lo, :z_hi, :obs, :model,
            :glacier_area_km2, :provenance, :label]
    for (name, t) in (("v_ref", a), ("v_member", b))
        absent = setdiff(need, propertynames(t))
        isempty(absent) || error("$name.matched lacks column(s) $(join(absent, ", ")); score it with wgms_validation")
    end
    ka, kb = _wgms_record_keys(a, "v_ref.matched"), _wgms_record_keys(b, "v_member.matched")
    index_a = Dict(k => i for (i, k) in zip(eachindex(ka), ka))
    rows_a, rows_b = Int[], Int[]
    for (j, k) in zip(eachindex(kb), kb)
        i = get(index_a, k, nothing)
        isnothing(i) || (push!(rows_a, i); push!(rows_b, j))
    end
    isempty(rows_a) && error("the two members share no matched record")
    A, B = a[rows_a, :], b[rows_b, :]
    all(isapprox.(A.obs, B.obs; atol=1e-9)) ||
        error("the members disagree on the observed balance of a shared record; they were not scored " *
              "against the same WGMS data")
    (A.glacier_id == B.glacier_id && isequal(A.glacier_area_km2, B.glacier_area_km2) &&
     A.provenance == B.provenance) ||
        error("the members disagree on the glacier or covariates of a shared record")
    out = A[:, setdiff(need, [:model, :label])]
    out.error_ref = A.model .- A.obs
    out.error_member = B.model .- B.obs
    return out, only(unique(A.label)), only(unique(B.label))
end

# Per-region sums of n, e_ref, e_member, e_ref², e_member² (rows) summarized as the region mean of
# |bias| and RMSE of each member: (|bias| ref, |bias| member, RMSE ref, RMSE member). A region with no
# record drops out of the mean.
function _wgms_error_means(S)
    out = zeros(4)
    k = 0
    for r in axes(S, 1)
        n = S[r, 1]
        n > 0 || continue
        out .+= (abs(S[r, 2]) / n, abs(S[r, 3]) / n, sqrt(S[r, 4] / n), sqrt(S[r, 5] / n))
        k += 1
    end
    return out ./ k
end

"""
    wgms_member_difference(v_ref, v_member; n_boot=1000, rng=Random.Xoshiro(0), level=0.95, min_n=30,
                           min_glaciers=5, representativeness_floor=:glacier_scatter,
                           by=(:all, :region, :size_class, :provenance), size_edges_km2=(1, 10),
                           elevation_bin_m=250) -> DataFrame

Whether two GEMB members differ against WGMS once the glaciers, not the records, are the sample.

The two [`wgms_validation`](@ref) results are paired on the records both matched (the same `source`,
`record_id` and `season`; they must carry the same observations). For each source × season and each
stratum of `by` (see [`wgms_scorecard_strata`](@ref); `:region` rows are the region × season × source
cells) the paired difference `member − ref` of |bias| and of RMSE is computed with regions weighted
equally, over the regions with at least `min_n` paired records in the stratum. Negative means the member
is closer to the observations.

A glacier is the resampling unit. Each of `n_boot` replicates draws, with replacement, as many glaciers
as the row holds, keeps all the records of a drawn glacier (every year and elevation) and recomputes the
difference. A glacier drawn twice counts twice and one not drawn not at all; a region absent from a
replicate drops out of that replicate's regional mean. Records of one glacier are strongly correlated
(shared forcing, calibration and geometry), and a few glaciers hold most records, so resampling records
would report a precision the data do not have. `rng` seeds the draws; the default makes a call
reproducible.

Columns: `label_ref`, `label_member`, `stratum_by`, `stratum`, `source`, `season`, `metric` (`:abs_bias`
or `:rmse`), `n`, `n_glaciers`, `n_glacier_years`, `n_regions`, `ref` and `member` (the metric of each),
`difference`, `lo`, `median` and `hi` (the `(1 - level) / 2`, 0.5 and `(1 + level) / 2` quantiles of the
replicates' differences), `floor` and `distinguishable`.

`distinguishable` is true only when all of these hold:
  * the interval `[lo, hi]` excludes zero;
  * `abs(median)` exceeds `floor`;
  * the row holds at least `min_glaciers` glaciers. A bootstrap over fewer cannot resolve the interval's
    tails: the number of distinct resamples of `G` glaciers is `binomial(2G - 1, G)`, 126 for `G = 5`
    and 35 for `G = 4`, and a 95 % interval needs at least 40.

`floor` [m w.e.] is the size below which a difference cannot be attributed to the member's physics,
because the model represents the region's glaciers by one tile-pooled series. A number
(`representativeness_floor`) applies to every row. The default, `:glacier_scatter`, is the scatter
between glaciers of the observed balance at the same elevation within the row's source, region and
season ([`wgms_glacier_scatter`](@ref), on the paired records): two glaciers of one region at one
elevation differ by about that much in what is measured, so a member difference smaller than it is
smaller than the difference between the measured glacier and the one the model stands for. A row spanning
regions uses the mean of its regions' scatters, skipping regions that have no elevation bin with two
glaciers; it is `NaN` (and the row not distinguishable) where there is none. Only differences larger
than a glacier-to-glacier difference pass; smaller ones stay in the table with their interval.
"""
function wgms_member_difference(v_ref, v_member; n_boot=1000, rng=Random.Xoshiro(0), level=0.95,
                                min_n=30, min_glaciers=5, representativeness_floor=:glacier_scatter,
                                by=(:all, :region, :size_class, :provenance), size_edges_km2=(1, 10),
                                elevation_bin_m=250)
    n_boot >= 1 || error("n_boot must be positive, got $n_boot")
    0 < level < 1 || error("level must be in (0, 1), got $level")
    tail = (1 - level) / 2
    n_boot * tail >= 1 ||
        error("n_boot = $n_boot leaves fewer than one replicate in each $(100tail) % tail of a " *
              "$(100level) % interval; use at least $(ceil(Int, 1 / tail))")
    (representativeness_floor === :glacier_scatter ||
     (representativeness_floor isa Real && representativeness_floor >= 0)) ||
        error("representativeness_floor must be :glacier_scatter or a non-negative number, got " *
              "$(repr(representativeness_floor))")

    paired, label_ref, label_member = _wgms_pair_members(v_ref.matched, v_member.matched)
    scatter = Dict((r.source, r.region, r.season) => r.glacier_scatter
                   for r in eachrow(wgms_glacier_scatter(paired; elevation_bin_m)))
    row_floor(source, regions, season) = representativeness_floor isa Real ? Float64(representativeness_floor) :
        (x = filter(isfinite, [scatter[(source, r, season)] for r in regions]); isempty(x) ? NaN : mean(x))

    rows = NamedTuple[]
    for kind in by
        class, levels = _wgms_strata(paired, kind, size_edges_km2)
        for level_name in levels, source in sort(unique(paired.source)), season in sort(unique(paired.season))
            sub = paired[(class .== level_name) .& (paired.source .=== source) .& (paired.season .=== season), :]
            counts = DataFrames.combine(groupby(sub, :region), nrow => :n)
            regions = sort(counts.region[counts.n .>= min_n])
            isempty(regions) && continue
            sub = sub[in.(sub.region, Ref(Set(regions))), :]

            # Per-glacier sums, from which any resample's per-region sums follow by addition.
            glaciers = groupby(sub, [:region, :glacier_id])
            G = length(glaciers)
            stats = zeros(G, 5)
            region_of = Vector{Int}(undef, G)
            for (i, g) in zip(axes(stats, 1), glaciers)
                e0, e1 = g.error_ref, g.error_member
                stats[i, :] .= (length(e0), sum(e0), sum(e1), sum(abs2, e0), sum(abs2, e1))
                region_of[i] = searchsortedfirst(regions, first(g.region))
            end
            S = zeros(length(regions), 5)
            for i in axes(stats, 1)
                S[region_of[i], :] .+= view(stats, i, :)
            end
            est = _wgms_error_means(S)

            draws = Matrix{Float64}(undef, n_boot, 2)
            for b in axes(draws, 1)
                fill!(S, 0)
                for _ in 1:G
                    i = rand(rng, Base.OneTo(G))
                    S[region_of[i], :] .+= view(stats, i, :)
                end
                d = _wgms_error_means(S)
                draws[b, 1], draws[b, 2] = d[2] - d[1], d[4] - d[3]
            end

            fl = row_floor(source, regions, season)
            for (k, (metric, ref, member)) in enumerate(((:abs_bias, est[1], est[2]), (:rmse, est[3], est[4])))
                lo, med, hi = quantile(view(draws, :, k), (tail, 0.5, 1 - tail))
                push!(rows, (; label_ref, label_member, stratum_by=kind, stratum=level_name, source, season,
                             metric, n=nrow(sub), n_glaciers=G,
                             n_glacier_years=length(unique(zip(sub.glacier_id, sub.year))),
                             n_regions=length(regions), ref, member, difference=member - ref, lo,
                             median=med, hi, floor=fl,
                             distinguishable=(lo > 0 || hi < 0) && abs(med) > fl && G >= min_glaciers))
            end
        end
    end
    return DataFrame(rows)
end

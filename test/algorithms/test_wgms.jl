"""
Tests for the WGMS seasonal validation (utilities_wgms.jl).

A synthetic FoG data package is written to a temporary directory and read back. Model series are
built from closed forms: a balance rate linear in elevation and constant in time, so every interpolated
value and every window balance is known exactly.
"""

using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using GlobalGlacierAnalysis.NCDatasets
using DataFrames
using Dates
using Statistics
using LinearAlgebra
using Random

const _GLACIER_HEADER = "country,short_name,names,id,latitude,longitude,gtng_region,glims_id,rgi50_ids,rgi60_ids,rgi70_ids,wgi_id,parent_glacier_id,references,remarks"
const _MB_HEADER = "country,glacier_name,glacier_id,outline_id,year,time_system,begin_date,begin_date_unc,midseason_date,midseason_date_unc,end_date,end_date_unc,winter_balance,winter_balance_unc,summer_balance,summer_balance_unc,annual_balance,annual_balance_unc,ela_position,ela,ela_unc,aar,area,investigators,agencies,references,remarks"
const _POINT_HEADER = "country,glacier_name,glacier_id,year,id,original_id,begin_date,begin_date_unc,end_date,end_date_unc,latitude,longitude,elevation,balance,balance_unc,density,density_unc,balance_code,remarks"
const _BAND_HEADER = "country,glacier_name,glacier_id,year,lower_elevation,upper_elevation,area,winter_balance,winter_balance_unc,summer_balance,summer_balance_unc,annual_balance,annual_balance_unc,remarks"

# Glacier 1 at (46.5 N, 8.1 E), glacier 2 at (64.8 N, 18.9 W). Glacier-year (1, 2010): winter
# 2009-10-01 → 2010-05-01 (212 d), summer 2010-05-01 → 2010-09-30 (152 d), floating dates.
const _GLACIERS = ["CH,G ONE,,1,46.5,8.1,11_central_europe,,,,,,,,", "IS,G TWO,,2,64.8,-18.9,06_iceland,,,,,,,,"]
const _MB = ["CH,G ONE,1,,2010,floating,2009-10-01,0,2010-05-01,0,2010-09-30,0,9.9,,-9.9,,0,,,,,,,,,,",
             "IS,G TWO,2,,2010,,2009-10-15,0,2010-05-20,182.5,2010-09-20,0,,,,,,,,,,,,,,,"]

# Point rows by id; each documents the rule it exercises. Columns: glacier, year, id, begin, begin_unc,
# end, end_unc, lat, lon, elevation, balance, balance_unc, density, code, remarks.
point_row(g, y, id, b, bu, e, eu, lat, lon, z, bal, balu, rho, code, rem="") =
    "XX,G,$g,$y,$id,,$b,$bu,$e,$eu,$lat,$lon,$z,$bal,$balu,$rho,,$code,$rem"
const _POINTS = [
    point_row(1, 2010, 1, "2009-10-01", 0, "2010-05-01", 0, 46.51, 8.11, 2550, 1.2, 0.1, 400, "winter"),  # kept
    point_row(1, 2010, 2, "2009-07-02", 182, "2010-05-01", 0, 46.51, 8.11, 2650, 1.4, "", 450, "winter"),  # begin substituted
    point_row(1, 2010, 3, "2010-05-01", 0, "2010-09-30", 0, "", "", 2550, -2.0, "", "", "summer"),        # glacier location
    point_row(1, 2010, 4, "2009-10-01", 0, "2010-09-30", 0, 46.51, 8.11, 2550, -0.8, "", "", "annual"),
    point_row(1, 2010, 5, "2009-10-01", 0, "2010-09-30", 0, 46.51, 8.11, 2550, -0.8, "", "", "index"),
    point_row(1, 2010, 6, "2009-10-01", 0, "2010-05-01", 0, 46.51, 8.11, "", 1.0, "", "", "winter"),     # no elevation
    point_row(1, 2010, 7, "2009-10-01", 0, "2010-05-01", 0, 47.50, 8.11, 2550, 1.0, "", "", "winter"),   # 111 km away
    point_row(2, 2010, 8, "2009-10-15", 0, "2010-05-20", 182.5, 64.81, -18.91, 1050, 2.0, "", "", "winter"),  # mb midseason uncertain too
    point_row(1, 2010, 9, "2009-10-01", 0, "2010-09-30", 0, 46.51, 8.11, 2550, 1.0, "", "", "winter"),   # 364-day "winter"
    point_row(1, 2010, 10, "2010-05-01", 0, "2010-09-30", 0, 46.51, 8.11, 2550, -1.0, "", "", "summer", "Debris-covered ice"),
    point_row(1, 2010, 11, "2010-05-01", 0, "2010-09-30", 0, 46.51, 8.11, 2550, -20.0, "", "", "summer"), # implausible
    point_row(1, 2010, 12, "2010-05-01", 0, "2010-09-30", 0, 46.51, 8.11, 2550, -1.0, 0.8, "", "summer"), # uncertainty
    point_row(1, 2010, 13, "2009-10-01", 0, "2010-05-01", 0, 46.51, 8.11, 2550, 1.0, "", 900, "winter"),  # density
    point_row(1, 2010, 14, "2009-10-01", 0, "2010-05-01", 0, 46.52, 8.12, 2600, 1.10, "", "", "winter"),  # agreeing pair
    point_row(1, 2010, 15, "2009-10-01", 0, "2010-05-01", 0, 46.52, 8.12, 2600, 1.105, "", "", "winter"),
    point_row(1, 2010, 16, "2009-10-01", 0, "2010-05-01", 0, 46.53, 8.13, 2700, 1.0, "", "", "winter"),   # conflicting pair
    point_row(1, 2010, 17, "2009-10-01", 0, "2010-05-01", 0, 46.53, 8.13, 2700, 1.6, "", "", "winter"),
    point_row(1, 2010, 18, "2009-10-01", 0, "2010-05-01", 0, 46.54, 8.14, 2750, 1.3, "", 550, "winter",
              "estimated density of snow 550 kg/m3"),
]

const _BANDS = [
    "CH,G ONE,1,2010,2500,2600,1000000,1.0,,-2.5,,-1.5,,",   # both seasons kept
    "CH,G ONE,1,2010,2600,2700,1000000,1.2,,,,,,",            # summer: no balance
    "CH,G ONE,1,2010,2700,3000,1000000,1.4,,-1.0,,0.4,,",     # 300 m wide
    "IS,G TWO,2,2010,1000,1100,1000000,2.0,,-1.5,,0.5,,",      # midseason uncertain
]

function write_fog(dir; glaciers=_GLACIERS, mb=_MB, points=_POINTS, bands=_BANDS,
                   point_header=_POINT_HEADER)
    mkpath(dir)
    write(joinpath(dir, "glacier.csv"), join(vcat(_GLACIER_HEADER, glaciers), "\n") * "\n")
    write(joinpath(dir, "mass_balance.csv"), join(vcat(_MB_HEADER, mb), "\n") * "\n")
    write(joinpath(dir, "mass_balance_point.csv"), join(vcat(point_header, points), "\n") * "\n")
    write(joinpath(dir, "mass_balance_band.csv"), join(vcat(_BAND_HEADER, bands), "\n") * "\n")
    return dir
end

# A glacier-year row of mass_balance.csv with floating dates and the given area [m²].
function mb_row(g, y, area)
    f = fill("", 27)
    f[3], f[5], f[6], f[7], f[8] = string(g), string(y), "floating", "$(y - 1)-10-01", "0"
    f[9], f[10], f[11], f[12] = "$y-05-01", "0", "$y-09-30", "0"
    f[23] = string(area)
    return join(f, ",")
end
band_row(g, y, lo, hi) = "XX,G,$g,$y,$lo,$hi,,,,,,,,"

# A model whose balance rate is α + β z [m w.e. d⁻¹] at every band and time: levels on daily steps.
function linear_rate_model(centers, t0, t1; α=0.01, β=-1e-5, refreeze_rate=0.001, height_rate=0.02)
    time = collect(DateTime(t0):Day(1):DateTime(t1))
    d = [Dates.value(t - first(time)) / 86_400_000 for t in time]
    rate = α .+ β .* centers
    return GGA.GEMBBandSMB(centers, time, d * rate', d * fill(refreeze_rate, length(centers))',
                           d * fill(height_rate, length(centers))')
end

# Matched records of one configuration: glacier specs carry `id`, `region`, `area` [km²], `source`,
# `provenance`, `years`, `elevations`, `offset` (added to the observed balance) and the model errors
# `err_ref` and `err_member` (model − obs). `err` picks which of the two the model carries.
function synthetic_matched(label, specs, err; season=:winter)
    rows = [(; label, source=g.source, record_id="$(g.id)/$y/$z", glacier_id=g.id, year=y, region=g.region,
             season, latitude=46.0, longitude=8.0, z_lo=Float64(z), z_hi=Float64(z), glacier_area_km2=g.area,
             provenance=g.provenance, obs=1 + 0.002 * (z - 2000) + 0.1 * (y - 2000) + g.offset,
             model=1 + 0.002 * (z - 2000) + 0.1 * (y - 2000) + g.offset + getproperty(g, err))
            for g in specs for y in g.years for z in g.elevations]
    return DataFrame(rows)
end
spec(id, region, area, source, provenance; years=2000:2001, elevations=(2000, 2100), offset=0.0,
     err_ref=0.0, err_member=0.0) =
    (; id, region, area, source, provenance, years, elevations, offset, err_ref, err_member)

@testset "WGMS seasonal validation" begin
    dir = write_fog(joinpath(mktempdir(), "data"))
    wgms = GGA.read_wgms(dir)

    @testset "reader" begin
        @test nrow(wgms.glaciers) == 2 && nrow(wgms.points) == length(_POINTS) && nrow(wgms.bands) == 4
        @test wgms.glaciers.region == [11, 6]
        @test wgms.points.begin_date[2] == Date(2009, 7, 2) && wgms.points.begin_date_unc[2] == 182
        @test ismissing(wgms.points.latitude[3]) && ismissing(wgms.glacier_years.time_system[2])
        @test !("winter_balance" in names(wgms.glacier_years))   # glacier-wide balances are not read
        @test all(ismissing, wgms.glacier_years.area) && all(ismissing, wgms.glaciers.area_km2)

        bad(f; kw...) = GGA.read_wgms(write_fog(joinpath(mktempdir(), "data"); kw...))
        @test_throws "lacks column(s) balance_code" bad(nothing; point_header=replace(_POINT_HEADER, "balance_code" => "code"))
        @test_throws "begin_date = \"2009-13-01\" is not a yyyy-mm-dd date" bad(nothing;
            points=[replace(_POINTS[1], "2009-10-01" => "2009-13-01")])
        @test_throws "balance_code = \"spring\" is not one of" bad(nothing; points=[replace(_POINTS[1], "winter" => "spring")])
        @test_throws "latitude = 95.0 is above 90" bad(nothing; points=[replace(_POINTS[1], "46.51" => "95.0")])
        @test_throws "begin_date_unc = -1.0 is below 0" bad(nothing;
            points=[point_row(1, 2010, 1, "2009-10-01", -1, "2010-05-01", 0, 46.5, 8.1, 2550, 1.0, "", "", "winter")])
        @test_throws "share the key id" bad(nothing; points=[_POINTS[1], _POINTS[1]])
        @test_throws "glacier id(s) absent from glacier.csv: 9" bad(nothing;
            points=[point_row(9, 2010, 1, "2009-10-01", 0, "2010-05-01", 0, 46.5, 8.1, 2550, 1.0, "", "", "winter")])
        @test_throws "glacier-years absent from mass_balance.csv" bad(nothing; bands=["CH,G ONE,1,2011,2500,2600,,1.0,,,,,,"])
        @test_throws "area = -5.0 is below 0" bad(nothing; mb=[mb_row(1, 2010, -5)])
        @test_throws "gtng_region = \"21_mars\" is not one of" bad(nothing; glaciers=[replace(_GLACIERS[1], "11_central_europe" => "21_mars"), _GLACIERS[2]])
    end

    @testset "geotile ids" begin
        @test GGA.wgms_geotile_id(46.51, 8.11) == "lat[+46+48]lon[+008+010]"
        @test GGA.wgms_geotile_id(64.81, -18.91) == "lat[+64+66]lon[-020-018]"
        @test GGA.wgms_geotile_id(-43.1, 171.4) == "lat[-44-42]lon[+170+172]"
        @test GGA.wgms_geotile_id(46.0, 8.0) == "lat[+46+48]lon[+008+010]"     # an edge belongs north and east
        @test GGA.wgms_geotile_id(90.0, 180.0) == "lat[+88+90]lon[+178+180]"
        @test_throws "latitude 91 is outside" GGA.wgms_geotile_id(91, 0)
    end

    @testset "point quality control" begin
        P = GGA.wgms_seasonal_points(wgms)
        reasons = Dict(r.id => r.reason for r in eachrow(P.excluded))
        @test reasons == Dict(4 => "annual", 5 => "index", 6 => "no elevation", 7 => "mislocated",
                              8 => "date uncertain", 9 => "window length", 10 => "remarks",
                              11 => "implausible balance", 12 => "balance uncertainty",
                              13 => "implausible density", 15 => "duplicate",
                              16 => "conflicting duplicate", 17 => "conflicting duplicate")
        @test sort(P.records.id) == [1, 2, 3, 14, 18]
        r(id) = P.records[findfirst(==(id), P.records.id), :]
        # The placeholder begin (2009-07-02 ± 182 d) holds the glacier-year begin 2009-10-01.
        @test r(2).begin_date == Date(2009, 10, 1) && r(2).begin_date_unc == 0 && r(2).date_source == "mixed"
        @test r(1).date_source == "point" && r(1).geotile == "lat[+46+48]lon[+008+010]"
        @test r(3).location_source == "glacier" && (r(3).latitude, r(3).longitude) == (46.5, 8.1)
        @test r(3).season === :summer && ismissing(r(3).density) && r(3).time_system == "floating"
        # Glacier-year (1, 2010) reports winter densities 400, 450, 550 and 900: more than one value.
        @test r(1).density_class == "measured_varied" && r(14).density_class == "not_reported"
        @test r(18).density_class == "assumed"
        @test (r(1).z_lo, r(1).z_hi, r(1).obs) == (2550.0, 2550.0, 1.2)
    end

    @testset "band quality control" begin
        B = GGA.wgms_seasonal_bands(wgms)
        @test nrow(B.records) == 3
        kept = Set(zip(B.records.z_lo, B.records.season))
        @test kept == Set([(2500.0, :winter), (2500.0, :summer), (2600.0, :winter)])
        x = Dict((r.lower_elevation, r.season) => r.reason for r in eachrow(B.excluded))
        @test x[(2600.0, :summer)] == "no balance"
        @test x[(2700.0, :winter)] == "band width" && x[(2700.0, :summer)] == "band width"
        @test x[(1000.0, :winter)] == "date uncertain" && x[(1000.0, :summer)] == "date uncertain"
        w = B.records[(B.records.z_lo .== 2500) .& (B.records.season .=== :winter), :]
        s = B.records[(B.records.z_lo .== 2500) .& (B.records.season .=== :summer), :]
        @test (w.begin_date[1], w.end_date[1]) == (Date(2009, 10, 1), Date(2010, 5, 1))
        @test (s.begin_date[1], s.end_date[1]) == (Date(2010, 5, 1), Date(2010, 9, 30))
        @test w.record_id[1] == "1/2010/2500-2600" && w.date_source[1] == "glacier_year"
    end

    @testset "elevation weights" begin
        c = [1050.0, 1150.0, 1250.0]
        v = [3.0, 5.0, 4.0]
        @test GGA.wgms_band_weights(c, 1150, 1150) == [0.0, 1.0, 0.0]
        @test GGA.wgms_band_weights(c, 1075, 1075) ⋅ v ≈ 3.5
        @test GGA.wgms_band_weights(c, 1020, 1020) ⋅ v == 3.0          # held for ≤ 50 m
        @test GGA.wgms_band_weights(c, 1290, 1300) ⋅ v == 4.0
        @test isnothing(GGA.wgms_band_weights(c, 990, 990))
        @test isnothing(GGA.wgms_band_weights(c, 1200, 1310))
        # Mean of a linear piece is its midpoint value; across the 1150 knot the profile runs 4 → 5 → 4.5
        # and the integral is two trapezoids: ((4+5)/2·50 + (5+4.5)/2·50) / 100 = 4.625.
        @test GGA.wgms_band_weights(c, 1050, 1100) ⋅ v ≈ 3.5
        @test GGA.wgms_band_weights(c, 1100, 1200) ⋅ v ≈ 4.625
        # Partly beyond the outermost centre: constant 3 over [1000, 1050], then 3 → 4 over [1050, 1100].
        @test GGA.wgms_band_weights(c, 1000, 1100) ⋅ v ≈ (3 * 50 + 3.5 * 50) / 100
        @test sum(GGA.wgms_band_weights(c, 1010, 1290)) ≈ 1
        @test GGA.wgms_band_weights([2000.0], 1960, 2040) == [1.0]
        @test_throws "z_lo = 10 is above z_hi = 5" GGA.wgms_band_weights(c, 10, 5)
    end

    @testset "model balances over each record's own window" begin
        centers = collect(2050.0:100:2950.0)
        α, β = 0.01, -1e-5
        model = linear_rate_model(centers, Date(2000), Date(2020); α, β)
        rec = GGA.wgms_seasonal_points(wgms).records
        m = GGA.wgms_model_balances(rec, id -> id == "lat[+46+48]lon[+008+010]" ? model : nothing)
        @test nrow(m.matched) == 5 && isempty(m.excluded)
        for r in eachrow(m.matched)
            days = Dates.value(r.end_date - r.begin_date)
            @test r.model ≈ (α + β * r.z_lo) * days rtol = 1e-12
            @test r.model_refreeze ≈ 0.001 * days rtol = 1e-12
            @test r.model_dh ≈ 0.02 * days rtol = 1e-12
            # Date uncertainty is 0 for these records, so the envelope collapses.
            @test r.model_date_lo ≈ r.model_date_hi ≈ r.model
        end

        # Envelope: ends moved by ±10 days give windows of days ∓ 20.
        one = copy(rec[1:1, :]); one.begin_date_unc .= 10; one.end_date_unc .= 10
        e = GGA.wgms_model_balances(one, _ -> model).matched[1, :]
        days = Dates.value(e.end_date - e.begin_date)
        rate = α + β * e.z_lo
        @test rate < 0      # so the wider window gives the lower balance
        @test e.model_date_lo ≈ rate * (days + 20) && e.model_date_hi ≈ rate * (days - 20)

        # A level between output times is interpolated linearly: weekly levels of a constant rate give
        # the exact window balance at any day.
        weekly = GGA.GEMBBandSMB(centers, model.time[1:7:end], model.smb[1:7:end, :])
        mw = GGA.wgms_model_balances(rec, _ -> weekly).matched
        @test mw.model ≈ m.matched.model rtol = 1e-12
        @test all(ismissing, mw.model_refreeze)

        short = linear_rate_model(centers, Date(2009, 11, 1), Date(2020))
        x = GGA.wgms_model_balances(rec, _ -> short)
        @test x.excluded.reason == fill("outside model record", 4)   # the four winter records start 2009-10-01
        high = linear_rate_model([3550.0, 3650.0], Date(2000), Date(2020))
        @test all(==("outside band range"), GGA.wgms_model_balances(rec, _ -> high).excluded.reason)
    end

    @testset "GEMBBandSMB checks its inputs" begin
        t = [DateTime(2000), DateTime(2000, 1, 2)]
        @test_throws "band centers are not strictly ascending" GGA.GEMBBandSMB([2.0, 1.0], t, zeros(2, 2))
        @test_throws DimensionMismatch GGA.GEMBBandSMB([1.0, 2.0], t, zeros(3, 2))
        @test_throws "smb holds non-finite values" GGA.GEMBBandSMB([1.0, 2.0], t, [0.0 NaN; 0.0 0.0])
        @test_throws "not strictly ascending" GGA.GEMBBandSMB([1.0, 2.0], reverse(t), zeros(2, 2))
    end

    @testset "tile file at fitted corrections, bilinear across the forcing nodes" begin
        tile_dir = mkpath(joinpath(mktempdir(), "tile_runs_v2"))
        path = joinpath(tile_dir, "N46_E008.nc")
        times = DateTime(2009, 9, 6, 23) .+ Week.(0:60)
        centers = [2450.0, 2550.0, 2650.0]
        dts, pss = [2.0, -1.0], [2.0, 1.0]       # descending on file; the reader sorts
        # Per interval (time, band, ΔT, pscale), linear in both corrections:
        #   P = 10·p, runoff = 3·(ΔT + 1) + 0.01·(z − 2450), EC = −0.5, refreeze = ΔT + 2.
        P = [10.0 * p for t in times, z in centers, d in dts, p in pss]
        Ro = [3.0 * (d + 1) + 0.01 * (z - 2450) for t in times, z in centers, d in dts, p in pss]
        EC = fill(-0.5, size(P))
        Rf = [d + 2.0 for t in times, z in centers, d in dts, p in pss]
        NCDataset(path, "c") do ds
            ds.attrib["geotile_id"] = "lat[+46+48]lon[+008+010]"
            defDim(ds, "time", length(times)); defDim(ds, "band", 3)
            defDim(ds, "delta_temperature", 2); defDim(ds, "precipitation_scaling", 2)
            defVar(ds, "time", times, ("time",); attrib=Dict("units" => "hours since 1950-01-01 00:00:00"))
            defVar(ds, "band_center", centers, ("band",))
            defVar(ds, "delta_temperature", dts, ("delta_temperature",))
            defVar(ds, "precipitation_scaling", pss, ("precipitation_scaling",))
            dimsv = ("time", "band", "delta_temperature", "precipitation_scaling")
            for (name, a) in ("precipitation" => P, "runoff" => Ro, "evaporation_condensation" => EC, "refreeze" => Rf)
                defVar(ds, name, Float32.(a), dimsv)
            end
        end
        p, d = 1.4, 0.5
        m = GGA.gemb_tile_band_smb(path; pscale=p, ΔT=d)
        @test m.band_centers == centers && m.time == times
        k = collect(1:length(times))
        for (j, z) in enumerate(centers)
            @test m.smb[:, j] ≈ k .* (10p - 3(d + 1) - 0.01(z - 2450) - 0.5) ./ 1000 rtol = 1e-6
            @test m.refreeze[:, j] ≈ k .* (d + 2) ./ 1000 rtol = 1e-6
        end
        @test isnothing(m.height)
        @test_throws "outside the forcing grid" GGA.gemb_tile_band_smb(path; pscale=2.5, ΔT=0.0)

        fits = DataFrame(id=["lat[+46+48]lon[+008+010]", "lat[+64+66]lon[-020-018]"], pscale=[p, 1.0], ΔT=[d, 0.0])
        model = GGA.wgms_tile_model(tile_dir, fits)
        @test model("lat[+46+48]lon[+008+010]").smb ≈ m.smb
        @test isnothing(model("lat[+64+66]lon[-020-018]"))     # fitted, but no tile file
        uncal = GGA.wgms_tile_model(tile_dir, fits; pscale=1.0, ΔT=0.0)
        @test uncal("lat[+46+48]lon[+008+010]").smb[end, 1] ≈ length(times) * (10 - 3 - 0.5) / 1000 rtol = 1e-6
    end

    @testset "statistics" begin
        # Two glaciers, two years, three elevations; obs = 1 + 0.002 (z − 2000).
        rows = [(source=:point, record_id="$g/$y/$z", glacier_id=g, year=y, region=11, season=s,
                 latitude=46.0 + g, longitude=8.0, z_lo=Float64(z), z_hi=Float64(z),
                 obs=(s === :winter ? 1 : -2) * (1 + 0.002 * (z - 2000)) + 0.1 * (y - 2010))
                for g in 1:2, y in 2010:2015, z in (2000, 2200, 2400), s in (:winter, :summer)]
        df = DataFrame(vec(rows))
        df.model = 2 .* df.obs .+ 0.5
        st = GGA.wgms_seasonal_statistics(df)
        for r in eachrow(st)
            o = df.obs[df.season .=== r.season]
            @test r.n == 36 && r.n_glaciers == 2 && r.n_glacier_years == 12
            @test r.bias ≈ mean(o) + 0.5
            @test r.rmse ≈ sqrt(mean((o .+ 0.5) .^ 2))
            @test r.r ≈ 1
        end
        eb = GGA.wgms_seasonal_statistics(df; elevation_bin_m=250)
        @test sort(unique(eb.z_bin)) == [2000.0, 2250.0] && sum(eb.n) == nrow(df)

        g = GGA.wgms_vertical_gradients(df)
        @test nrow(g) == 2 * 6 * 2 && all(g.span_m .== 400)
        @test all(g.z_min .== 2000) && all(g.z_max .== 2400)
        @test all(g.obs_gradient_se .< 1e-12)            # exactly linear in z
        @test all(g.obs_gradient[g.season .=== :winter] .≈ 0.2)
        @test all(g.model_gradient[g.season .=== :summer] .≈ -0.8)
        @test nrow(GGA.wgms_vertical_gradients(df; min_span_m=500)) == 0

        a = GGA.wgms_anomaly_correlation(df)
        @test all(a.n_sites .== 6) && all(a.r_anomaly .≈ 1)
        @test all(a.model_anomaly_sd .≈ 2 .* a.obs_anomaly_sd)
        @test isempty(GGA.wgms_anomaly_correlation(df; min_years=7))

        ratio = GGA.wgms_seasonal_ratio(df)
        w, s = df.season .=== :winter, df.season .=== :summer
        @test ratio.n_pairs == [36]
        @test ratio.obs_ratio[1] ≈ mean(df.obs[w]) / mean(-df.obs[s])
        @test ratio.model_ratio[1] ≈ mean(df.model[w]) / mean(-df.model[s])
    end

    @testset "snow densities from a column" begin
        dz, rho = [0.5, 0.5, 1.0], [300.0, 400.0, 500.0]
        @test GGA.top_mass_bulk_density(dz, rho, 150) ≈ 300
        @test GGA.top_mass_bulk_density(dz, rho, 350) ≈ 350 / (0.5 + 0.5)
        @test GGA.top_mass_bulk_density(dz, rho, 600) ≈ 600 / (1.0 + 250 / 500)
        @test_throws "holds 850.0 kg m-2, less than" GGA.top_mass_bulk_density(dz, rho, 900)
        age = [20.0, 100.0, 400.0]
        y = GGA.young_snow_bulk_density(dz, rho, age, 200)
        @test (y.mass, y.depth, y.density) == (350.0, 1.0, 350.0)
        @test isnan(GGA.young_snow_bulk_density(dz, rho, age, 10).density)
    end

    @testset "end to end: the true model scores zero error" begin
        centers = collect(2450.0:100:2750.0)
        # Construct the model from the records themselves: a rate per day that reproduces each
        # observation needs one rate per record, so use the band source on a single glacier-year and
        # rates fitted to its two winter bands: obs 1.0 at 2550 and 1.2 at 2650 over 212 days.
        rate = z -> (1.0 + 0.002 * (z - 2550)) / 212
        time = collect(DateTime(2009, 10, 1):Day(1):DateTime(2010, 5, 1))
        d = [Dates.value(t - first(time)) / 86_400_000 for t in time]
        model = GGA.GEMBBandSMB(centers, time, d * rate.(centers)')
        v = GGA.wgms_validation(id -> id == "lat[+46+48]lon[+008+010]" ? model : nothing; wgms, label="truth")
        bw = v.matched[(v.matched.source .=== :band) .& (v.matched.season .=== :winter), :]
        @test nrow(bw) == 2 && all(v.matched.label .== "truth")
        # Covariates: the fixture reports no glacier area; both bands lie inside the winter stakes'
        # elevation range [2550, 2750].
        @test all(ismissing, v.matched.glacier_area_km2)
        @test bw.provenance == [:band_in_stake_range, :band_in_stake_range]
        @test all(==(:stake), v.matched.provenance[v.matched.source .=== :point])
        @test bw.model ≈ bw.obs rtol = 1e-12
        @test v.summary[(v.summary.source .=== :band) .& (v.summary.season .=== :winter), :rmse][1] < 1e-12
        # Summer windows reach past the model series and are reported, not dropped silently.
        @test sum(v.excluded_model.n[v.excluded_model.reason .== "outside model record"]) == 2
        @test sum(v.excluded_qc.n[v.excluded_qc.source .=== :point]) == 13

        card = GGA.wgms_scorecard(v; min_n=1)
        bw_card = card[(card.source .=== :band) .& (card.season .=== :winter), :]
        @test nrow(bw_card) == 1 && bw_card.label[1] == "truth" && bw_card.n[1] == 2
        @test bw_card.mean_abs_bias[1] < 1e-12 && bw_card.mean_rmse[1] < 1e-12
        @test isnan(bw_card.median_r_anomaly[1])            # one year: no anomaly
        @test nrow(GGA.wgms_scorecard(v; min_n=100)) == 0
    end

    @testset "scorecard weights regions equally" begin
        summary = DataFrame(label="x", source=:point, region=[1, 2, 3], season=:winter, n=[100, 10, 40],
                            n_glaciers=[5, 1, 3], n_glacier_years=[50, 4, 20],
                            bias=[0.1, -5.0, -0.3], rmse=[0.2, 5.0, 0.6], r=[0.5, 0.9, 0.7])
        anomaly = DataFrame(label="x", source=:point, region=[1, 3], season=:winter, r_anomaly=[0.2, 0.6])
        grads = DataFrame(label="x", source=:point, region=[1, 3], season=:winter, median_difference=[0.1, 0.3])
        ratio = DataFrame(label="x", source=:point, region=[1, 3], obs_ratio=[1.0, 2.0], model_ratio=[exp(0.5), 2.0])
        v = (; summary, anomaly, gradient_summary=grads, ratio)
        c = GGA.wgms_scorecard(v)                           # region 2 has 10 < 30 records
        @test c.n_regions == [2] && c.n == [140]
        @test c.n_glaciers == [8] && c.n_glacier_years == [70] && c.restriction == ["all"]
        @test c.mean_abs_bias[1] ≈ 0.2 && c.mean_rmse[1] ≈ 0.4 && c.median_r[1] ≈ 0.6
        @test c.median_r_anomaly[1] ≈ 0.4 && c.median_gradient_difference[1] ≈ 0.2
        @test c.median_abs_log_ratio_error[1] ≈ 0.25

        grads.median_abs_difference = [0.1, 0.5]
        grads.obs_median = [0.4, 0.2]; grads.model_median = [0.1, 0.2]
        c = GGA.wgms_scorecard((; summary, anomaly, gradient_summary=grads, ratio))
        @test c.median_abs_gradient_difference[1] ≈ 0.3 && c.median_gradient_ratio[1] ≈ (0.25 + 1) / 2

        # Three configurations: y is best on bias and gradient, z on anomaly r.
        cards = DataFrame(label=["x", "y", "z"], source=:point, season=:winter,
                          mean_abs_bias=[0.5, 0.1, 0.5], median_abs_gradient_difference=[0.3, 0.1, NaN],
                          median_r_anomaly=[0.4, 0.5, 0.9])
        rk = GGA.wgms_rank(cards)
        @test rk.label == ["y", "z", "x"]
        @test rk.rank_mean_abs_bias == [1, 2, 2]                 # ties share the lower rank
        @test rk.rank_median_abs_gradient_difference == [1, 3, 2] # NaN ranks last
        @test rk.rank_median_r_anomaly == [2, 1, 3]
        @test rk.mean_rank ≈ [4 / 3, 2, 7 / 3]
        @test_throws "must be :lower or :higher" GGA.wgms_rank(cards; metrics=(mean_abs_bias=:best,))
    end

    @testset "record covariates: area, provenance and stake span" begin
        pts = [point_row(1, 2010, 1, "2009-10-01", 0, "2010-05-01", 0, 46.51, 8.11, 2400, 1.0, "", "", "winter"),
               point_row(1, 2010, 2, "2009-10-01", 0, "2010-05-01", 0, 46.51, 8.11, 2600, 1.2, "", "", "winter"),
               point_row(1, 2010, 3, "2009-10-01", 0, "2010-05-01", 0, 46.51, 8.11, 2300, "", "", "", "winter"),    # no balance
               point_row(1, 2010, 4, "2009-10-01", 0, "2010-05-01", 0, 46.51, 8.11, "", 1.0, "", "", "winter"),     # no elevation
               point_row(1, 2010, 5, "2009-10-01", 0, "2010-09-30", 0, 46.51, 8.11, 2900, 0.2, "", "", "annual"),   # not a seasonal stake
               point_row(1, 2010, 6, "2009-10-01", 0, "2010-05-01", 0, 46.51, 8.11, 2700, 12.0, "", "", "winter"),  # fails QC: implausible
               point_row(1, 2010, 7, "2010-05-01", 0, "2010-09-30", 0, 46.51, 8.11, 2500, -2.0, "", "", "summer"),
               point_row(2, 2010, 8, "2009-10-15", 0, "2010-05-20", 0, 64.81, -18.91, 1050, 12.0, "", "", "winter")] # fails QC
        mb = [mb_row(1, 2010, 1_000_000), mb_row(1, 2011, 3_000_000), mb_row(1, 2012, 2_000_000),
              mb_row(2, 2010, "")]
        w = GGA.read_wgms(write_fog(joinpath(mktempdir(), "data"); mb, points=pts,
                                    bands=[band_row(1, 2010, 2000, 2100), band_row(1, 2010, 2900, 3000)]))
        @test isequal(w.glaciers.area_km2, [2.0, missing])            # median of 1, 3 and 2 km²
        @test sort(GGA.wgms_seasonal_points(w).records.id) == [1, 2, 7]

        rec(source, g, y, season, lo, hi; id="$source/$g/$y/$season/$lo") =
            (; source, record_id=id, glacier_id=g, year=y, season, z_lo=Float64(lo), z_hi=Float64(hi))
        # Glacier-year-season (1, 2010, winter): stakes at 2400, 2600 (QC), 2700 (fails QC) → database
        # range [2400, 2700], QC range [2400, 2600]; the 50 m tolerance extends each side.
        records = DataFrame([
            rec(:band, 1, 2010, :winter, 2300, 2400),    # centre 2350: the lower edge, inclusive
            rec(:band, 1, 2010, :winter, 2200, 2300),    # 2250: beyond
            rec(:band, 1, 2010, :winter, 2700, 2800),    # 2750: the upper edge, inclusive (QC: beyond)
            rec(:band, 1, 2010, :winter, 2750, 2850),    # 2800: beyond
            rec(:band, 1, 2010, :summer, 2500, 2600),    # 2550 against the single summer stake at 2500
            rec(:band, 1, 2010, :summer, 2600, 2700),    # 2650: beyond
            rec(:band, 1, 2011, :winter, 2500, 2600),    # no stake in 2011
            rec(:band, 2, 2010, :winter, 1000, 1100),    # the only stake (1050) fails QC
            rec(:point, 1, 2010, :winter, 2400, 2400; id="1"),
            rec(:point, 1, 2010, :summer, 2500, 2500; id="7"),
            rec(:point, 2, 2010, :winter, 1050, 1050; id="8")])
        c = GGA.wgms_record_covariates(w, records)
        @test c.provenance == [:band_in_stake_range, :band_beyond_stakes, :band_in_stake_range,
                               :band_beyond_stakes, :band_in_stake_range, :band_beyond_stakes,
                               :band_no_stakes, :band_in_stake_range, :stake, :stake, :stake]
        @test all(ismissing, c.glacier_area_km2[c.glacier_id .== 2]) && all(==(2.0), c.glacier_area_km2[c.glacier_id .== 1])
        # Span over the band range [2000, 3000] of glacier 1: winter stakes 2400-2700 (database), the
        # single summer stake 0; glacier 2 has no bands.
        @test c.stake_span_fraction[9:10] ≈ [0.3, 0.0]
        @test ismissing(c.stake_span_fraction[11]) && all(ismissing, c.stake_span_fraction[1:8])
        @test names(c)[end-2:end] == ["glacier_area_km2", "provenance", "stake_span_fraction"]

        # Counting only the stakes that pass quality control narrows the range.
        q = GGA.wgms_record_covariates(w, records[1:8, :]; stakes=:qc)
        @test q.provenance[[3, 8]] == [:band_beyond_stakes, :band_no_stakes]
        qp = GGA.wgms_record_covariates(w, records[[9], :]; stakes=:qc)
        @test qp.stake_span_fraction ≈ [0.2]
        # Tolerance 0 narrows "in range" to the stakes' own extent.
        t0 = GGA.wgms_record_covariates(w, records[1:4, :]; tolerance_m=0)
        @test t0.provenance == [:band_beyond_stakes, :band_beyond_stakes, :band_beyond_stakes, :band_beyond_stakes]
        @test GGA.wgms_record_covariates(w, records[1:0, :]).provenance == Symbol[]

        @test_throws "stakes must be :database or :qc, got :all" GGA.wgms_record_covariates(w, records; stakes=:all)
        @test_throws "tolerance_m must be non-negative" GGA.wgms_record_covariates(w, records; tolerance_m=-1)
        @test_throws "records lack column(s) year" GGA.wgms_record_covariates(w, select(records, Not(:year)))
        @test_throws "source must be :point or :band, got :pit" GGA.wgms_record_covariates(
            w, transform(records[1:1, :], :source => (x -> [:pit]) => :source))
        @test_throws "glacier 9 is not in wgms.glaciers" GGA.wgms_record_covariates(
            w, transform(records[1:1, :], :glacier_id => (x -> [9]) => :glacier_id))
        @test_throws "point record 8 has no stake under stakes=:qc" GGA.wgms_record_covariates(
            w, records[[11], :]; stakes=:qc)
        @test_throws "no area_km2 column" GGA.wgms_record_covariates(
            (; glaciers=select(w.glaciers, Not(:area_km2)), points=w.points, bands=w.bands), records)
    end

    @testset "size classes" begin
        @test GGA.wgms_size_class(0.5) == "<1" && GGA.wgms_size_class(1) == "1-10"
        @test GGA.wgms_size_class(9.99) == "1-10" && GGA.wgms_size_class(10) == ">=10"
        @test GGA.wgms_size_class(missing) == "missing"
        @test GGA.wgms_size_class(3; edges_km2=(2,)) == ">=2" && GGA.wgms_size_class(1; edges_km2=(2,)) == "<2"
        @test GGA.wgms_size_class(30; edges_km2=(1, 5, 20)) == ">=20" && GGA.wgms_size_class(7; edges_km2=(1, 5, 20)) == "5-20"
        @test_throws "strictly ascending" GGA.wgms_size_class(1; edges_km2=(5, 1))
    end

    @testset "strata report counts and metrics on every row" begin
        # Constant model − obs error per glacier, so each region's bias is a known mean.
        specs = [spec(1, 1, 0.5, :point, :stake; err_ref=0.5), spec(2, 1, 5.0, :point, :stake; err_ref=-0.2),
                 spec(3, 1, 20.0, :band, :band_in_stake_range; err_ref=0.1),
                 spec(4, 2, missing, :band, :band_no_stakes; err_ref=0.3),
                 spec(5, 2, 2.0, :band, :band_beyond_stakes; err_ref=-0.4)]
        v = GGA._wgms_rescore(synthetic_matched("m", specs, :err_ref))
        s = GGA.wgms_scorecard_strata(v; min_n=1)
        row(by, level, source) = only(eachrow(s[(s.stratum_by .=== by) .& (s.stratum .== level) .& (s.source .=== source), :]))
        counts(r) = (r.n, r.n_glaciers, r.n_glacier_years)
        @test counts(row(:all, "all", :point)) == (8, 2, 4) && counts(row(:all, "all", :band)) == (12, 3, 6)
        @test counts(row(:size_class, "<1", :point)) == (4, 1, 2) && counts(row(:size_class, "1-10", :point)) == (4, 1, 2)
        @test counts(row(:size_class, "1-10", :band)) == (4, 1, 2)
        @test counts(row(:size_class, ">=10", :band)) == (4, 1, 2) && counts(row(:size_class, "missing", :band)) == (4, 1, 2)
        @test !any((s.stratum_by .=== :size_class) .& (s.stratum .== ">=10") .& (s.source .=== :point))   # no record: no row
        @test counts(row(:provenance, "stake", :point)) == (8, 2, 4) && counts(row(:provenance, "band_no_stakes", :band)) == (4, 1, 2)
        @test counts(row(:provenance, "band_in_stake_range", :band)) == (4, 1, 2)
        @test counts(row(:region, "2", :band)) == (8, 2, 4)
        # Every row carries the three counts, for every stratification.
        @test all(s.n .>= s.n_glacier_years .>= s.n_glaciers .>= 1) && all(s.n_scored .== s.n)
        @test Set(s.stratum_by) == Set([:all, :size_class, :provenance, :region])

        # Region-equal means of the pooled band rows: region 1 holds glacier 3 (bias 0.1), region 2 holds
        # glaciers 4 and 5 (mean error -0.05, RMS sqrt((0.09 + 0.16) / 2)).
        r = row(:all, "all", :band)
        @test r.n_regions == 2 && r.mean_abs_bias ≈ (0.1 + 0.05) / 2 && r.mean_rmse ≈ (0.1 + sqrt(0.125)) / 2
        @test row(:size_class, "<1", :point).mean_abs_bias ≈ 0.5

        # Regions with fewer than min_n records in a stratum are counted but not scored.
        s5 = GGA.wgms_scorecard_strata(v; min_n=5, by=(:all, :size_class))
        small = only(eachrow(s5[(s5.stratum .== "<1") .& (s5.source .=== :point), :]))
        @test small.n == 4 && small.n_scored == 0 && small.n_regions == 0 && isnan(small.mean_abs_bias)
        big = only(eachrow(s5[(s5.stratum .== "1-10") .& (s5.source .=== :band), :]))
        @test big.n == 4 && big.n_scored == 0
        @test_throws "unknown stratification :shape" GGA.wgms_scorecard_strata(v; by=(:shape,))
        @test_throws "lacks column(s) glacier_area_km2" GGA.wgms_scorecard_strata(
            (; matched=select(v.matched, Not(:glacier_area_km2))))

        # Scorecard restricted to anchored records: only the stake and in-range classes remain.
        anchored = GGA.wgms_scorecard(v; min_n=1, provenance=GGA.WGMS_ANCHORED_PROVENANCE)
        @test anchored.restriction == fill("stake+band_in_stake_range", 2)
        @test sort(anchored.n) == [4, 8]
        ab = anchored[anchored.source .=== :band, :]
        @test ab.n == [4] && ab.mean_abs_bias ≈ [0.1]
        @test_throws "unknown provenance class(es) :nonsense" GGA.wgms_scorecard(v; provenance=(:nonsense,))
        @test_throws "v.matched has no provenance column" GGA.wgms_scorecard((; summary=v.summary); provenance=(:stake,))
    end

    @testset "scatter between glaciers at matched elevations" begin
        # One 250 m bin; four glaciers with 2 records each, offsets 0..3 plus ± 0.5 around the offset.
        rows = [(; source=:band, region=1, season=:winter, glacier_id=g, z_lo=2000.0, z_hi=2000.0, obs=o + d)
                for (g, o) in enumerate((0.0, 1.0, 2.0, 3.0)) for d in (-0.5, 0.5)]
        # A second bin (2250-2500) with glaciers 1-3 at offsets 0, 2, 6 and one lone glacier in a third.
        append!(rows, [(; source=:band, region=1, season=:winter, glacier_id=g, z_lo=2300.0, z_hi=2400.0, obs=o)
                       for (g, o) in ((1, 0.0), (2, 2.0), (3, 6.0))])
        push!(rows, (; source=:band, region=1, season=:winter, glacier_id=9, z_lo=3000.0, z_hi=3000.0, obs=50.0))
        sc = GGA.wgms_glacier_scatter(DataFrame(rows))
        # Bin 1: sum of squares of [0, 1, 2, 3] about 1.5 is 5 with 3 df; bin 2: [0, 2, 6] about 8/3 is 18.67
        # with 2 df; the lone glacier adds none.
        @test nrow(sc) == 1 && sc.n_glaciers == [5] && sc.n_dof == [5]
        @test sc.glacier_scatter[1] ≈ sqrt((5 + 56 / 3) / 5)
        lone = GGA.wgms_glacier_scatter(DataFrame(rows[end:end]))
        @test lone.n_dof == [0] && isnan(lone.glacier_scatter[1])
        @test_throws "elevation_bin_m must be positive" GGA.wgms_glacier_scatter(DataFrame(rows); elevation_bin_m=0)
    end

    @testset "glacier-block bootstrap of a member difference" begin
        # Region 1: one glacier of 300 records favours the member (error 0 against the reference's
        # 0.5); five glaciers of two records favour the reference (0 against 0.5).
        big = spec(1, 1, 10.0, :band, :band_no_stakes; years=2000:2029, elevations=ntuple(i -> 2000 + 10i, 10),
                   err_ref=0.5, err_member=0.0)
        smalls = [spec(1 + g, 1, 2.0, :band, :band_no_stakes; years=2000:2001, elevations=(2000,),
                       err_ref=0.0, err_member=0.5) for g in 1:5]
        vs = Dict(l => GGA._wgms_rescore(synthetic_matched(l, [big; smalls], e))
                  for (l, e) in (("ref", :err_ref), ("member", :err_member)))
        d = GGA.wgms_member_difference(vs["ref"], vs["member"]; by=(:all,), representativeness_floor=0)
        ab = only(eachrow(d[d.metric .=== :abs_bias, :]))
        @test (ab.n, ab.n_glaciers, ab.n_regions) == (310, 6, 1) && ab.label_ref == "ref" && ab.label_member == "member"
        @test ab.ref ≈ 150 / 310 && ab.member ≈ 5 / 310 && ab.difference ≈ (5 - 150) / 310
        # Without the big glacier a third of the replicates have the member worse; the interval spans 0.
        @test ab.lo < 0 < ab.hi && !ab.distinguishable
        @test all(.!d.distinguishable) && all(d.lo .< 0 .< d.hi)

        # Resampling the 310 records would call the same difference significant.
        e0 = vs["ref"].matched.model .- vs["ref"].matched.obs
        e1 = vs["member"].matched.model .- vs["member"].matched.obs
        rng = Random.Xoshiro(1)
        record_level = [(i = rand(rng, 1:310, 310); abs(mean(e1[i])) - abs(mean(e0[i]))) for _ in 1:1000]
        @test quantile(record_level, 0.975) < 0

        # The default seed makes a call reproducible.
        again = GGA.wgms_member_difference(vs["ref"], vs["member"]; by=(:all,), representativeness_floor=0)
        @test isequal(d, again)
    end

    @testset "distinguishable needs the interval, the floor and enough glaciers" begin
        # Ten glaciers, five in each of two regions, every one favouring the member. Offsets to the observed
        # balance give the regions a scatter between glaciers.
        off1, off2 = 0.3 .* (0:4), 0.1 .* (0:4)
        specs = vcat([spec(g, 1, 3.0, :band, :band_in_stake_range; years=2000:2002, elevations=(2000, 2100),
                           offset=off1[g], err_ref=0.4 + 0.02g, err_member=0.1 + 0.01g) for g in 1:5],
                     [spec(g, 2, 3.0, :band, :band_in_stake_range; years=2000:2002, elevations=(2000, 2100),
                           offset=off2[g-5], err_ref=0.4 + 0.02g, err_member=0.1 + 0.01g) for g in 6:10])
        ref = GGA._wgms_rescore(synthetic_matched("ref", specs, :err_ref))
        mem = GGA._wgms_rescore(synthetic_matched("mem", specs, :err_member))
        diff(; kw...) = GGA.wgms_member_difference(ref, mem; by=(:all, :region), kw...)
        pick(d, by, level, metric=:abs_bias) = only(eachrow(d[(d.stratum_by .=== by) .& (d.stratum .== level) .& (d.metric .=== metric), :]))

        # |bias| is the mean error: region 1 holds g = 1:5, region 2 g = 6:10.
        mean_ref(g) = mean(0.4 .+ 0.02 .* g)
        mean_mem(g) = mean(0.1 .+ 0.01 .* g)
        all0 = pick(diff(representativeness_floor=0), :all, "all")
        @test all0.difference ≈ (mean_mem(1:5) + mean_mem(6:10)) / 2 - (mean_ref(1:5) + mean_ref(6:10)) / 2
        @test all0.n_regions == 2 && all0.n_glaciers == 10 && all0.n_glacier_years == 30
        @test all0.lo < all0.median < all0.hi < 0 && all0.distinguishable
        @test all0.difference ≈ -0.355

        # A number sets the floor of every row.
        @test pick(diff(representativeness_floor=0.3), :all, "all").distinguishable
        @test !pick(diff(representativeness_floor=0.4), :all, "all").distinguishable
        @test pick(diff(representativeness_floor=0.4), :all, "all").floor == 0.4
        # Fewer glaciers than min_glaciers: the interval is reported, the flag is not set.
        few = diff(representativeness_floor=0, min_glaciers=11)
        @test !any(few.distinguishable) && all(few.hi .< 0)

        # The default floor is the scatter between glaciers within each region and season.
        sd(x) = std(x)
        d = diff()
        r1, r2 = pick(d, :region, "1"), pick(d, :region, "2")
        @test r1.floor ≈ sd(off1) && r2.floor ≈ sd(off2)
        @test pick(d, :all, "all").floor ≈ (sd(off1) + sd(off2)) / 2
        @test !r1.distinguishable && r2.distinguishable
        @test abs(r1.median) < sd(off1) < 0.5     # region 1 scatters more than its members differ
        @test pick(d, :all, "all").distinguishable          # 0.355 exceeds the mean scatter 0.316
        # RMSE is reported alongside |bias|, with the same pairing.
        @test pick(d, :all, "all", :rmse).n == pick(d, :all, "all", :abs_bias).n

        # Reproducible by seed; another seed moves the interval, not the point estimate.
        a = diff(rng=Random.Xoshiro(7))
        @test isequal(a, diff(rng=Random.Xoshiro(7)))
        b = diff(rng=Random.Xoshiro(8))
        @test a.lo != b.lo && a.hi != b.hi && a.difference == b.difference

        # Members are paired on the records both matched (region 2 keeps 24 records, so min_n = 24).
        trimmed = (; matched=mem.matched[mem.matched.glacier_id .!= 10, :])
        t = GGA.wgms_member_difference(ref, trimmed; by=(:all,), representativeness_floor=0, min_n=24)
        @test all(t.n .== 54) && all(t.n_glaciers .== 9)

        @test_throws "leaves fewer than one replicate" GGA.wgms_member_difference(ref, mem; n_boot=39)
        @test_throws "level must be in (0, 1)" GGA.wgms_member_difference(ref, mem; level=1.0)
        @test_throws "representativeness_floor must be" GGA.wgms_member_difference(ref, mem; representativeness_floor=-1)
        @test_throws "representativeness_floor must be" GGA.wgms_member_difference(ref, mem; representativeness_floor=:median)
        @test_throws "unknown stratification :shape" GGA.wgms_member_difference(ref, mem; by=(:shape,))
        shifted = (; matched=transform(mem.matched, :obs => (x -> x .+ 0.1) => :obs))
        @test_throws "disagree on the observed balance" GGA.wgms_member_difference(ref, shifted)
        disjoint = (; matched=transform(mem.matched, :record_id => (x -> "z" .* x) => :record_id))
        @test_throws "share no matched record" GGA.wgms_member_difference(ref, disjoint)
        dup = (; matched=vcat(mem.matched, mem.matched[1:1, :]))
        @test_throws "v_member.matched holds duplicate" GGA.wgms_member_difference(ref, dup)
        @test_throws "v_ref.matched lacks column(s) provenance" GGA.wgms_member_difference(
            (; matched=select(ref.matched, Not(:provenance))), mem)
    end

    @testset "rankings by provenance restriction" begin
        # Two members on three glaciers: X is closer on the anchored classes, Y on bands without stakes.
        specs(ex, ey) = [spec(1, 1, 3.0, :band, :band_in_stake_range; err_ref=ex[1], err_member=ey[1]),
                         spec(2, 1, 3.0, :band, :band_no_stakes; err_ref=ex[2], err_member=ey[2]),
                         spec(3, 1, 3.0, :point, :stake; err_ref=ex[3], err_member=ey[3])]
        sp = specs((0.1, 0.5, 0.1), (0.3, 0.2, 0.3))
        X = GGA._wgms_rescore(synthetic_matched("X", sp, :err_ref))
        Y = GGA._wgms_rescore(synthetic_matched("Y", sp, :err_member))
        metrics = (mean_abs_bias=:lower,)
        r = GGA.wgms_rank_restricted([X, Y]; min_n=1, metrics)
        band = r[r.source .=== :band, :]
        pooled, anchored = band[band.restriction .== "all", :], band[band.restriction .== "stake+band_in_stake_range", :]
        @test pooled.label == ["Y", "X"] && pooled.rank_mean_abs_bias == [1, 2]     # X: (0.1 + 0.5) / 2 = 0.3; Y: 0.25
        @test anchored.label == ["X", "Y"] && anchored.rank_mean_abs_bias == [1, 2]  # 0.1 against 0.3
        @test anchored.n == [4, 4] && pooled.n == [8, 8]
        @test sort(unique(r.restriction)) == ["all", "stake+band_in_stake_range"]
        pt = r[r.source .=== :point, :]
        @test Set(pt.restriction) == Set(["all", "stake+band_in_stake_range"])

        # A configuration that matched fewer records is compared on the common ones.
        Ytrim = (; matched=Y.matched[Y.matched.glacier_id .!= 2, :])
        c = GGA.wgms_rank_restricted([X, Ytrim]; min_n=1, metrics)
        cb = c[(c.source .=== :band) .& (c.restriction .== "all"), :]
        @test cb.n == [4, 4] && cb.label == ["X", "Y"]

        # Cards of several restrictions rank within each restriction.
        cards = DataFrame(label=["a", "b", "a", "b"], restriction=["all", "all", "r", "r"], source=:band,
                          season=:winter, mean_abs_bias=[0.1, 0.2, 0.4, 0.3])
        rk = GGA.wgms_rank(cards; metrics)
        @test rk.label == ["a", "b", "b", "a"] && rk.restriction == ["all", "all", "r", "r"]
        @test rk.rank_mean_abs_bias == [1, 2, 1, 2]

        @test_throws "needs at least two configurations" GGA.wgms_rank_restricted([X])
        @test_throws "labels are not distinct" GGA.wgms_rank_restricted([X, X])
        @test_throws "share no matched record" GGA.wgms_rank_restricted(
            [X, (; matched=transform(Y.matched, :record_id => (x -> "z" .* x) => :record_id))])
        @test_throws "unknown provenance class(es) :nope" GGA.wgms_rank_restricted([X, Y]; restrictions=((:nope,),))
    end
end

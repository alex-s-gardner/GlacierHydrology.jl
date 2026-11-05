"""
Units of the regional GEMB components, end to end: GEMB layers in km³ of ice equivalent per geotile ->
`individual_geotile_synthesis_gembfit_dv` -> `_gembfit_dv` file -> `runs2rgi` -> `rgi_trends` ->
`region_fit_ref_and_err` -> `region_fits_table` and `mass_closure_residual`.

Every expected value is a BigFloat sum over the hand-set per-geotile rates.
"""

using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using DimensionalData
import DimensionalData as DD
using DataFrames
using Dates

@testset "Regional GEMB components are reported in Gt" begin
    # One geotile per RGI region 1-19, plus a 20th in region 1 that shares geotile 1's group.
    n = 20
    signed3(x) = (x < 0 ? "-" : "+") * lpad(abs(x), 3, '0')
    lon0 = [-180 + 2k for k in 0:n-1]
    ids = ["lat[+60+62]lon[$(signed3(l))$(signed3(l + 2))]" for l in lon0]
    rgi = [1:19; 1]
    group = [1:19; 1]
    k = 0:n-1

    # Per-geotile rates at the fitted forcing node, km³ ice equivalent yr⁻¹ (fac: km³ yr⁻¹), and
    # discharge in Gt yr⁻¹. Exact binary fractions, so the Float64 inputs equal their BigFloat values.
    rate = Dict(
        "acc" => 2.0 .+ 0.125 .* k,
        "melt" => 1.5 .+ 0.0625 .* k,
        "refreeze" => 0.25 .+ 0.015625 .* k,
        "ec" => 0.0625 .+ 0.0078125 .* k,
        "rain" => fill(0.1875, n),
        "fac" => -0.03125 .* (k .+ 1),
    )
    discharge_gtyr = 0.0625 .* k
    area_total = 10.0 .+ k

    ddate = Dim{:date}(Date(2000, 1, 15):Month(1):Date(2004, 12, 15))
    Δt = GGA.decimalyear.(ddate) .- GGA.decimalyear(first(ddate))
    dgeotile = Dim{:geotile}(ids)
    dheight = Dim{:height}(1050.0:100.0:1450.0)
    dpscale = Dim{:pscale}([1.0, 2.0])
    dΔT = Dim{:ΔT}([0.0, 1.0])
    pscale_fit = 2.0

    # Layers scale with pscale, so reading the wrong forcing node is visible.
    layer(r) = DimArray([r[g] * p * t for g in eachindex(ids), t in Δt, p in dpscale, _ in dΔT],
                        (dgeotile, ddate, dpscale, dΔT))
    gemb = DimStack((; (Symbol(v) => layer(rate[v]) for v in ("acc", "melt", "refreeze", "ec", "rain", "fac"))...))
    gemb = GGA.gemb_add_derived_vars!(gemb)
    gemb = DD.maplayers(A -> DD.rebuild(A; metadata=Dict("surface_mask" => "glacier")), gemb)

    centers_lon = lon0 .+ 1.0
    discharge = GGA.discharge2geotile(DataFrame(; longitude=centers_lon, latitude=fill(61.0, n),
                                                discharge_gtyr), dgeotile; mass2volume=true)
    DD.metadata(discharge)["surface_mask"] = "glacier"
    GGA.dv_adjust4discharge!(gemb, discharge)

    area_km2 = DimArray(repeat(area_total ./ length(dheight), 1, length(dheight)), (dgeotile, dheight))
    gemb_fit = DataFrame(; id=ids, pscale=fill(pscale_fit, n), ΔT=zeros(n), group)
    geotiles0 = DataFrame(id=ids)
    for i in 1:19
        geotiles0[!, "rgi$i"] = Float64.(rgi .== i)
    end

    mktempdir() do dir
        run = joinpath(dir, "glacier_dh_best_nmad5_v01_filled_ac_p1_synthesized.jld2")
        GGA.FileIO.save(run, Dict("dh_hyps" => zeros(dgeotile, ddate, dheight)))
        out = GGA.individual_geotile_synthesis_gembfit_dv(run, gemb, gemb_fit, discharge, geotiles0, area_km2)
        GGA.FileIO.save(replace(run, ".jld2" => "_gembfit_dv.jld2"), Dict("geotiles" => out))

        # The run has no physics members, so its synthesis is beside it.
        regional = GGA.runs2rgi([run]; physics_member=nothing)

        # Hand values. km³ ice equivalent -> Gt is ×917/1000.
        ρ = big"917"
        volume2mass = ρ / 1000
        members = Dict(r => findall(==(r), rgi) for r in 1:19)
        members[98] = findall(in(13:15), rgi)
        members[99] = collect(eachindex(rgi))
        km3(v, r) = pscale_fit * sum(big.(rate[v][members[r]]))   # regional rate, km³ i.e. yr⁻¹
        D(r) = sum(big.(discharge_gtyr[members[r]]))               # regional discharge, Gt yr⁻¹
        runoff_km3(r) = km3("melt", r) - km3("refreeze", r)
        smb_km3(r) = km3("acc", r) - km3("melt", r) + km3("refreeze", r) - km3("ec", r)
        expected_rate = Dict(
            "acc" => r -> km3("acc", r) * volume2mass,
            "melt" => r -> km3("melt", r) * volume2mass,
            "refreeze" => r -> km3("refreeze", r) * volume2mass,
            "ec" => r -> km3("ec", r) * volume2mass,
            "rain" => r -> km3("rain", r) * volume2mass,
            "runoff" => r -> runoff_km3(r) * volume2mass,
            "smb" => r -> smb_km3(r) * volume2mass,
            "discharge" => D,
            "dm" => r -> smb_km3(r) * volume2mass - D(r),
            "dm_altim" => r -> -km3("fac", r) * volume2mass,   # altimetry dv is zero here
            "fac" => r -> km3("fac", r),
            "dv" => r -> smb_km3(r) + km3("fac", r) - D(r) * 1000 / ρ,
        )
        series(v, r) = regional[v][run=At(run), rgi=At(r)]

        @testset "runs2rgi: masses in Gt, volumes in km³ ($v)" for v in keys(expected_rate)
            @test DD.metadata(regional[v])["units"] == (v in ("fac", "dv") ? "km3" : "Gt")
            for r in (1, 2, 13, 19, 98, 99)
                @test series(v, r) ≈ Float64.(expected_rate[v](r) .* big.(Δt)) rtol = 1e-12
            end
        end

        @testset "runs2rgi: a GEMB mass component is the km³ sum × δice/1000" begin
            # The same sum, built from the geotile columns of the file, which are km³ i.e.
            @test all(colmetadata(out, v, "units") == "km3 [i.e.]" for v in ("acc", "runoff", "smb", "discharge"))
            @test colmetadata(out, "dm", "units") == "Gt"
            for v in ("acc", "melt", "refreeze", "ec", "rain", "runoff", "smb")
                @test series(v, 99) ≈ sum(out[!, v]) .* (GGA.δice / 1000) rtol = 1e-12
            end
        end

        @testset "runs2rgi: the mass budget closes in Gt at every date" begin
            for r in (1, 7, 13, 98, 99)
                @test series("acc", r) .- series("runoff", r) .- series("ec", r) ≈ series("smb", r) rtol = 1e-12
                @test series("smb", r) .- series("discharge", r) ≈ series("dm", r) rtol = 1e-12
            end
        end

        drgi = Dim{:rgi}([1:19; 98; 99])
        derror = Dim{:error}([false, true])
        discharge_rgi = zeros(drgi, derror)
        for r in val(drgi)
            discharge_rgi[rgi=At(r), error=At(false)] = Float64(D(r))
            discharge_rgi[rgi=At(r), error=At(true)] = 0.15 * Float64(D(r))
        end

        fits = GGA.rgi_trends(regional, discharge_rgi, (first(ddate), last(ddate)))
        region_fits = GGA.region_fit_ref_and_err(fits, run; discharge=discharge_rgi)
        trend(v, r) = region_fits[varname=At(v), rgi=At(r), parameter=At("trend"), error=At(false)]

        @testset "trends: Gt yr⁻¹, closure, freshwater flux and GSI" begin
            for v in keys(expected_rate), r in (1, 98, 99)
                @test trend(v, r) ≈ Float64(expected_rate[v](r)) rtol = 1e-8
            end
            for r in (1, 5, 98, 99)
                @test abs(GGA.mass_closure_residual(region_fits; rgi=r)) < 1e-8
            end
            # Freshwater to rivers and oceans: runoff plus discharge, both Gt yr⁻¹.
            freshwater = trend("runoff", 99) + discharge_rgi[rgi=At(99), error=At(false)]
            @test freshwater ≈ Float64(runoff_km3(99) * volume2mass + D(99)) rtol = 1e-8
            net_acc = (km3("acc", 99) - km3("ec", 99)) * volume2mass
            @test trend("net_acc", 99) ≈ Float64(net_acc) rtol = 1e-8
            @test trend("gsi", 99) ≈ Float64((runoff_km3(99) * volume2mass + D(99)) / net_acc) rtol = 1e-8
        end

        @testset "region_fits_table labels masses Gt and volumes km³" begin
            table = GGA.region_fits_table(region_fits; varnames=["runoff", "acc", "smb", "dm", "fac", "dv"],
                                          rgi_regions=[1, 99], digits=4)
            row = findfirst(==(99), table.rgi)
            value(label) = parse(Float64, first(split(table[row, label], " ± ")))
            @test value("runoff_trend_[Gt/yr]") == round(Float64(runoff_km3(99) * volume2mass); digits=4)
            @test value("acc_trend_[Gt/yr]") == round(Float64(km3("acc", 99) * volume2mass); digits=4)
            @test value("smb_trend_[Gt/yr]") == round(Float64(smb_km3(99) * volume2mass); digits=4)
            @test value("fac_trend_[km³/yr]") == round(Float64(km3("fac", 99)); digits=4)
            @test "dv_trend_[km³/yr]" in names(table)
        end

        @testset "runs2rgi throws on an unexpected units tag" begin
            bad = joinpath(dir, "glacier_dh_best_nmad3_v01_filled_ac_p1_synthesized.jld2")
            out_bad = copy(out)
            colmetadata!(out_bad, "runoff", "units", "m", style=:note)
            GGA.FileIO.save(replace(bad, ".jld2" => "_gembfit_dv.jld2"), Dict("geotiles" => out_bad))
            @test_throws "mass column runoff has units \"m\"" GGA.runs2rgi([bad]; physics_member=nothing)
        end

        @testset "runs2rgi_physics: GEMB physics members as a dimension beside the runs" begin
            # Two members' syntheses of the same run, the second with twice the accumulation.
            members = ["B-kfit-table-w2-g2.5", "C-kfit-table-w2-g2.5"]
            twice = deepcopy(out)
            foreach(x -> x .*= 2, twice.acc)
            for (m, synthesis) in zip(members, (out, twice))
                file = GGA.gemb_physics_path(replace(run, ".jld2" => "_gembfit_dv.jld2"), m)
                mkpath(dirname(file))
                GGA.FileIO.save(file, Dict("geotiles" => synthesis))
            end
            # By default the reference member is read, which is the first here.
            @test members[1] == GGA.gemb_reference_physics_member
            @test GGA.runs2rgi([run])["acc"] == regional["acc"]
            both = GGA.runs2rgi_physics([run], members; vars2sum=["acc", "melt", "dv"])
            @test sort(collect(keys(both))) == ["acc", "dv", "melt"]
            A = both["acc"]
            @test dims(A, :physics) == Dim{:physics}(members)
            @test DD.name.(dims(A)) == (:run, :physics, :rgi, :date)
            @test DD.metadata(A)["units"] == "Gt"
            @test DD.metadata(A)["gemb_physics_members"] == members
            @test A[physics=At(members[1])] == regional["acc"]
            @test A[physics=At(members[2])] ≈ 2 .* regional["acc"]
            @test both["melt"][physics=At(members[2])] == regional["melt"]
            @test_throws "physics members repeat" GGA.runs2rgi_physics([run], [members[1], members[1]])
            @test_throws "is not a GEMB physics member id" GGA.runs2rgi_physics([run], ["A_kfit_table"])
        end
    end
end

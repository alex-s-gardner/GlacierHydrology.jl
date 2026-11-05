"""
Tests for reading GEMB_GlacierSims tile NetCDF (utilities_gemb_tiles.jl).

Expected values are computed from the fixture arrays in BigFloat, band by band, independently of the
reader's matrix-vector formulation.
"""

using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using GlobalGlacierAnalysis.NCDatasets
using DimensionalData
using Dates
using Random
using DataFrames
using TOML

# Write a one-forcing-node tile file, at temperature offset `ΔT` and precipitation scaling `pscale`.
# Fluxes are (time, band) in kg m⁻², per interval; `firn` is the non-mass height change dh - dh_mass in
# metres.
function write_test_tile(path, times, band_centers; P, Rn, Ro, Rf, Mt, EC, firn, ρi=917.0,
                         geotile_id="lat[+60+62]lon[-142-140]", run_parameters=Dict{String,Any}(),
                         ΔT=0.0, pscale=1.0)
    n_time, n_band = size(P)
    dh_mass = cumsum(P .- Ro .+ EC; dims=1) ./ ρi
    dh = dh_mass .+ firn
    NCDataset(path, "c") do ds
        defDim(ds, "time", n_time)
        defDim(ds, "band", n_band)
        defDim(ds, "delta_temperature", 1)
        defDim(ds, "precipitation_scaling", 1)
        ds.attrib["geotile_id"] = geotile_id
        ds.attrib["model_density_ice"] = ρi
        # GEMB_GlacierSims' run-parameter convention: a roster attribute naming one attribute per value.
        for (k, v) in run_parameters
            ds.attrib[k] = v
        end
        isempty(run_parameters) || (ds.attrib["run_parameters"] = join(sort(collect(keys(run_parameters))), " "))
        t = defVar(ds, "time", Float64, ("time",); attrib=Dict("units" => "days since 1950-01-01 00:00:00"))
        t[:] = times
        defVar(ds, "band_center", band_centers, ("band",))
        defVar(ds, "delta_temperature", [Float64(ΔT)], ("delta_temperature",))
        defVar(ds, "precipitation_scaling", [Float64(pscale)], ("precipitation_scaling",))
        for (name, a) in ("dh" => dh, "dh_mass" => dh_mass, "precipitation" => P, "rain" => Rn,
                          "runoff" => Ro, "refreeze" => Rf, "melt" => Mt, "evaporation_condensation" => EC)
            defVar(ds, name, Float32.(reshape(a, n_time, n_band, 1, 1)),
                   ("time", "band", "delta_temperature", "precipitation_scaling"))
        end
    end
    return path
end

# Hypsometry on the project height bins with `km2` at the given bin centers and zero elsewhere.
function test_hypsometry(km2_at)
    _, height_centers = GGA.project_height_bins()
    area = DimArray(zeros(length(height_centers)), Dim{:height}(collect(height_centers)))
    for (h, a) in km2_at
        area[height=At(h)] = a
    end
    return area
end

days_since_epoch(times) = Dates.value.(times .- DateTime(1950)) ./ 86_400_000

# Cumulative series referenced to its first element, the reader's output convention.
referenced(per_interval) = (c = cumsum(per_interval); c .- first(c))

@testset "GEMB tile reader" begin
    @testset "rain is split from runoff per band, before the area weighting" begin
        times = [DateTime(2000, 1, 31, 23), DateTime(2000, 2, 29, 23), DateTime(2000, 3, 31, 23)]
        band_centers = [1050.0, 1150.0]
        # (time, band), kg m⁻² per interval. In interval 3 band 1 keeps 20 of its 30 rain and band 2
        # sheds all its 20 rain plus 5 of meltwater.
        P  = [100.0 80.0; 50.0 40.0; 80.0 60.0]
        Rn = [10.0 0.0; 0.0 0.0; 30.0 20.0]
        Ro = [20.0 15.0; 0.0 0.0; 10.0 25.0]
        Rf = [5.0 3.0; 0.0 0.0; 5.0 4.0]
        Mt = [15.0 18.0; 0.0 0.0; 0.0 9.0]
        EC = [-1.0 -2.0; 2.0 1.0; -3.0 0.0]
        firn = [0.01 0.02; 0.03 0.01; -0.02 0.04]

        # 10 km² at band 1, 30 km² at band 2, and 5 km² at 1250 m, beyond the bands and held at band 2.
        area = test_hypsometry((1050.0 => 10.0, 1150.0 => 30.0, 1250.0 => 5.0))

        path = write_test_tile(joinpath(mktempdir(), "tile.nc"), times, band_centers; P, Rn, Ro, Rf, Mt, EC, firn)
        tile = GGA.read_gemb_tile(path, area, days_since_epoch(times); n_pscale=1, n_ΔT=1)
        L = Dict(k => vec(v) for (k, v) in tile.layers)

        # Hand calculation in BigFloat: band weights 10 and 35 km² = 0.010 and 0.035 km³ per metre.
        ρ = big"917"
        w = [big"10", big"35"] .* big"1e-3"
        b(a) = big.(a)
        retained = max.(b(Rn) .- b(Ro), 0)
        melt_runoff = max.(b(Ro) .- b(Rn), 0)
        acc_hand = referenced((b(P) .- b(Rn) .+ retained) * w ./ ρ)
        runoff_hand = referenced(melt_runoff * w ./ ρ)
        rain_hand = referenced((b(Rn) .- retained) * w ./ ρ)
        refreeze_hand = referenced(b(Rf) * w ./ ρ)
        ec_hand = -referenced(b(EC) * w ./ ρ)
        smb_hand = referenced((b(P) .- b(Ro) .+ b(EC)) * w ./ ρ)
        fac_hand = (b(firn) * w) .- first(b(firn) * w)
        dv_hand = smb_hand .+ fac_hand

        @test tile.extrapolated_fraction ≈ 5 / 45
        @test L[:acc] ≈ Float64.(acc_hand) rtol = 1e-6
        @test L[:melt] .- L[:refreeze] ≈ Float64.(runoff_hand) rtol = 1e-6
        @test L[:rain] ≈ Float64.(rain_hand) rtol = 1e-6
        @test L[:refreeze] ≈ Float64.(refreeze_hand) rtol = 1e-6
        @test L[:ec] ≈ Float64.(ec_hand) rtol = 1e-6
        @test L[:dv_mass] ≈ Float64.(smb_hand) rtol = 1e-6
        @test L[:dv] ≈ Float64.(dv_hand) rtol = 1e-6
        @test L[:fac] ≈ Float64.(fac_hand) rtol = 1e-5
        @test L[:acc] .- L[:melt] .+ L[:refreeze] .- L[:ec] ≈ Float64.(smb_hand) rtol = 1e-6
        @test tile.closure_km3 < 1e-8

        # Over intervals 2-3 the only meltwater runoff is band 2's 5 kg m⁻² in interval 3:
        # 5 × 0.035 / 917 km³. Accumulation is (0.5 + 1.4) + (0.7 + 1.4) = 4.0 km³ m⁻¹·kg m⁻² over 917.
        @test last(runoff_hand) ≈ 5 * big"0.035" / ρ rtol = 1e-60
        @test last(acc_hand) ≈ big"4.0" / ρ rtol = 1e-60
        @test round(last(L[:melt] .- L[:refreeze]); digits=6) == 0.000191
        @test round(last(L[:acc]); digits=6) == 0.004362

        # Splitting the area-weighted tile totals gives no meltwater runoff at all here
        # (rain 1.0 > runoff 0.975 km³ m⁻¹·kg m⁻² in interval 3); the reader must not do that.
        Rt, Rot = b(Rn) * w, b(Ro) * w
        @test last(referenced(max.(Rot .- Rt, 0) ./ ρ)) == 0
        @test last(L[:melt] .- L[:refreeze]) > 1e-4
    end

    @testset "per-band split: SMB, dv and fac do not depend on it; acc + rain is the precipitation" begin
        rng = Xoshiro(20261015)
        n_time, n_band = 12, 5
        times = [DateTime(2000, 1, 1) + Day(7k) for k in 1:n_time]
        band_centers = collect(1050.0:100.0:1450.0)
        Rn = 40 .* rand(rng, n_time, n_band)
        P = Rn .+ 60 .* rand(rng, n_time, n_band)   # precipitation includes the rain
        Ro = 60 .* rand(rng, n_time, n_band)
        Rf = 10 .* rand(rng, n_time, n_band)
        Mt = 50 .* rand(rng, n_time, n_band)
        EC = 4 .* rand(rng, n_time, n_band) .- 2
        firn = 0.05 .* randn(rng, n_time, n_band)
        area = test_hypsometry(h => a for (h, a) in zip(band_centers, (3.0, 8.0, 12.0, 6.0, 1.0)))

        path = write_test_tile(joinpath(mktempdir(), "tile.nc"), times, band_centers; P, Rn, Ro, Rf, Mt, EC, firn)
        tile = GGA.read_gemb_tile(path, area, days_since_epoch(times); n_pscale=1, n_ΔT=1)
        L = Dict(k => vec(v) for (k, v) in tile.layers)

        ρ = big"917"
        w = big.([3.0, 8.0, 12.0, 6.0, 1.0]) .* big"1e-3"
        b(a) = big.(Float32.(a))   # the file stores Float32
        smb_hand = referenced((b(P) .- b(Ro) .+ b(EC)) * w ./ ρ)
        precipitation_hand = referenced(b(P) * w ./ ρ)
        runoff_tile_split = referenced(max.(b(Ro) * w .- b(Rn) * w, 0) ./ ρ)

        @test L[:acc] .- L[:melt] .+ L[:refreeze] .- L[:ec] ≈ Float64.(smb_hand) rtol = 1e-6
        @test L[:dv_mass] ≈ Float64.(smb_hand) rtol = 1e-6
        @test L[:fac] ≈ L[:dv] .- L[:dv_mass] rtol = 1e-12
        @test L[:acc] .+ L[:rain] ≈ Float64.(precipitation_hand) rtol = 1e-6
        # Σ max(x_b, 0) ≥ max(Σ x_b, 0): the per-band meltwater runoff is never below the tile-level one.
        @test all(L[:melt] .- L[:refreeze] .>= Float64.(runoff_tile_split) .- 1e-12)
        @test all(diff(L[:melt] .- L[:refreeze]) .>= 0)
        @test all(diff(L[:rain]) .>= 0)
        @test all(diff(L[:acc]) .>= 0)
    end

    @testset "tile provenance: read from the tiles, written beside the calibration, round-tripped" begin
        times = [DateTime(2000, 1, 31, 23), DateTime(2000, 2, 29, 23), DateTime(2000, 3, 31, 23)]
        z = zeros(3, 1)
        flux = (; P=z, Rn=z, Ro=z, Rf=z, Mt=z, EC=z, firn=z)
        software = Dict{String,Any}(
            "software_gemb_tree" => "1111111111111111111111111111111111111111",
            "software_gemb_climateforcing_tree" => "2222222222222222222222222222222222222222",
            "software_gemb_glaciersims_tree" => "3333333333333333333333333333333333333333+dirty:0a1b2c3d",
            "software_gemb_version" => "0.9.0")
        ids = ["lat[+60+62]lon[-142-140]", "lat[+60+62]lon[-140-138]"]
        # A model parameter resolved per tile differs between the two; it is not sweep-wide.
        tile_parameters = [merge(software, Dict{String,Any}("model_albedo_ice" => a)) for a in (0.40, 0.45)]

        sweep = mkpath(joinpath(mktempdir(), "tile_runs_v2"))
        paths = [write_test_tile(joinpath(sweep, "tile_$i.nc"), times, [1050.0]; flux...,
                                 geotile_id=ids[i], run_parameters=tile_parameters[i]) for i in 1:2]

        provenance = GGA.gemb_tile_provenance(sweep, paths)
        @test provenance["tile_dir"] == realpath(sweep)
        @test provenance["tile_run_parameters"] == software

        # process_gemb_tiles records the provenance of the files it read.
        _, height_centers = GGA.project_height_bins()
        dgeotile = Dim{:geotile}(ids)
        area = DimArray(zeros(2, length(height_centers)), (dgeotile, Dim{:height}(collect(height_centers))))
        area[:, findfirst(==(1050.0), height_centers)] .= 10.0
        gemb = GGA.process_gemb_tiles(DataFrame(id=ids), area; tile_dir=sweep,
                                      date_center=[Date(2000, 2, 15), Date(2000, 3, 15)], show_stats=false)
        @test DimensionalData.metadata(gemb)["tile_provenance"] == provenance

        # The calibration output's sidecar, in GEMB_GlacierSims' format.
        gembfit = joinpath(mktempdir(), "glacier_rgi7_dh_cop30_v2_cc_nmad5_v01_filled_ac_p6_synthesized_gembfit.arrow")
        sidecar = GGA.write_gembfit_provenance(gembfit, provenance)
        @test sidecar == replace(gembfit, ".arrow" => "_provenance.toml")
        @test sidecar == GGA.gembfit_provenance_file(gembfit)
        stored = TOML.parsefile(sidecar)
        @test stored == Dict{String,Any}("tile_dir" => realpath(sweep), "tile_run_parameters" => software)

        # What GEMB_GlacierSims' `assert_corrections_fitted_on` checks, against each tile's own record
        # (read here straight from the attributes): the directory, the three trees, and every recorded key.
        @test realpath(stored["tile_dir"]) == realpath(sweep)
        @test all(haskey(stored["tile_run_parameters"], k) for k in GGA.GEMB_TILE_PROVENANCE_REQUIRED)
        for path in paths
            attributes = GGA.NCDatasets.NCDataset(path, "r") do ds
                Dict(k => ds.attrib[k] for k in split(ds.attrib["run_parameters"]))
            end
            @test all(attributes[k] == v for (k, v) in stored["tile_run_parameters"])
        end

        @test_throws "not an .arrow calibration output" GGA.gembfit_provenance_file(replace(gembfit, ".arrow" => ".jld2"))

        # Tiles of two sweeps, or a tile that does not record its code, are refused before any read.
        mixed = mkpath(joinpath(mktempdir(), "tile_runs_mixed"))
        other = merge(software, Dict{String,Any}("software_gemb_tree" => "4444444444444444444444444444444444444444"))
        mixed_paths = [write_test_tile(joinpath(mixed, "tile_$i.nc"), times, [1050.0]; flux...,
                                       geotile_id=ids[i], run_parameters=(software, other)[i]) for i in 1:2]
        @test_throws "GEMB tile files disagree on software_gemb_tree" GGA.gemb_tile_provenance(mixed, mixed_paths)
        @test_throws "GEMB tile files disagree on software_gemb_tree" GGA.process_gemb_tiles(
            DataFrame(id=ids), area; tile_dir=mixed, date_center=[Date(2000, 2, 15)], show_stats=false)

        untraced = mkpath(joinpath(mktempdir(), "tile_runs_untraced"))
        partial = Dict{String,Any}(k => v for (k, v) in software if k != "software_gemb_glaciersims_tree")
        untraced_path = write_test_tile(joinpath(untraced, "tile.nc"), times, [1050.0]; flux...,
                                        run_parameters=partial)
        @test_throws "records no software_gemb_glaciersims_tree" GGA.gemb_tile_provenance(untraced, [untraced_path])
        @test_throws "tile files outside" GGA.gemb_tile_provenance(sweep, [untraced_path])
    end
end

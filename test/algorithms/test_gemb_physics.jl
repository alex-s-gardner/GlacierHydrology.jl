"""
Tests for GEMB physics members (utilities_gemb_physics.jl and the member checks in the tile reader,
ensemble and calibration): member ids, the per-member output folder, and the refusal to mix tiles or
fits of different members.

Uses `write_test_tile` and `test_hypsometry` from test_gemb_tiles.jl.
"""

using Test
using GlacierHydrology
import GlacierHydrology
using DimensionalData
using Dates
using DataFrames
using GlacierHydrology.NCDatasets
using TOML

const _GRID_TILE_SOFTWARE = Dict{String,Any}(
    "software_gemb_tree" => "1111111111111111111111111111111111111111",
    "software_gemb_climateforcing_tree" => "2222222222222222222222222222222222222222",
    "software_gemb_glaciersims_tree" => "3333333333333333333333333333333333333333")

# A tile on a (ΔT, pscale) grid: per interval, at node (ΔT, pscale), precipitation `P(t, p)` and runoff
# `Ro(t, p)` [kg m⁻²] in every band, everything else zero; `run_parameters` beside the software trees.
function write_grid_tile(path, id, dts, pss; times, bands=[1050.0], P=(t, p) -> 100.0 * p,
                         Ro=(t, p) -> 10.0 * (1 + t), run_parameters=Dict{String,Any}())
    n, nb = length(times), length(bands)
    Pa = [P(t, p) for _ in 1:n, _ in 1:nb, t in dts, p in pss]
    Ra = [Ro(t, p) for _ in 1:n, _ in 1:nb, t in dts, p in pss]
    z = zeros(size(Pa))
    dh_mass = cumsum(Pa .- Ra; dims=1) ./ 917.0
    attributes = merge(_GRID_TILE_SOFTWARE, run_parameters)
    NCDataset(path, "c") do ds
        defDim(ds, "time", n); defDim(ds, "band", nb)
        defDim(ds, "delta_temperature", length(dts)); defDim(ds, "precipitation_scaling", length(pss))
        ds.attrib["geotile_id"] = id
        ds.attrib["model_density_ice"] = 917.0
        for (k, v) in attributes
            ds.attrib[k] = v
        end
        ds.attrib["run_parameters"] = join(sort(collect(keys(attributes))), " ")
        t = defVar(ds, "time", Float64, ("time",); attrib=Dict("units" => "days since 1950-01-01 00:00:00"))
        t[:] = times
        defVar(ds, "band_center", Float64.(bands), ("band",))
        defVar(ds, "delta_temperature", Float64.(dts), ("delta_temperature",))
        defVar(ds, "precipitation_scaling", Float64.(pss), ("precipitation_scaling",))
        for (name, a) in ("dh" => dh_mass, "dh_mass" => dh_mass, "precipitation" => Pa, "rain" => z,
                          "runoff" => Ra, "refreeze" => z, "melt" => Ra, "evaporation_condensation" => z)
            defVar(ds, name, Float32.(a), ("time", "band", "delta_temperature", "precipitation_scaling"))
        end
    end
    return path
end

@testset "GEMB physics members" begin
    @testset "member ids and folders" begin
        @test GlacierHydrology.check_gemb_physics_member("C-kfit-table") == "C-kfit-table"
        @test GlacierHydrology.check_gemb_physics_member(GlacierHydrology.gemb_reference_physics_member) == "S-kfit-table"
        @test GlacierHydrology.gemb_physics_member_parts(GlacierHydrology.gemb_reference_physics_member) ==
              (; tile_member="S-kfit-table", precipitation_gradient=0.0)
        for bad in ("C_ks35_table", "C-ks35", "C/ks35-table", "", :C, "C-kfit-table-g2.5-g5",
                    "C-kfit-g2.5", "C-g2-table")
            @test_throws "is not a GEMB physics member id" GlacierHydrology.check_gemb_physics_member(bad)
        end
        # A precipitation gradient is a fourth token, in one spelling, and never zero.
        @test GlacierHydrology.check_gemb_physics_member("B-kshaw-table-g2.5") == "B-kshaw-table-g2.5"
        @test GlacierHydrology.gemb_physics_member_parts("B-kshaw-table-g2.5") == (; tile_member="B-kshaw-table", precipitation_gradient=2.5)
        @test GlacierHydrology.gemb_physics_member_parts("B-kfit-table") == (; tile_member="B-kfit-table", precipitation_gradient=0.0)
        @test GlacierHydrology.gemb_physics_member_id("B-kfit-table", 5) == "B-kfit-table-g5"
        @test GlacierHydrology.gemb_physics_member_id("B-kfit-table", 0) == "B-kfit-table"
        @test_throws "the gradient-free member is \"B-kfit-table\"" GlacierHydrology.check_gemb_physics_member("B-kfit-table-g0")
        @test_throws "write the gradient as \"g5\"" GlacierHydrology.check_gemb_physics_member("B-kfit-table-g5.0")
        @test_throws "already carries a precipitation gradient" GlacierHydrology.gemb_physics_member_id("B-kfit-table-g5", 2.5)
        # A further GEMB_GlacierSims level (a forcing correction, say) is a fourth tile-sweep token, and a
        # gradient still follows it.
        @test GlacierHydrology.gemb_physics_member_parts("D1-kfit-table-w2") == (; tile_member="D1-kfit-table-w2", precipitation_gradient=0.0)
        @test GlacierHydrology.gemb_physics_member_parts("D1-kfit-table-w2-g2.5") == (; tile_member="D1-kfit-table-w2", precipitation_gradient=2.5)
        @test GlacierHydrology.gemb_physics_member_id("D1-kfit-table-w2", 5) == "D1-kfit-table-w2-g5"

        root = joinpath(GlacierHydrology.pathlocal[:data_dir], "binned", "2deg")
        file = joinpath(root, "glacier_rgi7_dh_cop30_v2_cc_nmad5_v01_filled_ac_p6_synthesized_gembfit.arrow")
        moved = GlacierHydrology.gemb_physics_path(file, "C-kfit-table")
        @test moved == joinpath(root, "gemb_physics", "C-kfit-table", basename(file))
        synthesized = replace(file, "_gembfit.arrow" => ".jld2")
        @test GlacierHydrology.gembfit_path(synthesized, "C-kfit-table") == moved
        @test GlacierHydrology.gembfit_dv_path(synthesized, "C-kfit-table") ==
              joinpath(root, "gemb_physics", "C-kfit-table", replace(basename(synthesized), ".jld2" => "_gembfit_dv.jld2"))
        # The altimetry member's file-name rules see the same name in a member folder.
        @test GlacierHydrology.binned_filled_fileparts(moved) == GlacierHydrology.binned_filled_fileparts(file)
        @test_throws "is not a GEMB physics member id" GlacierHydrology.gemb_physics_path(file, "C_ks35_table")
    end

    @testset "tiles of different members are not mixed" begin
        times = [DateTime(2000, 1, 31, 23), DateTime(2000, 2, 29, 23), DateTime(2000, 3, 31, 23)]
        z = zeros(3, 1)
        flux = (; P=z, Rn=z, Ro=z, Rf=z, Mt=z, EC=z, firn=z)
        software = Dict{String,Any}(
            "software_gemb_tree" => "1111111111111111111111111111111111111111",
            "software_gemb_climateforcing_tree" => "2222222222222222222222222222222222222222",
            "software_gemb_glaciersims_tree" => "3333333333333333333333333333333333333333")
        member(m) = Dict{String,Any}("physics_member" => m, "physics_grain" => string(first(m)),
                                     "physics_band_decoupling" => split(m, '-')[2],
                                     "physics_impurity" => split(m, '-')[3])
        ids = ["lat[+60+62]lon[-142-140]", "lat[+60+62]lon[-140-138]"]
        write_sweep(dir, members; kwargs...) = [
            write_test_tile(joinpath(dir, "tile_$i.nc"), times, [1050.0]; flux..., geotile_id=ids[i],
                            run_parameters=merge(software, member(members[i])), kwargs...) for i in 1:2]

        sweep = mkpath(joinpath(mktempdir(), "tile_runs_v2", "C-kfit-table"))
        paths = write_sweep(sweep, ["C-kfit-table", "C-kfit-table"])
        tiles = collect(values(GlacierHydrology.gemb_tile_coverage(sweep)))
        p = GlacierHydrology.gemb_tile_provenance(sweep, tiles; physics_member="C-kfit-table")
        @test p["tile_run_parameters"] == merge(software, member("C-kfit-table"))
        @test_throws "tiles of GEMB physics member \"C-kfit-table\", not \"B-kfit-table\"; members must not be mixed" GlacierHydrology.gemb_tile_provenance(
            sweep, tiles; physics_member="B-kfit-table")

        # A WGMS score pairs tiles with fits only of the member it names.
        fits = DataFrame(id=ids, pscale=[1.0, 1.0], ΔT=[0.0, 0.0])
        @test GlacierHydrology.wgms_tile_model(sweep, fits; physics_member="C-kfit-table") isa Function
        @test_throws "tiles of GEMB physics member \"C-kfit-table\", not \"B-kfit-table\"" GlacierHydrology.wgms_tile_model(
            sweep, fits; physics_member="B-kfit-table")

        # The calibration sidecar records the member beside the code, so GEMB_GlacierSims'
        # `assert_corrections_fitted_on` can compare it with the reference sweep's tiles.
        gembfit = GlacierHydrology.gemb_physics_path(joinpath(mktempdir(), "glacier_rgi7_dh_cop30_v2_cc_nmad5_v01_filled_ac_p6_synthesized_gembfit.arrow"),
                                        "C-kfit-table")
        mkpath(dirname(gembfit))
        stored = TOML.parsefile(GlacierHydrology.write_gembfit_provenance(gembfit, p))
        @test stored["tile_run_parameters"]["physics_member"] == "C-kfit-table"
        @test stored["tile_run_parameters"]["physics_grain"] == "C"

        mixed = mkpath(joinpath(mktempdir(), "tile_runs_mixed"))
        mixed_paths = write_sweep(mixed, ["C-kfit-table", "B-kfit-table"])
        @test_throws "members must not be mixed" GlacierHydrology.gemb_tile_provenance(
            mixed, collect(values(GlacierHydrology.gemb_tile_coverage(mixed))); physics_member="C-kfit-table")

        # process_gemb_tiles reads only the member it is asked for.
        _, height_centers = GlacierHydrology.project_height_bins()
        area = DimArray(zeros(2, length(height_centers)), (Dim{:geotile}(ids), Dim{:height}(collect(height_centers))))
        area[:, findfirst(==(1050.0), height_centers)] .= 10.0
        dates = [Date(2000, 2, 15), Date(2000, 3, 15)]
        gemb = GlacierHydrology.process_gemb_tiles(DataFrame(id=ids), area; tile_dir=sweep, date_center=dates,
                                      show_stats=false, physics_member="C-kfit-table")
        @test DimensionalData.metadata(gemb)["tile_provenance"] == p
        @test_throws "members must not be mixed" GlacierHydrology.process_gemb_tiles(
            DataFrame(id=ids), area; tile_dir=sweep, date_center=dates, show_stats=false,
            physics_member="B-kfit-table")

        # Tiles of one sweep must share the forcing grid, value for value, not only its size.
        shifted = mkpath(joinpath(mktempdir(), "tile_runs_shifted", "C-kfit-table"))
        shifted_paths = [write_test_tile(joinpath(shifted, "tile_$i.nc"), times, [1050.0]; flux...,
                                         geotile_id=ids[i], run_parameters=merge(software, member("C-kfit-table")),
                                         ΔT=(0.0, 0.5)[i]) for i in 1:2]
        @test_throws "but the ensemble's is ΔT [0.0] × pscale [1.0]; every tile of a sweep must share it" GlacierHydrology.process_gemb_tiles(
            DataFrame(id=ids), area; tile_dir=shifted, date_center=dates, show_stats=false,
            physics_member="C-kfit-table")
    end

    @testset "the pilot ensemble reads its nodes from tiles that carry more" begin
        times = [DateTime(2000, 1, 31, 23), DateTime(2000, 2, 29, 23), DateTime(2000, 3, 31, 23)]
        ids = ["lat[+60+62]lon[-142-140]", "lat[+60+62]lon[-140-138]"]
        dir = mktempdir()
        # A dense tile, its axes out of order, beside a tile on the reduced grid alone.
        pilot = Dict{String,Any}("physics_member" => "B-kfit-table")
        write_grid_tile(joinpath(dir, "dense.nc"), ids[1], [1.0, -1.0, 0.0], [2.0, 0.5, 1.0]; times, run_parameters=pilot)
        write_grid_tile(joinpath(dir, "reduced.nc"), ids[2], [0.0, 1.0], [1.0, 2.0]; times, run_parameters=pilot)
        _, height_centers = GlacierHydrology.project_height_bins()
        area = DimArray(zeros(2, length(height_centers)), (Dim{:geotile}(ids), Dim{:height}(collect(height_centers))))
        area[:, findfirst(==(1050.0), height_centers)] .= 10.0
        dates = [Date(2000, 2, 15), Date(2000, 3, 15)]
        ensemble(grid) = GlacierHydrology.process_gemb_tiles(DataFrame(id=ids), area; tile_dir=dir, date_center=dates,
                                            show_stats=false, forcing_grid=grid, physics_member="B-kfit-table")

        gemb = ensemble((; pscale=[2.0, 1.0], ΔT=[1.0, 0.0]))
        @test val(dims(gemb, :pscale)) == [1.0, 2.0] && val(dims(gemb, :ΔT)) == [0.0, 1.0]
        # Identical fluxes at every shared node, so both geotiles read the same values there.
        @test gemb[:acc][geotile=At(ids[1])] == gemb[:acc][geotile=At(ids[2])]
        @test gemb[:melt][geotile=At(ids[1])] == gemb[:melt][geotile=At(ids[2])]
        # Each node's values come from that node: accumulation scales with pscale, runoff with 1 + ΔT.
        acc = gemb[:acc][geotile=At(ids[1])]
        melt = gemb[:melt][geotile=At(ids[1])]
        @test acc[date=2, pscale=At(2.0), ΔT=At(0.0)] ≈ 2 * acc[date=2, pscale=At(1.0), ΔT=At(0.0)] > 0
        @test melt[date=2, pscale=At(1.0), ΔT=At(1.0)] ≈ 2 * melt[date=2, pscale=At(1.0), ΔT=At(0.0)] > 0

        @test_throws "reduced.nc has no ΔT node -1.0" ensemble((; pscale=[1.0], ΔT=[-1.0, 0.0]))
        @test_throws "forcing_grid repeats a node" ensemble((; pscale=[1.0, 1.0], ΔT=[0.0]))
        # Without a grid the tiles must still agree, so mixed grids are an error.
        @test_throws "forcing grid is" ensemble(nothing)
    end

    @testset "precipitation-gradient members read each band at its own pscale" begin
        times = [DateTime(2000, 1, 31, 23), DateTime(2000, 2, 29, 23), DateTime(2000, 3, 31, 23)]
        ids = ["lat[+60+62]lon[-142-140]"]
        member(m) = Dict{String,Any}("physics_member" => m, "physics_grain" => "B",
                                     "physics_band_decoupling" => "kfit", "physics_impurity" => "table")
        dir = mkpath(joinpath(mktempdir(), "B-kfit-table"))
        # Two bands; precipitation linear in pscale, runoff 5 pscale², so interpolation in pscale is visible.
        path = write_grid_tile(joinpath(dir, "tile.nc"), ids[1], [0.0], [1.0, 2.0, 4.0]; times,
                               bands=[1050.0, 1150.0], P=(t, p) -> 100.0 * p, Ro=(t, p) -> 5.0 * p^2,
                               run_parameters=member("B-kfit-table"))
        _, height_centers = GlacierHydrology.project_height_bins()
        area = DimArray(zeros(1, length(height_centers)), (Dim{:geotile}(ids), Dim{:height}(collect(height_centers))))
        area[1, findfirst(==(1050.0), height_centers)] = 10.0
        area[1, findfirst(==(1150.0), height_centers)] = 30.0
        days = Dates.value.(times .- DateTime(1950)) ./ 86_400_000
        tile_read(g, ps) = GlacierHydrology.read_gemb_tile(path, area[geotile=At(ids[1])], days;
                                         grid=(; pscale=ps, ΔT=[0.0]), subset=true, precipitation_gradient=g)

        # At g = 0 the read is the band sum node by node, exactly.
        full = GlacierHydrology.read_gemb_tile(path, area[geotile=At(ids[1])], days; grid=(; pscale=[1.0, 2.0, 4.0], ΔT=[0.0]))
        flat = tile_read(0.0, [1.0, 2.0, 4.0])
        @test all(flat.layers[k] == full.layers[k] for k in keys(full.layers))

        # Mean-centred on the band weights (10 and 30 km² at 1050 and 1150 m: z̄ = 1125 m), 5 % per 100 m:
        # m = 0.9625 and 1.0125, so at pscale 2 the bands sit at 1.925 and 2.025.
        w = [10.0, 30.0] .* 1e-3
        m = GlacierHydrology.gemb_band_pscale_multipliers([1050.0, 1150.0], w, 5.0)
        @test m ≈ [0.9625, 1.0125]
        @test sum(w .* m) ≈ sum(w)
        @test all(GlacierHydrology.gemb_band_pscale_multipliers([1050.0, 1150.0], w, 0.0) .== 1)
        g5 = tile_read(5.0, [2.0])
        n = length(times) - 1
        # Precipitation is linear in pscale, so the gradient moves it between bands and keeps the total.
        @test g5.layers[:acc][end] ≈ n * sum(w .* 100 .* 2 .* m) / 917 rtol = 1e-6
        @test g5.layers[:acc][end] ≈ flat.layers[:acc][end, 2] rtol = 1e-6
        # Runoff, 5 p² at the nodes, interpolated linearly between them: 5 (0.075·1 + 0.925·4) and
        # 5 (0.9875·4 + 0.0125·16).
        ro = [5 * (0.075 * 1 + 0.925 * 4), 5 * (0.9875 * 4 + 0.0125 * 16)]
        @test g5.layers[:melt][end] ≈ n * sum(w .* ro) / 917 rtol = 1e-6
        @test !(g5.layers[:melt][end] ≈ flat.layers[:melt][end, 2])

        # A band pushed off the tile's pscale axis is refused, with the range the tile supports.
        @test all(GlacierHydrology.gemb_gradient_pscale_range(m, w, [1.0, 2.0, 4.0]) .≈ (1 / 0.9625, 4 / 1.0125))
        @test_throws "outside the tile's axis 1.0..4.0" tile_read(5.0, [4.0])
        @test_throws "supports geotile-mean pscale 1.039..3.951 only" tile_read(5.0, [1.0])
        @test GlacierHydrology.gemb_gradient_pscale_range([-0.1, 1.5], w, [1.0, 4.0]) == (Inf, -Inf)

        # A gradient member covers only the geotiles that support every node of its grid.
        paths = GlacierHydrology.gemb_tile_paths(dir; covering=(Date(2000, 2, 15), Date(2000, 3, 15)))
        @test GlacierHydrology.gemb_gradient_supported_geotiles(ids, area, paths, 5.0, [2.0, 3.0]).supported == ids
        out = GlacierHydrology.gemb_gradient_supported_geotiles(ids, area, paths, 5.0, [1.0, 2.0])
        @test isempty(out.supported)
        @test out.excluded.id == ids
        @test out.excluded.pscale_lo[1] ≈ 1 / 0.9625 && out.excluded.pscale_hi[1] ≈ 4 / 1.0125
        @test GlacierHydrology.gemb_gradient_supported_geotiles(ids, area, paths, 2.5, [1.0, 2.0]).excluded.id == ids

        # The ensemble of a gradient member names it, so its fits cannot pass as the tile member's.
        gemb = GlacierHydrology.process_gemb_tiles(DataFrame(id=ids), area; tile_dir=dir, date_center=[Date(2000, 2, 15), Date(2000, 3, 15)],
                                      show_stats=false, physics_member="B-kfit-table-g5",
                                      forcing_grid=(; pscale=[2.0], ΔT=[0.0]))
        p = DimensionalData.metadata(gemb)["tile_provenance"]["tile_run_parameters"]
        @test p["physics_member"] == "B-kfit-table-g5"
        @test p["physics_precipitation_gradient"] == 5.0
        @test_throws "members must not be mixed" GlacierHydrology._check_gemb_physics_provenance(p, "B-kfit-table", "fit")
        @test_throws "is not a GEMB physics member id" GlacierHydrology.process_gemb_tiles(
            DataFrame(id=ids), area; tile_dir=dir, date_center=[Date(2000, 2, 15)], show_stats=false,
            physics_member="B-kfit-table-g0")

        # WGMS scoring reads each band at the same pscale, centred on the same hypsometry.
        fits = DataFrame(id=ids, pscale=[2.0], ΔT=[0.0])
        band = GlacierHydrology.wgms_tile_model(dir, fits; physics_member="B-kfit-table-g5", area_km2=area)(ids[1])
        # The band series accumulate from zero before the first output time, over all three intervals.
        @test band.smb[end, :] ≈ 3 .* (100 .* 2 .* m .- ro) ./ 1000 rtol = 1e-6
        plain = GlacierHydrology.wgms_tile_model(dir, fits; physics_member="B-kfit-table")(ids[1])
        @test plain.smb[end, :] ≈ 3 .* fill(200 - 20, 2) ./ 1000 rtol = 1e-6
        @test_throws "pass area_km2" GlacierHydrology.wgms_tile_model(dir, fits; physics_member="B-kfit-table-g5")

        # gemb_info: a gradient member reads its tile member's sweep, on the gradient grid, into a file
        # of its own; production has none yet.
        i5 = GlacierHydrology.gemb_info(; gemb_run_id=10, physics_member="B-kfit-table-g5")
        i0 = GlacierHydrology.gemb_info(; gemb_run_id=10, physics_member="B-kfit-table")
        @test i5.tile_dir == i0.tile_dir
        @test i5.filename_gemb_combined != i0.filename_gemb_combined
        @test i5.forcing_grid == GlacierHydrology.gemb_pilot_gradient_forcing_grid
        @test_throws "defined for the pilot sweep (GEMB run 10) only" GlacierHydrology.gemb_info(; gemb_run_id=9, physics_member="B-kfit-table-g5")
    end

    @testset "the pilot sweep: own tiles, complete groups, no pscale windows" begin
        times = [DateTime(2000, 1, 31, 23), DateTime(2000, 2, 29, 23), DateTime(2000, 3, 31, 23)]
        z = zeros(3, 1)
        flux = (; P=z, Rn=z, Ro=z, Rf=z, Mt=z, EC=z, firn=z)
        dir = mktempdir()
        for id in ("lat[+60+62]lon[-142-140]", "lat[+60+62]lon[-140-138]")
            write_test_tile(joinpath(dir, "tile_$(hash(id)).nc"), times, [1050.0]; flux..., geotile_id=id)
        end
        ids = ["lat[+58+60]lon[-138-136]", "lat[+60+62]lon[-140-138]", "lat[+60+62]lon[-142-140]"]
        coverage = GlacierHydrology.gemb_tile_coverage(dir)
        own(covering) = GlacierHydrology.gemb_tile_own_geotiles(ids, GlacierHydrology.gemb_tile_paths(coverage; covering))
        @test own((Date(2000, 2, 15), Date(2000, 3, 15))) == ids[2:3]
        @test isempty(own((Date(2000, 2, 15), Date(2001, 3, 15))))

        # A group is fitted only when all its geotiles are in the ensemble.
        groups = DataFrame(id=["a", "b", "c", "d"], group=[1, 1, 2, 3])
        present = DataFrame(id=["a", "c", "d"], group=[1, 2, 3], extra=[10, 20, 30])
        kept = GlacierHydrology.complete_geotile_groups(present, groups)
        @test kept.id == ["c", "d"]
        @test kept.extra == [20, 30]

        # The sanity windows can be skipped for the pilot only.
        run = joinpath(GlacierHydrology.pathlocal[:data_dir], "binned", "2deg",
                       "glacier_rgi7_dh_best_cc_nmad3_v01_filled_ac_p6_synthesized.jld2")
        dims4 = (Dim{:geotile}(["lat[+60+62]lon[-142-140]"]), Dim{:date}([Date(2000, 1, 15)]),
                 Dim{:pscale}([1.0]), Dim{:ΔT}([0.0]))
        dv(id) = DimArray(zeros(1, 1, 1, 1), dims4;
                          metadata=Dict("surface_mask" => "glacier_rgi7", "gemb_run_id" => id))
        @test_throws "pscale_sanity_check can be turned off only for a sweep over a subset of tiles (GEMB run 10)" GlacierHydrology.gemb_calibration(
            [run], dv(9); pscale_sanity_check=false)
        # For the pilot the guard passes, and the call fails later on the absent tile provenance.
        @test_throws "dv_gemb records no tile_provenance" GlacierHydrology.gemb_calibration([run], dv(10))
    end

    @testset "the calibrated synthesis reads one member's ensembles" begin
        run = joinpath(GlacierHydrology.pathlocal[:data_dir], "binned", "2deg",
                       "glacier_rgi7_dh_best_cc_nmad3_v01_filled_ac_p6_synthesized.jld2")
        dims4 = (Dim{:geotile}(["lat[+60+62]lon[-142-140]"]), Dim{:date}([Date(2000, 1, 15)]),
                 Dim{:pscale}([1.0]), Dim{:ΔT}([0.0]))
        stack(member) = DimStack((; dv=DimArray(zeros(1, 1, 1, 1), dims4;
                                                 metadata=Dict("surface_mask" => "glacier_rgi7",
                                                               "gemb_physics_member" => member))))
        discharge = DimArray(zeros(1), dims4[1]; metadata=Dict("surface_mask" => "glacier_rgi7"))
        @test GlacierHydrology.gemb_physics_member(stack("C-kfit-table")[:dv]) == "C-kfit-table"
        @test_throws "records no gemb_physics_member" GlacierHydrology.gemb_physics_member(DimArray(zeros(1), dims4[1]))
        run_rgi6 = replace(run, "glacier_rgi7" => "glacier")
        @test_throws "the GEMB ensembles passed are of different physics members (\"B-kfit-table\", \"C-kfit-table\")" GlacierHydrology.geotile_synthesis_gembfit_dv(
            [run, run_rgi6],
            Dict("glacier_rgi7" => discharge,
                 "glacier" => DimArray(zeros(1), dims4[1]; metadata=Dict("surface_mask" => "glacier"))),
            Dict("glacier_rgi7" => stack("B-kfit-table"),
                 "glacier" => DimStack((; dv=DimArray(zeros(1, 1, 1, 1), dims4;
                                                       metadata=Dict("surface_mask" => "glacier",
                                                                     "gemb_physics_member" => "C-kfit-table")))));
            geotile_width=2, force_remake_before=nothing)
    end
end

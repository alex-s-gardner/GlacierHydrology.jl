using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
import GlacierVolumeChange as GVC
using DimensionalData
using Dates
using DataFrames
using Statistics

@testset "Project configuration" begin
    @testset "gemb_info" begin
        # Every point-run entry must be constructible: the pscale and Δheight grids carry one
        # label per value, and filename_gemb_combined resolves its interpolated id.
        for gemb_run_id in 1:6
            g = GGA.gemb_info(; gemb_run_id)

            @test length(g.precipitation_scale) == length(dims(g.precipitation_scale, :pscale))
            @test length(g.elevation_delta) == length(dims(g.elevation_delta, :Δheight))

            @test !occursin("\$", g.filename_gemb_combined)
            @test occursin(g.file_uniqueid, g.filename_gemb_combined)
            @test endswith(g.filename_gemb_combined, ".jld2")

            @test g.modify_melt_only isa Bool
        end

        # The tile path has a different shape: its forcing axes live in the data, so it declares a
        # directory to read instead of a precipitation grid and an elevation-class list. Runs 9 and 10
        # hold one sweep per physics member. The reference carries a precipitation gradient, which only
        # the pilot sweep (run 10) defines, so run 9 is read at the reference's tile member.
        ref = GGA.gemb_reference_physics_member
        ref_tiles = GGA.gemb_physics_member_parts(ref).tile_member
        @test GGA.gemb_physics_member_parts(ref).precipitation_gradient > 0
        tile_info(gemb_run_id) = gemb_run_id == 9 ? GGA.gemb_info(; gemb_run_id, physics_member=ref_tiles) :
                                 gemb_run_id == 10 ? GGA.gemb_info(; gemb_run_id, physics_member=ref) :
                                 GGA.gemb_info(; gemb_run_id)
        for gemb_run_id in (7, 8, 9, 10)
            tiles = tile_info(gemb_run_id)
            @test tiles.tile_dir isa AbstractString
            @test isnothing(tiles.precipitation_scale)
            @test !occursin("\$", tiles.filename_gemb_combined)
            @test endswith(tiles.filename_gemb_combined, ".jld2")
            @test tiles.modify_melt_only isa Bool
        end

        # Each tile run reads its own directory and writes its own combined file.
        tile_runs = [tile_info(gemb_run_id) for gemb_run_id in (7, 8, 9, 10)]
        @test allunique(getproperty.(tile_runs, :tile_dir))
        @test allunique(getproperty.(tile_runs, :filename_gemb_combined))
        @test splitpath(tile_info(9).tile_dir)[end-1:end] == ["tile_runs_v2", ref_tiles]
        @test splitpath(tile_info(10).tile_dir)[end-1:end] == ["tile_runs_pilot", ref_tiles]
        @test tile_info(9).physics_member == ref_tiles
        @test tile_info(10).physics_member == ref
        @test occursin(ref, tile_info(10).filename_gemb_combined)
        @test_throws "defined for the pilot sweep (GEMB run 10) only" GGA.gemb_info(; gemb_run_id=9, physics_member=ref)
        # The pilot is read on its own grid, which every pilot tile carries, and a gradient member on the
        # geotile-mean nodes its bands stay on; production on the tiles' own.
        @test GGA.gemb_info(; gemb_run_id=10, physics_member=ref_tiles).forcing_grid == GGA.gemb_pilot_forcing_grid
        @test tile_info(10).forcing_grid == GGA.gemb_pilot_gradient_forcing_grid
        @test isnothing(tile_info(9).forcing_grid)
        @test 1.0 in GGA.gemb_pilot_forcing_grid.pscale && 0.0 in GGA.gemb_pilot_forcing_grid.ΔT

        # Each physics member of a run has its own directory and combined file, and a member is
        # required where there are members and refused where there are none.
        members = ["B-kfit-table", "C-kfit-table", "C-kshaw-table"]
        infos = [GGA.gemb_info(; gemb_run_id=9, physics_member=m) for m in members]
        @test allunique(getproperty.(infos, :tile_dir))
        @test allunique(getproperty.(infos, :filename_gemb_combined))
        @test all(occursin(m, i.filename_gemb_combined) for (m, i) in zip(members, infos))
        @test_throws "GEMB run 9 holds one tile sweep per physics member; pass physics_member" GGA.gemb_info(; gemb_run_id=9)
        @test_throws "GEMB run 8 has no physics members" GGA.gemb_info(; gemb_run_id=8, physics_member=ref)
        @test_throws "is not a GEMB physics member id" GGA.gemb_info(; gemb_run_id=9, physics_member="C_ks35_table")
        @test_throws "is not a GEMB physics member id" GGA.gemb_info(; gemb_run_id=9, physics_member="../x-y-z")

        @test_throws "unrecognized gemb_run_id" GGA.gemb_info(; gemb_run_id=11)
    end

    @testset "file_is_current" begin
        absent = joinpath(mktempdir(), "not_written.jld2")

        # A missing file is never current, and must not compare an mtime against `nothing`.
        @test GGA.file_is_current(absent, nothing) == false
        @test GGA.file_is_current(absent, DateTime(2020, 1, 1)) == false

        present = tempname()
        write(present, "x")
        try
            # No rebuild date means any existing file is reusable.
            @test GGA.file_is_current(present, nothing) == true

            # Written after the rebuild date -> reusable; written before -> rebuild.
            @test GGA.file_is_current(present, DateTime(2000, 1, 1)) == true
            @test GGA.file_is_current(present, DateTime(2999, 1, 1)) == false
        finally
            rm(present; force=true)
        end
    end

    @testset "project_paths covers every product" begin
        products = GVC.project_products(project_id=:v01)
        paths = GVC.project_paths(project_id=:v01)
        @test keys(paths) == keys(products)
        for mission in keys(products)
            @test occursin(string(products[mission].name), paths[mission].geotile)
        end
    end

    @testset "GEMB calibration cost penalty" begin
        res = DimArray(sin.(2π .* (1:120) ./ 12) .+ 0.1 .* (1:120) ./ 120, Dim{:date}(DateTime(2010, 1, 15):Month(1):DateTime(2019, 12, 15)))
        kw = (seasonality_weight=0.85, distance_from_origin_penalty=GGA.distance_from_origin_penalty,
            ΔT_to_pscale_weight=GGA.ΔT_to_pscale_weight)
        prior = GGA.gemb_forcing_prior

        # both penalties are a factor ≥ 1 that is exactly 1 at their centre
        bare = GGA.model_fit_cost_function(res, 1.0, 0.0; kw..., origin_penalty_mode=:legacy)
        @test GGA.model_fit_cost_function(res, prior.pscale, prior.ΔT; kw...) ≈ bare
        @test GGA.model_fit_cost_function(res, 1.0, 0.0; kw...) > bare
        @test GGA.model_fit_cost_function(res, prior.pscale, prior.ΔT; kw..., origin_penalty_mode=:legacy) > bare

        # the prior distance is a Mahalanobis distance in (log pscale, ΔT)
        @test GGA._forcing_prior_distance(prior.pscale * exp(prior.log_pscale_sd), prior.ΔT, (; prior..., corr=0.0)) ≈ 1
        @test GGA._forcing_prior_distance(prior.pscale, prior.ΔT + prior.ΔT_sd, (; prior..., corr=0.0)) ≈ 1

        @test_throws "origin_penalty_mode must be :legacy or :prior" GGA.model_fit_cost_function(res, 1.0, 0.0; kw..., origin_penalty_mode=:additive)
        @test_throws "assumes a multiplicative pscale and an additive ΔT" GGA.model_fit_cost_function(res, 1.0, 0.0; kw...,
            is_scaling_factor=Dict("pscale" => true, "ΔT" => true))
    end

    @testset "GEMB calibration grid search" begin
        pscale_grid = exp.(range(log(0.25), log(4), 41))
        ΔT_grid = collect(-3:0.25:6)
        f(x) = (log(x[1]) - log(1.7))^2 + 0.1 * (x[2] - 2.3)^2
        pscale, ΔT, cost = GGA._gemb_grid_minimize(f, pscale_grid, ΔT_grid)
        @test pscale ≈ 1.7 rtol = 0.01
        @test ΔT ≈ 2.3 atol = 0.03
        @test cost < 1e-4
    end
end

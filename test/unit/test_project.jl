using Test
using GlacierHydrology
import GlacierHydrology
import GlacierVolumeChange as GVC
using DimensionalData
using Dates
using DataFrames
using Statistics

@testset "Project configuration" begin
    @testset "gemb_info" begin
        # Runs 9 and 10 hold one sweep per physics member. The reference member carries no
        # precipitation gradient; a gradient member (only the pilot sweep, run 10, defines them) reads
        # its tile member's directory under a combined file of its own.
        ref = GlacierHydrology.gemb_reference_physics_member
        @test GlacierHydrology.gemb_physics_member_parts(ref) == (; tile_member=ref, precipitation_gradient=0.0)
        grad = GlacierHydrology.gemb_physics_member_id(ref, 2.5)
        for gemb_run_id in (9, 10)
            tiles = GlacierHydrology.gemb_info(; gemb_run_id, physics_member=ref)
            @test tiles.tile_dir isa AbstractString
            @test !occursin("\$", tiles.filename_gemb_combined)
            @test endswith(tiles.filename_gemb_combined, ".jld2")
            @test tiles.physics_member == ref
        end

        # Each tile run reads its own directory and writes its own combined file.
        run9 = GlacierHydrology.gemb_info(; gemb_run_id=9, physics_member=ref)
        run10 = GlacierHydrology.gemb_info(; gemb_run_id=10, physics_member=ref)
        run10_grad = GlacierHydrology.gemb_info(; gemb_run_id=10, physics_member=grad)
        @test allunique(getproperty.((run9, run10, run10_grad), :filename_gemb_combined))
        @test run9.tile_dir != run10.tile_dir
        @test splitpath(run9.tile_dir)[end-1:end] == ["tile_runs_v2", ref]
        @test splitpath(run10.tile_dir)[end-1:end] == ["tile_runs_pilot", ref]
        @test run10_grad.tile_dir == run10.tile_dir
        @test run10_grad.physics_member == grad
        @test occursin(grad, run10_grad.filename_gemb_combined)
        @test_throws "defined for the pilot sweep (GEMB run 10) only" GlacierHydrology.gemb_info(; gemb_run_id=9, physics_member=grad)
        # The pilot is read on its own grid, which every pilot tile carries, and a gradient member on the
        # geotile-mean nodes its bands stay on; production on the tiles' own.
        @test run10.forcing_grid == GlacierHydrology.gemb_pilot_forcing_grid
        @test run10_grad.forcing_grid == GlacierHydrology.gemb_pilot_gradient_forcing_grid
        @test isnothing(run9.forcing_grid)
        @test 1.0 in GlacierHydrology.gemb_pilot_forcing_grid.pscale && 0.0 in GlacierHydrology.gemb_pilot_forcing_grid.ΔT

        # Each physics member of a run has its own directory and combined file, and a member is required.
        members = ["B-kfit-table", "C-kfit-table", "C-kshaw-table"]
        infos = [GlacierHydrology.gemb_info(; gemb_run_id=9, physics_member=m) for m in members]
        @test allunique(getproperty.(infos, :tile_dir))
        @test allunique(getproperty.(infos, :filename_gemb_combined))
        @test all(occursin(m, i.filename_gemb_combined) for (m, i) in zip(members, infos))
        @test_throws UndefKeywordError GlacierHydrology.gemb_info(; gemb_run_id=9)
        @test_throws "is not a GEMB physics member id" GlacierHydrology.gemb_info(; gemb_run_id=9, physics_member="C_ks35_table")
        @test_throws "is not a GEMB physics member id" GlacierHydrology.gemb_info(; gemb_run_id=9, physics_member="../x-y-z")

        @test_throws "unrecognized gemb_run_id" GlacierHydrology.gemb_info(; gemb_run_id=8, physics_member=ref)
        @test_throws "unrecognized gemb_run_id" GlacierHydrology.gemb_info(; gemb_run_id=11, physics_member=ref)

        # The run's settings: the pilot covers a subset of tiles, and neither run has a forcing prior yet.
        @test run9.gemb_run_id == 9 && !run9.partial_coverage
        @test run10.partial_coverage
        @test all(isnothing(i.forcing_prior) && iszero(i.distance_from_origin_penalty) for i in (run9, run10))
    end

    @testset "file_is_current" begin
        absent = joinpath(mktempdir(), "not_written.jld2")

        # A missing file is never current, and must not compare an mtime against `nothing`.
        @test GlacierHydrology.file_is_current(absent, nothing) == false
        @test GlacierHydrology.file_is_current(absent, DateTime(2020, 1, 1)) == false

        present = tempname()
        write(present, "x")
        try
            # No rebuild date means any existing file is reusable.
            @test GlacierHydrology.file_is_current(present, nothing) == true

            # Written after the rebuild date -> reusable; written before -> rebuild.
            @test GlacierHydrology.file_is_current(present, DateTime(2000, 1, 1)) == true
            @test GlacierHydrology.file_is_current(present, DateTime(2999, 1, 1)) == false
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
        kw = (seasonality_weight=0.85, distance_from_origin_penalty=0.35,
            ΔT_to_pscale_weight=GlacierHydrology.ΔT_to_pscale_weight)
        prior = (pscale=1.60, log_pscale_sd=0.45, ΔT=2.20, ΔT_sd=1.65, corr=0.50)

        # both penalties are a factor ≥ 1 that is exactly 1 at their centre
        bare = GlacierHydrology.model_fit_cost_function(res, 1.0, 0.0; kw..., origin_penalty_mode=:legacy)
        @test GlacierHydrology.model_fit_cost_function(res, prior.pscale, prior.ΔT; kw..., forcing_prior=prior) ≈ bare
        @test GlacierHydrology.model_fit_cost_function(res, 1.0, 0.0; kw..., forcing_prior=prior) > bare
        @test GlacierHydrology.model_fit_cost_function(res, prior.pscale, prior.ΔT; kw..., origin_penalty_mode=:legacy) > bare

        # the prior distance is a Mahalanobis distance in (log pscale, ΔT)
        @test GlacierHydrology._forcing_prior_distance(prior.pscale * exp(prior.log_pscale_sd), prior.ΔT, (; prior..., corr=0.0)) ≈ 1
        @test GlacierHydrology._forcing_prior_distance(prior.pscale, prior.ΔT + prior.ΔT_sd, (; prior..., corr=0.0)) ≈ 1

        # without a prior the penalty must be zero, and is then no penalty
        @test GlacierHydrology.model_fit_cost_function(res, 3.0, 4.0; kw..., distance_from_origin_penalty=0) ≈ bare
        @test_throws "a forcing-prior penalty needs a forcing_prior" GlacierHydrology.model_fit_cost_function(res, 1.0, 0.0; kw...)

        @test_throws "origin_penalty_mode must be :legacy or :prior" GlacierHydrology.model_fit_cost_function(res, 1.0, 0.0; kw..., origin_penalty_mode=:additive)
        @test_throws "assumes a multiplicative pscale and an additive ΔT" GlacierHydrology.model_fit_cost_function(res, 1.0, 0.0; kw..., forcing_prior=prior,
            is_scaling_factor=Dict("pscale" => true, "ΔT" => true))
    end

    @testset "GEMB calibration grid search" begin
        pscale_grid = exp.(range(log(0.25), log(4), 41))
        ΔT_grid = collect(-3:0.25:6)
        f(x) = (log(x[1]) - log(1.7))^2 + 0.1 * (x[2] - 2.3)^2
        pscale, ΔT, cost = GlacierHydrology._gemb_grid_minimize(f, pscale_grid, ΔT_grid)
        @test pscale ≈ 1.7 rtol = 0.01
        @test ΔT ≈ 2.3 atol = 0.03
        @test cost < 1e-4
    end
end

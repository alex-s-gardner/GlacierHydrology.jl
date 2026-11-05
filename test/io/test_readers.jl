"""
Tests for external data readers in utilities_readers.jl

These tests validate parsing and unit conversions for external datasets including
glacier discharge, GRACE satellite gravity, and GlaMBIE mass balance data.
"""

using Test
using GlacierHydrology
import GlacierHydrology
using CSV
using DataFrames
using Statistics
using MAT
using Dates

# Include test fixtures
include("../fixtures/mock_data_files.jl")

@testset "External Data Readers" begin

    @testset "Glacier discharge reading" begin
        mktempdir() do temp_dir
            # Create mock discharge CSV file
            discharge_file = joinpath(temp_dir, "discharge_nh.csv")
            create_mock_discharge_csv(discharge_file; n_glaciers=8)

            # Read the CSV manually to test format.
            # The fixture writes 14 metadata lines, so the column header is line 15 and data starts
            # on line 16. (`header=14` pointed at the last "#" comment, so every column name came
            # back wrong and the row count was off.)
            df = CSV.read(discharge_file, DataFrame; header=15, skipto=16)

            # Test: Should have required columns
            @test hasproperty(df, :RGIId)
            @test hasproperty(df, :Year)
            @test hasproperty(df, :Discharge_Gt_yr)
            @test hasproperty(df, :Uncertainty_Gt_yr)

            # Test: Data dimensions
            n_years = 21  # 2000-2020
            n_glaciers = 8
            @test nrow(df) == n_years * n_glaciers

            # Test: RGI IDs should be properly formatted
            @test all(startswith.(df.RGIId, "RGI60-"))

            # Test: Years should be in expected range
            @test minimum(df.Year) >= 2000
            @test maximum(df.Year) <= 2020

            # Test: Discharge values should be positive
            @test all(df.Discharge_Gt_yr .>= 0)
            @test all(df.Uncertainty_Gt_yr .>= 0)

            # Test: Uncertainties should be smaller than discharge values (reasonable)
            @test mean(df.Uncertainty_Gt_yr) < mean(df.Discharge_Gt_yr)
        end
    end

    @testset "GRACE data reading" begin
        mktempdir() do temp_dir
            # Create mock GRACE .mat file
            grace_file = joinpath(temp_dir, "grace_rgi.mat")
            create_mock_grace_mat(grace_file; n_regions=5, n_times=240)

            # Read the .mat file
            grace = matread(grace_file)

            # Test: Should have required variables
            @test haskey(grace, "region_codes")
            @test haskey(grace, "time")
            @test haskey(grace, "mass_change_Gt")
            @test haskey(grace, "uncertainty_Gt")

            # Test: Dimensions should match
            n_regions = length(grace["region_codes"])
            n_times = length(grace["time"])
            @test n_regions == 5
            @test n_times == 240

            @test size(grace["mass_change_Gt"]) == (n_regions, n_times)
            @test size(grace["uncertainty_Gt"]) == (n_regions, n_times)

            # Test: Region codes should be strings.
            # NOTE: this previously read `grace["region_codes"] .isa String for _ in ...`, which is
            # not valid Julia -- the whole file failed to parse, so none of its tests ever ran.
            @test all(code isa String for code in grace["region_codes"])

            # Test: Time should be decimal years
            @test all(grace["time"] .>= 2003.0)
            @test all(grace["time"] .<= 2025.0)

            # Test: Time should be monotonically increasing
            @test all(diff(grace["time"]) .> 0)

            # Test: Mass change should have reasonable values (negative for most regions)
            # Most glacier regions are losing mass
            final_mass_changes = grace["mass_change_Gt"][:, end]
            @test sum(final_mass_changes .< 0) >= 3  # At least 3/5 regions losing mass

            # Test: Uncertainties should be positive
            @test all(grace["uncertainty_Gt"] .>= 0)

            # Test: Uncertainties should be reasonable (not larger than absolute mass changes)
            for i in 1:n_regions
                mean_uncertainty = mean(grace["uncertainty_Gt"][i, :])
                mean_abs_change = mean(abs.(grace["mass_change_Gt"][i, :]))
                @test mean_uncertainty < mean_abs_change * 0.5  # Uncertainty < 50% of signal
            end
        end
    end

    @testset "GlaMBIE 2024 reading" begin
        mktempdir() do temp_dir
            # Create mock GlaMBIE CSV
            glambie_file = joinpath(temp_dir, "glambie_2024.csv")
            create_mock_glambie_csv(glambie_file; n_years=25)

            # Read the CSV
            df = CSV.read(glambie_file, DataFrame)

            # Test: Should have required columns
            @test hasproperty(df, :Year)
            @test hasproperty(df, :RGI_Region)
            @test hasproperty(df, :Mass_Balance_Gt)
            @test hasproperty(df, :Uncertainty_Gt)

            # Test: Years should span 2000-2024
            years_present = unique(df.Year)
            @test minimum(years_present) == 2000
            @test maximum(years_present) == 2024

            # Test: RGI regions should include 1-19, 98, 99
            regions_present = unique(df.RGI_Region)
            @test 1 in regions_present
            @test 19 in regions_present
            @test 98 in regions_present  # Global
            @test 99 in regions_present  # Global excluding Greenland/Antarctica

            # Test: Each region should have data for all years
            n_years = length(years_present)
            for region in regions_present
                region_data = filter(row -> row.RGI_Region == region, df)
                @test nrow(region_data) == n_years
            end

            # Test: Mass balance should be cumulative (generally decreasing)
            # Check region 98 (global)
            global_data = filter(row -> row.RGI_Region == 98, df)
            sort!(global_data, :Year)

            # Cumulative mass balance should trend negative
            @test global_data.Mass_Balance_Gt[end] < global_data.Mass_Balance_Gt[1]

            # Test: Uncertainties should be positive and growing with time
            @test all(global_data.Uncertainty_Gt .>= 0)
            @test global_data.Uncertainty_Gt[end] > global_data.Uncertainty_Gt[1]

            # Test: Global (98) should be roughly sum of regions
            # (This is a rough check since we have synthetic data)
            year_2020 = filter(row -> row.Year == 2020, df)
            global_2020 = filter(row -> row.RGI_Region == 98, year_2020)[1, :Mass_Balance_Gt]
            regions_sum_2020 = sum(filter(row -> 1 <= row.RGI_Region <= 19, year_2020).Mass_Balance_Gt)

            # Global should be roughly in the same order of magnitude as sum of regions
            @test abs(global_2020) > abs(regions_sum_2020) * 0.5
            @test abs(global_2020) < abs(regions_sum_2020) * 2.0
        end
    end

    @testset "Unit conversions - Gt to mm SLE" begin
        # Test standard unit conversions used in readers

        # 1 Gt = 1e12 kg = 1e9 m³ of water (density = 1000 kg/m³)
        mass_Gt = 1.0
        volume_m3 = mass_Gt * 1e12 / 1000  # kg to m³
        @test volume_m3 ≈ 1e9

        # Spread over ocean area (362.5 million km²)
        ocean_area_m2 = 362.5e12  # Convert to m²
        sea_level_m = volume_m3 / ocean_area_m2
        sea_level_mm = sea_level_m * 1000

        # 1 Gt should equal ~0.00276 mm SLE
        @test sea_level_mm ≈ 0.00276 rtol=0.01

        # Test: 360 Gt loss per year = ~1 mm/yr SLE
        annual_loss_Gt = 360.0
        sle_mm_yr = annual_loss_Gt * 0.00276
        @test sle_mm_yr ≈ 1.0 rtol=0.05
    end

    @testset "Data consistency checks" begin
        mktempdir() do temp_dir
            # Create mock files
            paths = create_mock_data_directory(temp_dir)

            # Test: All files should exist and be readable
            @test isfile(paths[:discharge_nh])
            @test isfile(paths[:grace_rgi])
            @test isfile(paths[:glambie_2024])

            # Test: Files should have content
            @test filesize(paths[:discharge_nh]) > 100  # Has header + data
            @test filesize(paths[:grace_rgi]) > 100
            @test filesize(paths[:glambie_2024]) > 100

            # Test: Files should be parseable
            @test_nowarn CSV.read(paths[:discharge_nh], DataFrame; header=15, skipto=16)
            @test_nowarn matread(paths[:grace_rgi])
            @test_nowarn CSV.read(paths[:glambie_2024], DataFrame)
        end
    end

    @testset "Edge case - empty regions" begin
        # Test handling of regions with no data
        mktempdir() do temp_dir
            # Create GlaMBIE with missing data for some regions
            glambie_file = joinpath(temp_dir, "glambie_sparse.csv")

            # Create data with gaps
            records = []
            for region in [1, 3, 5]  # Only some regions
                for year in 2000:2010
                    push!(records, (
                        Year = year,
                        RGI_Region = region,
                        Mass_Balance_Gt = -10.0 * region * (year - 2000),
                        Uncertainty_Gt = 2.0
                    ))
                end
            end

            df = DataFrame(records)
            CSV.write(glambie_file, df)

            # Read and check
            df_read = CSV.read(glambie_file, DataFrame)
            @test nrow(df_read) == 3 * 11  # 3 regions × 11 years

            # Test: Only specified regions should be present
            regions = unique(df_read.RGI_Region)
            @test regions == [1, 3, 5]
            @test !(2 in regions)
            @test !(4 in regions)
        end
    end

    @testset "Temporal coverage checks" begin
        mktempdir() do temp_dir
            # Test different temporal ranges for different datasets

            # GRACE: ~2003-2023
            grace_file = joinpath(temp_dir, "grace_short.mat")
            create_mock_grace_mat(grace_file; n_regions=3, n_times=120)  # 10 years
            grace = matread(grace_file)

            @test length(grace["time"]) == 120
            @test maximum(grace["time"]) - minimum(grace["time"]) ≈ 120/12 rtol=0.1  # ~10 years

            # GlaMBIE: 2000-2024
            glambie_file = joinpath(temp_dir, "glambie_full.csv")
            create_mock_glambie_csv(glambie_file; n_years=25)
            df_glambie = CSV.read(glambie_file, DataFrame)

            years_glambie = unique(df_glambie.Year)
            @test length(years_glambie) == 25
            @test minimum(years_glambie) == 2000
            @test maximum(years_glambie) == 2024
        end
    end

    @testset "RGI ID format validation" begin
        mktempdir() do temp_dir
            # Test that RGI IDs follow expected format
            discharge_file = joinpath(temp_dir, "discharge_format_test.csv")
            create_mock_discharge_csv(discharge_file; n_glaciers=5)

            df = CSV.read(discharge_file, DataFrame; header=15, skipto=16)

            # Test: RGI IDs should match pattern RGI60-XX.XXXXX
            for rgi_id in df.RGIId
                @test occursin(r"^RGI60-\d{2}\.\d{5,}$", rgi_id)
            end

            # Test: Region numbers should be valid (01-19)
            for rgi_id in df.RGIId
                region_str = split(rgi_id, '-')[2][1:2]
                region_num = parse(Int, region_str)
                @test 1 <= region_num <= 19
            end
        end
    end

    @testset "Uncertainty propagation" begin
        mktempdir() do temp_dir
            # Test that uncertainties are properly formatted

            # Discharge uncertainties
            discharge_file = joinpath(temp_dir, "discharge_uncert.csv")
            create_mock_discharge_csv(discharge_file; n_glaciers=3)
            df_discharge = CSV.read(discharge_file, DataFrame; header=15, skipto=16)

            # Test: Relative uncertainty should be reasonable (5-30%)
            rel_uncertainty = df_discharge.Uncertainty_Gt_yr ./ df_discharge.Discharge_Gt_yr
            @test all(0.01 .< rel_uncertainty .< 0.50)  # 1-50% uncertainty

            # GRACE uncertainties
            grace_file = joinpath(temp_dir, "grace_uncert.mat")
            create_mock_grace_mat(grace_file; n_regions=2, n_times=60)
            grace = matread(grace_file)

            # Test: Uncertainties should be positive
            @test all(grace["uncertainty_Gt"] .>= 0)

            # Test: Uncertainties should be smaller than signal amplitude
            for i in 1:size(grace["mass_change_Gt"], 1)
                signal_range = maximum(grace["mass_change_Gt"][i, :]) - minimum(grace["mass_change_Gt"][i, :])
                mean_uncert = mean(grace["uncertainty_Gt"][i, :])
                @test mean_uncert < signal_range
            end
        end
    end

end

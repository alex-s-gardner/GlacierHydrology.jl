"""
Tests for regional aggregation and postprocessing in utilities_postprocessing.jl

These tests validate aggregation from geotile to RGI regions, trend fitting,
ensemble statistics, and regional summary generation.
"""

using Test
using GlacierHydrology
import GlacierHydrology
import GlacierVolumeChange as GVC
using DimensionalData
import DimensionalData as DD
using Dates
using Statistics
using Random

# Include test fixtures
include("../fixtures/synthetic_timeseries.jl")

@testset "Postprocessing and Aggregation" begin

    @testset "Geotile to RGI aggregation - simple" begin
        # Create synthetic geotile time series with known trends
        n_dates = 24
        n_heights = 5
        n_geotiles = 3

        dates = [DateTime(2018,1,1) + Month(i) for i in 0:n_dates-1]
        heights = collect(range(1500, 2500, length=n_heights))
        geotiles = ["lat+60+62lon-050-048",
                   "lat+62+64lon-050-048",
                   "lat+64+66lon-050-048"]

        # Create elevation change data: each geotile has different trend
        dh_data = zeros(n_geotiles, n_dates, n_heights)
        area_data = zeros(n_geotiles, n_heights)

        trends = [-0.3, -0.5, -0.4]  # m/yr per geotile
        for i in 1:n_geotiles
            t = collect(0:n_dates-1) ./ 12.0  # Years
            for j in 1:n_dates
                for k in 1:n_heights
                    dh_data[i, j, k] = trends[i] * t[j] + 0.05 * randn()
                end
            end
            # Equal area per height bin
            area_data[i, :] .= 10.0  # km² per bin
        end

        # Create DimArrays
        geotile_dim = DD.Dim{:geotile}(geotiles)
        date_dim = DD.Dim{:date}(dates)
        height_dim = DD.Dim{:height}(heights)

        dh = DimArray(dh_data, (geotile_dim, date_dim, height_dim))
        area = DimArray(area_data, (geotile_dim, height_dim))

        # Compute volume change per geotile
        dv = GVC.dh2dv_geotile(dh, area)

        # Aggregate across geotiles (sum)
        dv_regional = sum(dv; dims=:geotile)[1, :]

        # Test: Regional volume should be sum of geotile volumes
        for j in 1:n_dates
            expected = sum(dv[i, j] for i in 1:n_geotiles)
            @test dv_regional[j] ≈ expected rtol=1e-10
        end

        # Test: Regional trend should be area-weighted average of geotile trends
        total_area = sum(area_data)
        expected_regional_trend = sum(trends[i] * sum(area_data[i, :]) for i in 1:n_geotiles) / total_area

        # Fit trend to regional data
        t = collect(0:n_dates-1) ./ 12.0
        # Simple linear regression: dh/dt = (Σt·dh - n·mean(t)·mean(dh)) / (Σt² - n·mean(t)²)
        # For volume: dV/dt ≈ trend * total_area / 1000
        expected_dv_trend = expected_regional_trend * total_area / 1000  # km³/yr

        # Actual trend from data
        actual_dv_trend = (dv_regional[end] - dv_regional[1]) / (t[end] - t[1])

        @test actual_dv_trend ≈ expected_dv_trend rtol=0.15  # Allow noise
    end

    @testset "Regional trend fitting - linear" begin
        # Test trend fitting to time series
        n_years = 10
        dates, values = generate_trend_seasonal_ts(;
            n_points=n_years*4,  # Quarterly
            trend=-0.5,  # m/yr
            amplitude=0.2,
            noise_sigma=0.05,
            start_year=2015
        )

        # Fit linear trend (simple least squares)
        t = collect(0:length(values)-1) ./ 4.0  # Years from start

        # y = a + b*t
        n = length(t)
        sum_t = sum(t)
        sum_y = sum(values)
        sum_tt = sum(t .^ 2)
        sum_ty = sum(t .* values)

        slope = (n * sum_ty - sum_t * sum_y) / (n * sum_tt - sum_t^2)
        intercept = (sum_y - slope * sum_t) / n

        # Test: Recovered slope should match input trend
        @test slope ≈ -0.5 rtol=0.15  # Within 15% given noise and seasonality

        # Test: Detrended data should have reduced variance
        fitted = intercept .+ slope .* t
        residuals = values .- fitted
        @test std(residuals) < std(values)
    end

    @testset "Ensemble quantile calculation" begin
        # Create synthetic ensemble with known distribution
        Random.seed!(123)
        n_ensemble = 100
        n_times = 20

        # Each ensemble member: trend + random offset
        ensemble_data = zeros(n_ensemble, n_times)
        true_trend = -0.4

        for i in 1:n_ensemble
            offset = 0.2 * randn()
            t = collect(0:n_times-1) ./ 12.0
            ensemble_data[i, :] = true_trend .* t .+ offset .+ 0.05 .* randn(n_times)
        end

        # Compute quantiles
        p05 = [quantile(ensemble_data[:, j], 0.05) for j in 1:n_times]
        p50 = [quantile(ensemble_data[:, j], 0.50) for j in 1:n_times]
        p95 = [quantile(ensemble_data[:, j], 0.95) for j in 1:n_times]

        # Test: Median should be close to true trend line
        t = collect(0:n_times-1) ./ 12.0
        expected_median = true_trend .* t
        @test mean(abs.(p50 .- expected_median)) < 0.15

        # Test: 90% of data should fall within p05-p95
        for j in 1:n_times
            in_range = sum((ensemble_data[:, j] .>= p05[j]) .&&
                          (ensemble_data[:, j] .<= p95[j]))
            @test in_range / n_ensemble > 0.85
            @test in_range / n_ensemble < 0.95
        end

        # Test: spread is stable in time. Each member here is `trend*t + fixed offset + noise`,
        # with a per-member offset drawn once and homoscedastic noise, so the ensemble spread does
        # not grow with time -- asserting that it does was testing a property the construction
        # deliberately lacks. A random-walk error model would be needed for growing spread.
        spread_start = p95[1] - p05[1]
        spread_end = p95[end] - p05[end]
        @test spread_end ≈ spread_start rtol=0.35
        @test spread_start > 0 && spread_end > 0
    end

    @testset "Error propagation in aggregation" begin
        # Test that errors combine correctly when aggregating regions
        n_regions = 5
        values = [10.0, 15.0, 8.0, 12.0, 20.0]  # Gt/yr per region
        errors = [1.0, 1.5, 0.8, 1.2, 2.0]      # Uncertainty per region

        # Aggregate: total = sum(values)
        total = sum(values)
        @test total ≈ 65.0

        # Error propagation for independent measurements:
        # σ_total = sqrt(Σ σ_i²)
        error_total = sqrt(sum(errors .^ 2))
        expected_error = sqrt(1.0^2 + 1.5^2 + 0.8^2 + 1.2^2 + 2.0^2)
        @test error_total ≈ expected_error rtol=1e-10

        # Test: Relative uncertainty decreases with aggregation
        rel_errors = errors ./ values
        rel_error_total = error_total / total
        @test rel_error_total < mean(rel_errors)
    end

    @testset "Area-weighted averaging" begin
        # Test area-weighted mean calculation
        values = [1.0, 2.0, 3.0, 4.0]  # m/yr at different locations
        areas = [10.0, 20.0, 15.0, 5.0]  # km² at each location

        # Area-weighted mean: Σ(value_i * area_i) / Σ(area_i)
        weighted_mean = sum(values .* areas) / sum(areas)
        expected = (1.0*10 + 2.0*20 + 3.0*15 + 4.0*5) / (10 + 20 + 15 + 5)
        @test weighted_mean ≈ expected rtol=1e-10
        @test weighted_mean ≈ 2.3 rtol=1e-6

        # Test: Weighted mean should differ from simple mean
        simple_mean = mean(values)
        @test weighted_mean != simple_mean

        # Test: If all areas equal, weighted mean = simple mean
        equal_areas = fill(10.0, 4)
        weighted_equal = sum(values .* equal_areas) / sum(equal_areas)
        @test weighted_equal ≈ simple_mean rtol=1e-10
    end

    @testset "Temporal centering to reference period" begin
        # Test centering time series to a reference period mean
        n_points = 60
        dates, values = generate_trend_seasonal_ts(;
            n_points=n_points,
            trend=0.5,
            amplitude=1.0,
            noise_sigma=0.1,
            start_year=2000
        )

        # Reference period: 2005-2010 (indices 20-40)
        ref_indices = 21:40  # Year 5-10
        ref_mean = mean(values[ref_indices])

        # Center to reference period
        centered_values = values .- ref_mean

        # Test: Reference period mean should be zero
        @test mean(centered_values[ref_indices]) ≈ 0.0 atol=1e-10

        # Test: Shape preserved (only offset changed)
        @test std(centered_values) ≈ std(values) rtol=1e-10
        @test diff(centered_values) ≈ diff(values)
    end

    @testset "Multi-geotile aggregation with missing data" begin
        # Test aggregation when some geotiles have gaps
        n_dates = 12
        n_geotiles = 4

        dates = [DateTime(2020,1,1) + Month(i) for i in 0:n_dates-1]
        geotiles = ["gt1", "gt2", "gt3", "gt4"]

        # Create data with gaps
        data = randn(n_geotiles, n_dates) .* 5.0
        # Geotile 2 missing months 5-8
        data[2, 5:8] .= NaN
        # Geotile 4 missing months 3-6
        data[4, 3:6] .= NaN

        geotile_dim = DD.Dim{:geotile}(geotiles)
        date_dim = DD.Dim{:date}(dates)
        dv = DimArray(data, (geotile_dim, date_dim))

        # Aggregate (sum, ignoring NaN)
        dv_agg = [sum(filter(!isnan, dv[:, j])) for j in 1:n_dates]

        # Test: Aggregation should handle NaN correctly
        # Month 1: all 4 geotiles present
        @test !isnan(dv_agg[1])
        @test dv_agg[1] ≈ sum(data[:, 1])

        # Month 5: geotile 2 (months 5-8) *and* geotile 4 (months 3-6) are both missing, so only
        # geotiles 1 and 3 contribute
        @test !isnan(dv_agg[5])
        expected_m5 = data[1, 5] + data[3, 5]
        @test dv_agg[5] ≈ expected_m5

        # Test: No time should have NaN in aggregation (at least one geotile always present)
        @test all(.!isnan.(dv_agg))
    end

    @testset "Trend significance testing" begin
        # Test whether a trend is statistically significant
        Random.seed!(456)

        # Significant trend: -0.5 m/yr with small noise
        n = 40
        t = collect(0:n-1) ./ 12.0
        significant_trend = -0.5 .* t .+ 0.05 .* randn(n)

        # Insignificant trend: -0.05 m/yr with large noise
        insignificant_trend = -0.05 .* t .+ 0.5 .* randn(n)

        # Fit trends
        fit_sig = sum(t .* significant_trend) / sum(t .^ 2)
        fit_insig = sum(t .* insignificant_trend) / sum(t .^ 2)

        # Compute residuals and standard errors
        residuals_sig = significant_trend .- (fit_sig .* t)
        residuals_insig = insignificant_trend .- (fit_insig .* t)

        rmse_sig = sqrt(mean(residuals_sig .^ 2))
        rmse_insig = sqrt(mean(residuals_insig .^ 2))

        # Test: Significant trend has smaller RMSE relative to trend magnitude
        @test rmse_sig / abs(fit_sig) < rmse_insig / abs(fit_insig)

        # Test: Recovered trends should be close to true trends
        @test fit_sig ≈ -0.5 rtol=0.2
        @test abs(fit_insig) < 0.15  # Noise dominates
    end

    @testset "Regional summary statistics" begin
        # Test calculation of regional summary metrics
        Random.seed!(789)

        # Synthetic regional time series
        n_dates = 36
        t = collect(0:n_dates-1) ./ 12.0

        # Region 1: strong mass loss
        region1 = -0.6 .* t .+ 0.1 .* randn(n_dates)

        # Region 2: weak mass loss
        region2 = -0.2 .* t .+ 0.1 .* randn(n_dates)

        # Region 3: near balance
        region3 = -0.05 .* t .+ 0.1 .* randn(n_dates)

        # Compute summary statistics
        summary = Dict(
            :region1 => (mean=mean(region1), std=std(region1),
                        trend=(region1[end]-region1[1])/(t[end]-t[1])),
            :region2 => (mean=mean(region2), std=std(region2),
                        trend=(region2[end]-region2[1])/(t[end]-t[1])),
            :region3 => (mean=mean(region3), std=std(region3),
                        trend=(region3[end]-region3[1])/(t[end]-t[1]))
        )

        # Test: Regions should rank correctly by mass loss rate
        @test abs(summary[:region1].trend) > abs(summary[:region2].trend)
        @test abs(summary[:region2].trend) > abs(summary[:region3].trend)

        # Test: All regions should have negative trends
        @test summary[:region1].trend < 0
        @test summary[:region2].trend < 0
        # Region 3 might be slightly positive due to noise, but small
        @test abs(summary[:region3].trend) < 0.15
    end

    @testset "Edge case - single time point" begin
        # Test handling of single time point (no trend can be fit)
        n_geotiles = 3
        n_heights = 5

        dates = [DateTime(2020,1,1)]
        heights = collect(range(1500, 2500, length=n_heights))
        geotiles = ["gt1", "gt2", "gt3"]

        # Single snapshot
        dh_data = randn(n_geotiles, 1, n_heights)
        area_data = fill(10.0, n_geotiles, n_heights)

        geotile_dim = DD.Dim{:geotile}(geotiles)
        date_dim = DD.Dim{:date}(dates)
        height_dim = DD.Dim{:height}(heights)

        dh = DimArray(dh_data, (geotile_dim, date_dim, height_dim))
        area = DimArray(area_data, (geotile_dim, height_dim))

        # Should be able to compute volume
        dv = GVC.dh2dv_geotile(dh, area)
        @test length(dv[1, :]) == 1
        @test !isnan(dv[1, 1])

        # Cannot fit trend with single point
        # (This should be handled gracefully in actual code)
    end

    @testset "Edge case - constant time series" begin
        # Test trend fitting on constant time series (zero trend)
        n = 20
        values = fill(5.0, n)  # Constant
        t = collect(0:n-1) ./ 12.0

        # Fit trend
        n_points = length(t)
        sum_t = sum(t)
        sum_y = sum(values)
        sum_tt = sum(t .^ 2)
        sum_ty = sum(t .* values)

        slope = (n_points * sum_ty - sum_t * sum_y) / (n_points * sum_tt - sum_t^2)

        # Test: Trend should be zero
        @test abs(slope) < 1e-10
    end

    @testset "Minimum fractional error keeps the larger error" begin
        dvar = Dim{:varname}(["runoff", "dm"]); drgi = Dim{:rgi}([1, 2])
        dpar = Dim{:parameter}(["trend", "amplitude"]); derr = Dim{:error}([false, true])
        rf = zeros(dvar, drgi, dpar, derr)
        rf[varname=At("runoff"), rgi=At(1), parameter=At("trend"), error=At(false)] = 100.0
        rf[varname=At("runoff"), rgi=At(1), parameter=At("trend"), error=At(true)] = 5.0     # below 20%: raised
        rf[varname=At("runoff"), rgi=At(2), parameter=At("trend"), error=At(false)] = -50.0
        rf[varname=At("runoff"), rgi=At(2), parameter=At("trend"), error=At(true)] = 30.0    # above 20%: kept
        rf[varname=At("dm"), rgi=At(1), parameter=At("trend"), error=At(true)] = 7.0         # not in varnames
        GlacierHydrology.region_min_frac_error!(rf; fractional_error_min=0.2, varnames=["runoff"], parameters=["trend"])
        @test rf[varname=At("runoff"), rgi=At(1), parameter=At("trend"), error=At(true)] ≈ 20.0
        @test rf[varname=At("runoff"), rgi=At(2), parameter=At("trend"), error=At(true)] ≈ 30.0
        @test rf[varname=At("dm"), rgi=At(1), parameter=At("trend"), error=At(true)] ≈ 7.0
    end
end

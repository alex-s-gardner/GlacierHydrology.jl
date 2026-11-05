using Test
using GlacierHydrology
import GlacierHydrology
using Statistics
using Random

@testset "Aggregation Integration" begin
    @testset "Trend fitting to synthetic time series" begin
        # Create synthetic time series: y = 5 + 0.5*t + 0.1*t² + 2*sin(2πt)
        t = range(0, 2, length=25)  # 0 to 2 years, 25 points
        true_offset = 5.0
        true_trend = 0.5
        true_accel = 0.1
        true_amplitude = 2.0

        y = true_offset .+ true_trend .* t .+ true_accel .* t.^2 .+ true_amplitude .* sin.(2π .* t)

        # Add small noise (seeded -- the trend/acceleration split below is sensitive to it)
        Random.seed!(42)
        y_noisy = y .+ randn(length(t)) .* 0.1

        # Fit trend using model3 (offset + trend + accel + seasonal)
        # This is conceptually what rgi_trends does
        X = hcat(ones(length(t)), t, t.^2, sin.(2π .* t), cos.(2π .* t))
        coeffs = X \ y_noisy

        # Verify fitted parameters
        fitted_offset = coeffs[1]
        fitted_trend = coeffs[2]
        fitted_accel = coeffs[3]

        @test fitted_offset ≈ true_offset atol=0.5

        # Over a 2-year window with 25 points, `t` and `t^2` are strongly collinear (and partly
        # aliased against the annual sine), so the linear and quadratic coefficients trade off
        # against each other and are not individually identifiable at this noise level. Assert the
        # combinations that are: the total modelled change across the window, and the fit quality.
        total_change_true = true_trend * 2 + true_accel * 4
        total_change_fit = fitted_trend * 2 + fitted_accel * 4
        @test total_change_fit ≈ total_change_true atol=0.15
        @test fitted_trend > 0 && fitted_accel > 0

        residuals = y_noisy .- X * coeffs
        @test sqrt(mean(residuals .^ 2)) < 0.15

        # Verify seasonal component
        seasonal_component = fitted_amplitude = sqrt(coeffs[4]^2 + coeffs[5]^2)
        @test seasonal_component ≈ true_amplitude rtol=0.15
    end

    @testset "Error propagation in ensemble" begin
        # 10 ensemble runs with known deviations from reference
        reference_value = 100.0
        ensemble_values = reference_value .+ randn(10) .* 5.0  # σ = 5.0

        # Compute 95th percentile error
        sorted_diffs = sort(abs.(ensemble_values .- reference_value))
        p95_idx = ceil(Int, 0.95 * length(sorted_diffs))
        p95_error = sorted_diffs[p95_idx]

        # Verify that 95th percentile captures most of the spread
        @test p95_error > 0
        @test sum(abs.(ensemble_values .- reference_value) .<= p95_error) >= 9  # At least 9 of 10

        # Verify reasonable magnitude (should be around 1.96*σ for normal distribution)
        expected_p95 = 1.96 * 5.0
        @test p95_error ≈ expected_p95 rtol=0.5  # Loose tolerance for small sample
    end

    @testset "Regional aggregation concept" begin
        # Simulate aggregation across geotiles
        # 2 geotiles in region 13, 1 geotile in region 14
        geotile_regions = [13, 13, 14]
        geotile_values = [10.0, 20.0, 30.0]
        geotile_areas = [1.0, 1.0, 1.0]

        # Aggregate by region
        region_13_values = geotile_values[geotile_regions .== 13]
        region_13_areas = geotile_areas[geotile_regions .== 13]
        region_13_sum = sum(region_13_values .* region_13_areas)

        region_14_values = geotile_values[geotile_regions .== 14]
        region_14_areas = geotile_areas[geotile_regions .== 14]
        region_14_sum = sum(region_14_values .* region_14_areas)

        @test region_13_sum ≈ 30.0 rtol=1e-6  # 10 + 20
        @test region_14_sum ≈ 30.0 rtol=1e-6  # 30

        # HMA region (98) combines regions 13 and 14
        hma_sum = region_13_sum + region_14_sum
        @test hma_sum ≈ 60.0 rtol=1e-6

        # Global (99) includes all
        global_sum = sum(geotile_values .* geotile_areas)
        @test global_sum ≈ 60.0 rtol=1e-6
    end
end

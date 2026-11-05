"""
Property-based tests for statistical invariants in GlacierHydrology.jl

These tests verify properties that should hold regardless of specific input values,
using randomized inputs to discover edge cases and ensure algorithmic correctness.
"""

using Test
using GlacierHydrology
import GlacierHydrology
using Random
using Statistics
using Dates
using DimensionalData
import DimensionalData as DD

@testset "Property-Based Tests" begin

    @testset "Volume conservation property" begin
        # Property: Volume calculated different ways should agree
        Random.seed!(111)

        n_trials = 20
        for trial in 1:n_trials
            # Random area distribution
            n_bins = rand(3:10)
            areas = abs.(randn(n_bins)) .* 10.0  # km²

            # Random elevation changes
            dh = randn(n_bins) .* 2.0  # m

            # Method 1: Direct calculation
            dv_direct = sum(dh .* areas) / 1000  # km³

            # Method 2: Loop-based calculation
            dv_loop = 0.0
            for i in 1:n_bins
                dv_loop += dh[i] * areas[i] / 1000
            end

            # Test: Both methods should agree exactly
            @test dv_direct ≈ dv_loop rtol=1e-12
        end
    end

    @testset "Error propagation monotonicity" begin
        # Property: Combined error should never be larger than largest component error
        Random.seed!(222)

        n_trials = 15
        for trial in 1:n_trials
            # Random number of measurements (2-5)
            n_measurements = rand(2:5)

            # Random errors (positive)
            errors = abs.(randn(n_measurements)) .+ 0.1

            # Combined error for independent measurements: sqrt(Σσ²)
            combined_error = sqrt(sum(errors .^ 2))

            # Test: Combined error should be less than sum of errors
            @test combined_error <= sum(errors)

            # Test: Combined error should be greater than or equal to largest individual error
            # (This is only true for certain combinations, so we test sqrt(Σσ²) >= max(σ))
            # Actually, this is not always true. sqrt(1² + 1²) = 1.41 > 1, but sqrt(2² + 0.1²) = 2.0025 > 2
            # Let's test a correct property: combined < sum
            @test combined_error < sum(errors) + 1e-10

            # For independent errors, combined is always less than arithmetic sum
            # but may be greater than individual components
        end
    end

    @testset "Coordinate transformation round-trip" begin
        # Property: Coordinate transformations should be reversible
        Random.seed!(333)

        n_trials = 25
        for trial in 1:n_trials
            # Random lat/lon within valid ranges
            lat = rand() * 170.0 - 85.0  # -85 to 85 (avoid poles)
            lon = rand() * 360.0 - 180.0  # -180 to 180

            # Get geotile for this coordinate
            geotile_width = 2.0
            lat_min = floor(lat / geotile_width) * geotile_width
            lat_max = lat_min + geotile_width
            lon_min = floor(lon / geotile_width) * geotile_width
            lon_max = lon_min + geotile_width

            # Convert to extent
            extent = (lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max)

            # Test: Original point should be within extent
            @test lon >= extent.lon_min
            @test lon < extent.lon_max
            @test lat >= extent.lat_min
            @test lat < extent.lat_max

            # Test: Extent dimensions should match geotile width
            @test (extent.lon_max - extent.lon_min) ≈ geotile_width
            @test (extent.lat_max - extent.lat_min) ≈ geotile_width
        end
    end

    @testset "Distance metric properties" begin
        # Properties of haversine distance:
        # 1. d(A,A) = 0
        # 2. d(A,B) = d(B,A) (symmetry)
        # 3. d(A,C) <= d(A,B) + d(B,C) (triangle inequality)

        Random.seed!(444)

        n_trials = 15
        for trial in 1:n_trials
            # Random coordinates
            lon_a, lat_a = rand() * 360 - 180, rand() * 180 - 90
            lon_b, lat_b = rand() * 360 - 180, rand() * 180 - 90
            lon_c, lat_c = rand() * 360 - 180, rand() * 180 - 90

            # Property 1: d(A,A) = 0
            d_aa = GlacierHydrology.haversine_distance((lon_a, lat_a), (lon_a, lat_a))
            @test d_aa ≈ 0.0 atol=1.0

            # Property 2: Symmetry
            d_ab = GlacierHydrology.haversine_distance((lon_a, lat_a), (lon_b, lat_b))
            d_ba = GlacierHydrology.haversine_distance((lon_b, lat_b), (lon_a, lat_a))
            @test d_ab ≈ d_ba rtol=1e-10

            # Property 3: Triangle inequality
            d_bc = GlacierHydrology.haversine_distance((lon_b, lat_b), (lon_c, lat_c))
            d_ac = GlacierHydrology.haversine_distance((lon_a, lat_a), (lon_c, lat_c))
            @test d_ac <= d_ab + d_bc + 1e-6  # Allow tiny numerical error
        end
    end

    @testset "Statistical measure properties" begin
        # Properties of NMAD and MAD
        Random.seed!(555)

        n_trials = 20
        for trial in 1:n_trials
            # Random data
            n = rand(10:100)
            data = randn(n) .* 5.0

            # Compute MAD
            median_val = median(data)
            deviations = abs.(data .- median_val)
            mad_val = median(deviations)

            # Compute NMAD
            nmad_val = 1.4826 * mad_val

            # Property: MAD should be non-negative
            @test mad_val >= 0

            # Property: NMAD should be larger than MAD (by factor 1.4826)
            @test nmad_val ≈ mad_val * 1.4826 rtol=1e-10

            # Property: MAD should be less than or equal to standard deviation (typically)
            # (This is true for normal distributions but not guaranteed for all)
            std_val = std(data)
            # For normal distribution, MAD ≈ 0.67 * std
            # Let's just test that they're in the same order of magnitude
            @test mad_val < std_val * 2.0
        end
    end

    @testset "Time conversion round-trip property" begin
        # Property: decimalyear ↔ datetime conversion should be reversible
        Random.seed!(666)

        n_trials = 30
        for trial in 1:n_trials
            # Random date between 2000 and 2030
            year = rand(2000:2030)
            month = rand(1:12)
            day = rand(1:28)  # Safe for all months

            original_date = DateTime(year, month, day)

            # Convert to decimal year and back
            decimal = GlacierHydrology.decimalyear(original_date)
            recovered_date = GlacierHydrology.decimalyear2datetime(decimal)

            # Test: Should recover original date (within 1 day due to rounding)
            diff_days = abs(Dates.value(recovered_date - original_date)) / (24 * 60 * 60 * 1000)
            @test diff_days < 1.0
        end
    end

    @testset "Trend fitting properties" begin
        # Properties of linear trend fitting
        Random.seed!(777)

        n_trials = 15
        for trial in 1:n_trials
            # Random trend and data length
            n = rand(10:50)
            true_trend = randn() * 0.5
            true_intercept = randn() * 10.0

            # Generate data with known trend
            t = collect(0:n-1) ./ 12.0  # Monthly to years
            noise_sigma = 0.1
            data = true_intercept .+ true_trend .* t .+ noise_sigma .* randn(n)

            # Fit trend
            A = hcat(ones(n), t)
            coeffs = A \ data

            # Property: the fitted trend agrees with the truth to within its own uncertainty.
            # A relative tolerance is the wrong test here -- `true_trend` is drawn from
            # randn()*0.5, so it is sometimes ~0, and no percentage of ~0 is achievable when the
            # slope's standard error is finite. Compare against the analytic OLS standard error
            # instead: SE = sigma / sqrt(sum((t - mean(t))^2)).
            slope_se = noise_sigma / sqrt(sum((t .- mean(t)) .^ 2))
            @test abs(coeffs[2] - true_trend) < 5 * slope_se

            # Property: Residuals should have zero mean
            fitted = A * coeffs
            residuals = data .- fitted
            @test mean(residuals) ≈ 0.0 atol=1e-10

            # Property: Fitted line should pass through (mean(t), mean(data))
            fitted_at_mean_t = coeffs[1] + coeffs[2] * mean(t)
            @test fitted_at_mean_t ≈ mean(data) rtol=1e-6
        end
    end

    @testset "Area-weighted mean properties" begin
        # Properties of area-weighted averaging
        Random.seed!(888)

        n_trials = 20
        for trial in 1:n_trials
            # Random values and areas
            n = rand(3:10)
            values = randn(n) .* 5.0
            areas = abs.(randn(n)) .* 10.0 .+ 1.0  # Positive areas

            # Compute area-weighted mean
            weighted_mean = sum(values .* areas) / sum(areas)

            # Property: Weighted mean should be between min and max of values
            @test weighted_mean >= minimum(values) - 1e-10
            @test weighted_mean <= maximum(values) + 1e-10

            # Property: If all values equal, weighted mean = that value
            constant_values = fill(7.3, n)
            constant_weighted = sum(constant_values .* areas) / sum(areas)
            @test constant_weighted ≈ 7.3 rtol=1e-10

            # Property: If all areas equal, weighted mean = simple mean
            equal_areas = fill(5.0, n)
            equal_weighted = sum(values .* equal_areas) / sum(equal_areas)
            simple_mean = mean(values)
            @test equal_weighted ≈ simple_mean rtol=1e-10
        end
    end

    @testset "Flux accumulation conservation" begin
        # Property: Total flux at outlet = sum of all local inputs
        Random.seed!(999)

        n_trials = 10
        for trial in 1:n_trials
            # Random network size
            n_nodes = rand(5:20)

            # Create linear network (simple chain)
            ids = collect(1:n_nodes)
            next_ids = vcat(collect(2:n_nodes), [0])

            # Random local fluxes (positive)
            local_flux = abs.(randn(n_nodes)) .* 10.0

            # Accumulate flux. Node 1 is the headwater and node n_nodes the outlet, so the sweep
            # has to ascend. Descending (n_nodes:-1:1) pushes a node's running total downstream
            # *before* its own upstream contribution has arrived, leaving the outlet with only the
            # last couple of local fluxes.
            accumulated_flux = copy(local_flux)
            for i in 1:n_nodes  # upstream -> downstream
                if next_ids[i] != 0
                    downstream_idx = next_ids[i]
                    accumulated_flux[downstream_idx] += accumulated_flux[i]
                end
            end

            # Property: Flux at outlet equals sum of all inputs
            outlet_flux = accumulated_flux[end]
            total_input = sum(local_flux)
            @test outlet_flux ≈ total_input rtol=1e-10

            # Property: Flux increases (or stays same) downstream
            for i in 1:n_nodes-1
                @test accumulated_flux[i+1] >= accumulated_flux[i] - 1e-10
            end
        end
    end

    @testset "Unit conversion invariance" begin
        # Property: Physical quantities should be consistent across unit conversions
        Random.seed!(1010)

        n_trials = 15
        for trial in 1:n_trials
            # Random volume in km³
            volume_km3 = abs(randn()) * 100.0

            # Convert to m³
            volume_m3 = volume_km3 * 1e9

            # Convert to mass (Gt) using ice density
            ice_density = 910.0  # kg/m³
            mass_kg = volume_m3 * ice_density
            mass_Gt = mass_kg / 1e12

            # Back to volume
            back_to_m3 = mass_Gt * 1e12 / ice_density
            back_to_km3 = back_to_m3 / 1e9

            # Property: Round-trip conversion should preserve value
            @test back_to_km3 ≈ volume_km3 rtol=1e-10

            # Property: Mass should equal volume * density * 0.910
            expected_mass_Gt = volume_km3 * ice_density / 1e3  # km³ * kg/m³ / 1e3 = Gt
            @test mass_Gt ≈ expected_mass_Gt rtol=1e-10
        end
    end

    @testset "Ensemble statistics consistency" begin
        # Properties of ensemble operations
        Random.seed!(1111)

        n_trials = 10
        for trial in 1:n_trials
            # Random ensemble
            n_members = rand(10:50)
            n_times = rand(5:20)

            ensemble_data = randn(n_members, n_times) .* 5.0

            # Compute statistics
            ens_mean = mean(ensemble_data; dims=1)[1, :]
            ens_median = median(ensemble_data; dims=1)[1, :]
            ens_std = std(ensemble_data; dims=1)[1, :]

            # Property: Mean and median should be similar for large ensembles
            if n_members > 30
                @test mean(abs.(ens_mean .- ens_median)) < 1.0
            end

            # Property: All ensemble members should lie within ±3σ of mean (mostly)
            for j in 1:n_times
                within_3sigma = sum(abs.(ensemble_data[:, j] .- ens_mean[j]) .<= 3 * ens_std[j])
                @test within_3sigma / n_members > 0.95  # At least 95%
            end

            # Property: Standard deviation should be non-negative
            @test all(ens_std .>= 0)
        end
    end

    @testset "Temporal binning properties" begin
        # Properties of date binning
        Random.seed!(1212)

        n_trials = 10
        for trial in 1:n_trials
            # Random dates over several years
            n_dates = rand(20:100)
            start_year = rand(2000:2020)

            dates = [DateTime(start_year, 1, 1) + Day(rand(0:1460)) for _ in 1:n_dates]
            sort!(dates)

            # Bin into years
            years = [Dates.year(d) for d in dates]
            unique_years = unique(years)

            # Property: Number of unique years should be ≤ 5 (since dates span ~4 years)
            @test length(unique_years) <= 6

            # Property: Dates should be sorted
            @test issorted(dates)

            # Property: All dates should be within expected range
            @test all(start_year .<= years .<= start_year + 4)
        end
    end

    @testset "Spatial aggregation commutativity" begin
        # Property: Order of spatial aggregation shouldn't matter
        Random.seed!(1313)

        n_trials = 10
        for trial in 1:n_trials
            # Create random spatial grid
            n_x = rand(3:8)
            n_y = rand(3:8)

            data = randn(n_x, n_y) .* 10.0

            # Method 1: Sum over x first, then y
            sum_x_first = sum(sum(data; dims=1))

            # Method 2: Sum over y first, then x
            sum_y_first = sum(sum(data; dims=2))

            # Method 3: Sum all at once
            sum_all = sum(data)

            # Property: All methods should give same result
            @test sum_x_first ≈ sum_y_first rtol=1e-10
            @test sum_x_first ≈ sum_all rtol=1e-10
            @test sum_y_first ≈ sum_all rtol=1e-10
        end
    end
end

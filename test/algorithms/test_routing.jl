"""
Tests for river network routing in utilities_routing.jl

These tests validate river network traversal, flux accumulation,
hydrologic routing, and population exposure calculations.
"""

using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using DataFrames
using Statistics
using DimensionalData
import DimensionalData as DD

# Include test fixtures
include("../fixtures/synthetic_network.jl")

@testset "Routing Operations" begin
    @testset "trace_downstream" begin
        # Simple linear network: 1 → 2 → 3 → 0 (ocean)
        ids = [1, 2, 3]
        next_ids = [2, 3, 0]

        # Trace from node 1
        path = GGA.trace_downstream(1, ids, next_ids)
        @test path == [1, 2, 3]

        # Trace from node 2
        path2 = GGA.trace_downstream(2, ids, next_ids)
        @test path2 == [2, 3]

        # Trace from ocean terminus
        path3 = GGA.trace_downstream(3, ids, next_ids)
        @test path3 == [3]
    end

    @testset "trace_downstream - branching network" begin
        # Branching network: 1 → 3, 2 → 3 → 4 → 0
        ids = [1, 2, 3, 4]
        next_ids = [3, 3, 4, 0]

        # Trace from branch 1
        path1 = GGA.trace_downstream(1, ids, next_ids)
        @test path1 == [1, 3, 4]

        # Trace from branch 2
        path2 = GGA.trace_downstream(2, ids, next_ids)
        @test path2 == [2, 3, 4]

        # Trace from confluence
        path3 = GGA.trace_downstream(3, ids, next_ids)
        @test path3 == [3, 4]

        # All branches eventually reach node 4
        @test 4 in path1
        @test 4 in path2
        @test 4 in path3
    end

    @testset "haversine_distance" begin
        # Same point
        d_same = GGA.haversine_distance((0.0, 0.0), (0.0, 0.0))
        @test d_same ≈ 0.0 atol=1.0

        # Equator separation (1 degree ≈ 111.32 km)
        d_eq = GGA.haversine_distance((0.0, 0.0), (1.0, 0.0))
        @test d_eq ≈ 111320.0 rtol=1e-2

        # Antipodal points (half Earth circumference ≈ 20,000 km)
        d_antipodal = GGA.haversine_distance((0.0, 0.0), (180.0, 0.0))
        @test d_antipodal ≈ 20037000.0 rtol=1e-2

        # Pole to pole
        d_pole = GGA.haversine_distance((0.0, 90.0), (0.0, -90.0))
        @test d_pole ≈ 20004000.0 rtol=1e-2
    end

    # `flux_accumulate!` takes five positional arguments -- not a DataFrame with keyword column
    # names -- and mutates a DimArray whose *last* dimension is indexed by node id:
    #
    #   flux_accumulate!(river_inputs, id, nextdown_id, headbasin, majorbasin_id)
    #
    # (see the real callers in utilities_routing.jl and land_surface_model_routing.jl). `headbasin`
    # flags headwater nodes -- those nothing else drains into -- and `majorbasin_id` groups nodes
    # into independently-processed basins.
    function accumulate_flux(ids, next_ids, local_flux)
        flux = DimArray(reshape(collect(float.(local_flux)), 1, length(ids)),
                        (Dim{:Ti}(1:1), Dim{:id}(ids)))
        headbasin = [!(id in next_ids) for id in ids]
        GGA.flux_accumulate!(flux, ids, next_ids, headbasin, fill(1, length(ids)))
        return flux
    end

    @testset "flux_accumulate! - linear network" begin
        # Linear network: 1 → 2 → 3 → ocean
        network = create_linear_network(3)
        local_flux = [10.0, 5.0, 3.0]  # Gt/yr at each node

        flux = accumulate_flux(network.COMID, network.NextDownID, local_flux)

        # Node 3 (terminus) accumulates everything upstream
        @test flux[1, DD.At(3)] ≈ sum(local_flux)

        # Node 2 carries its own flux plus node 1
        @test flux[1, DD.At(2)] ≈ 10.0 + 5.0

        # Node 1 (headwater) keeps only its own flux
        @test flux[1, DD.At(1)] ≈ 10.0
    end

    @testset "flux_accumulate! - branching network" begin
        # Branching network: 1 → 3 ← 2, then 3 → 4 → 5 → ocean
        network = create_branching_network()
        local_flux = [8.0, 6.0, 2.0, 1.0, 0.5]

        flux = accumulate_flux(network.COMID, network.NextDownID, local_flux)

        # Confluence receives both branches plus its own local input
        @test flux[1, DD.At(3)] ≈ 8.0 + 6.0 + 2.0 rtol=1e-6

        # Outlet carries the whole network
        @test flux[1, DD.At(5)] ≈ sum(local_flux) rtol=1e-6

        # Flux never decreases downstream
        for i in 1:nrow(network)
            if network.NextDownID[i] != 0
                @test flux[1, DD.At(network.NextDownID[i])] >= flux[1, DD.At(network.COMID[i])]
            end
        end
    end

    @testset "linear_reservoir_impulse_response_monthly" begin
        # Test hydrologic routing with known residence time
        # `linear_reservoir_impulse_response_monthly(Tb)` takes a single scalar -- the baseflow
        # residence time in *days* -- and returns the monthly impulse response as a vector of
        # fractions. It is not a convolution routine: applying a response to a time series is
        # `apply_vector_impulse_resonse!`. Fractions are rounded to two digits and truncated at
        # the first zero, so they sum to 1 only to within that rounding.
        response = GGA.linear_reservoir_impulse_response_monthly(45)

        @test response isa AbstractVector
        @test !isempty(response)

        # Mass is (near) conserved
        @test sum(response) ≈ 1.0 atol=0.02

        # Every entry is a valid fraction
        @test all(0 .<= response .<= 1)

        # The tail decays: once past the peak the response is non-increasing
        peak = argmax(response)
        @test all(diff(response[peak:end]) .<= 1e-12)

        # The last month is a small remainder
        @test response[end] <= 0.05
    end

    @testset "linear_reservoir_impulse_response_monthly - residence time scaling" begin
        fast = GGA.linear_reservoir_impulse_response_monthly(15)
        slow = GGA.linear_reservoir_impulse_response_monthly(45)

        # A longer residence time spreads the response over more months...
        @test length(slow) > length(fast)

        # ...and lowers the fraction released in the first month
        @test slow[1] < fast[1]

        # both still conserve mass to within rounding
        @test sum(fast) ≈ 1.0 atol=0.02
        @test sum(slow) ≈ 1.0 atol=0.02
    end

    @testset "flux conservation in routing" begin
        # Test that total flux is conserved through routing
        network = create_dendritic_network()   # 1→4, 2→4, 3→5, 4→5, 5→6, 6→ocean

        local_flux = abs.(randn(nrow(network))) .* 5.0
        # nodes 1, 2 and 3 have nothing flowing into them
        @test [!(id in network.NextDownID) for id in network.COMID] ==
              [true, true, true, false, false, false]

        flux = accumulate_flux(network.COMID, network.NextDownID, local_flux)

        # Total flux at the outlet equals the sum of all local inputs
        outlet_id = network.COMID[findfirst(network.NextDownID .== 0)]
        @test flux[1, DD.At(outlet_id)] ≈ sum(local_flux) rtol=1e-10
    end

    @testset "network topology validation" begin
        # Test detection of circular references (should not occur)

        # Valid network
        valid_ids = [1, 2, 3]
        valid_next = [2, 3, 0]
        @test_nowarn GGA.trace_downstream(1, valid_ids, valid_next)

        # Circular network: 1 → 2 → 1 (invalid). `trace_downstream` guards against this with an
        # iteration cap rather than by raising: the kwarg is `maxiters` (not `max_steps`), and a
        # cycle simply terminates once the cap is reached.
        circular_ids = [1, 2]
        circular_next = [2, 1]

        cycle_path = GGA.trace_downstream(1, circular_ids, circular_next)
        @test length(cycle_path) <= length(circular_ids) + 1   # bounded, does not hang

        bounded = GGA.trace_downstream(1, circular_ids, circular_next; maxiters=10)
        @test length(bounded) <= 11
        @test all(in(circular_ids), bounded)
    end

    @testset "multiple outlets handling" begin
        # Test network with multiple ocean outlets (endorheic + coastal)

        # Network: 1→2→0, 3→4→0 (two separate drainage systems)
        ids = [1, 2, 3, 4]
        next_ids = [2, 0, 4, 0]

        # Trace from first system
        path1 = GGA.trace_downstream(1, ids, next_ids)
        @test path1 == [1, 2]
        @test 3 ∉ path1  # Should not include other drainage

        # Trace from second system
        path2 = GGA.trace_downstream(3, ids, next_ids)
        @test path2 == [3, 4]
        @test 1 ∉ path2  # Should not include other drainage
    end

    @testset "distance calculations along network" begin
        # Test cumulative distance calculation
        network = create_linear_network(5)

        # Each segment is 10 km
        @test all(network.lengthkm .≈ 10.0)

        # Calculate cumulative distance from headwater
        total_distance = sum(network.lengthkm)
        @test total_distance ≈ 50.0  # 5 segments × 10 km

        # Verify haversine distances match network distances (roughly)
        for i in 1:(nrow(network)-1)
            lon1, lat1 = network.lon[i], network.lat[i]
            lon2, lat2 = network.lon[i+1], network.lat[i+1]

            haversine_km = GGA.haversine_distance((lon1, lat1), (lon2, lat2)) / 1000
            network_km = network.lengthkm[i]

            # Should be similar (within 20% for straight line vs river meander)
            @test haversine_km ≈ network_km rtol=0.3
        end
    end

    @testset "flux routing with seasonal variation" begin
        # Routing a time series is a convolution of the series with the impulse response, done by
        # `apply_vector_impulse_resonse!`. `linear_reservoir_impulse_response_monthly` only builds
        # the response kernel from a residence time.
        n_months = 24

        # Seasonal input (sine wave), shaped (1, 1, n_months) as apply_vector_impulse_resonse!
        # expects a 3D array whose third axis is time
        t = collect(0:n_months-1) ./ 12.0  # Years
        input_flux = 10.0 .+ 5.0 .* sin.(2π .* t)  # Mean 10, amplitude 5

        M = reshape(copy(input_flux), 1, 1, n_months)
        response = GGA.linear_reservoir_impulse_response_monthly(60)
        GGA.apply_vector_impulse_resonse!(M, response)
        output_flux = vec(M[1, 1, :])

        # Test: Output should be smoother than input (reservoir dampens)
        @test std(output_flux[6:end]) < std(input_flux)

        # Test: Mean should be roughly preserved (the kernel sums to ~1)
        @test mean(output_flux[6:end]) ≈ mean(input_flux[6:end]) rtol=0.15

        # Test: the seasonal peak is delayed relative to the input
        @test argmax(output_flux[1:12]) > argmax(input_flux[1:12])
    end

    @testset "Edge case - zero flux" begin
        # Test handling of zero flux
        network = create_linear_network(3)
        flux = accumulate_flux(network.COMID, network.NextDownID, zeros(3))

        # All fluxes should remain zero
        @test all(parent(flux) .== 0.0)
    end

    @testset "Edge case - single node network" begin
        # Single-node network (direct to ocean): nothing to accumulate
        network = DataFrame(COMID=[1], NextDownID=[0], lengthkm=[10.0], lon=[-120.0], lat=[45.0])
        flux = accumulate_flux(network.COMID, network.NextDownID, [15.0])

        # Flux should be unchanged
        @test flux[1, DD.At(1)] == 15.0
    end

    @testset "Large network scaling" begin
        # Test with larger network to ensure performance
        n_nodes = 100

        # Create chain network
        ids = collect(1:n_nodes)
        next_ids = vcat(collect(2:n_nodes), [0])

        # Should complete quickly
        local flux
        @test_nowarn flux = accumulate_flux(ids, next_ids, ones(n_nodes))

        # Test: Outlet should have sum of all inputs
        @test flux[1, DD.At(n_nodes)] ≈ n_nodes rtol=1e-10

        # and flux grows by exactly one unit per step down the chain
        @test [flux[1, DD.At(i)] for i in ids] ≈ collect(1.0:n_nodes) rtol=1e-10
    end

    @testset "glacier river products are named by glacier inventory" begin
        # The RGI6 names are the ones the existing products were written under.
        @test GGA.glacier_rivers_file("glacier") ==
              joinpath(GGA.pathlocal[:river], "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier.arrow")
        @test GGA.glacier_routing_file("glacier") ==
              replace(GGA.pathlocal[:glacier_individual], ".gpkg" => "_routing.arrow")
        @test GGA.glacier_rivers_file("glacier_rgi7") ==
              joinpath(GGA.pathlocal[:river], "riv_pfaf_MERIT_Hydro_v07_Basins_v01_glacier_rgi7.arrow")
        @test GGA.glacier_routing_file("glacier_rgi7") ==
              replace(GGA.pathlocal[:glacier_rgi7_individual], ".gpkg" => "_routing.arrow")
        @test GGA.glacier_routing_file("glacier_rgi7") != GGA.glacier_routing_file("glacier")
    end
end

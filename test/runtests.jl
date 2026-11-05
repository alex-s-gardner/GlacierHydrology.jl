using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA

@testset "GlobalGlacierAnalysis.jl Tests" begin
    # Include test helpers
    include("test_helpers.jl")

    # Unit tests - Pure mathematical functions
    @testset "Unit Tests" begin
        include("unit/test_array_ops.jl")
        include("unit/test_project.jl")
    end

    # Algorithm tests - Complex operations with synthetic data
    @testset "Algorithm Tests" begin
        include("algorithms/test_gemb.jl")
        include("algorithms/test_gemb_tiles.jl")
        include("algorithms/test_gemb_physics.jl")
        include("algorithms/test_routing.jl")
        include("algorithms/test_discharge.jl")
        include("algorithms/test_postprocessing.jl")
        include("algorithms/test_regional_units.jl")
        include("algorithms/test_wgms.jl")
    end

    # I/O tests - Data persistence and round-trip
    @testset "I/O Tests" begin
        include("io/test_dimarray.jl")
        include("io/test_jld2.jl")
        include("io/test_netcdf.jl")
        include("io/test_readers.jl")
    end

    # Integration tests - End-to-end workflows
    @testset "Integration Tests" begin
        include("integration/test_aggregation.jl")
        include("integration/test_properties.jl")
    end
end

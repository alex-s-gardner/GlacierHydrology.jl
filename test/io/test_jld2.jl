using Test
using GlacierHydrology
import GlacierHydrology
using JLD2
using DimensionalData
import DimensionalData as DD

@testset "JLD2 I/O" begin
    @testset "DimArray save and load" begin
        mktempdir() do tmpdir
            filepath = joinpath(tmpdir, "test_dimarray.jld2")

            # Create DimArray with custom dimensions
            data = randn(3, 4, 5)
            geotiles = ["lat+45+47lon-123-121", "lat+47+49lon-123-121", "lat+49+51lon-123-121"]
            dates = [DateTime(2018,1,1), DateTime(2018,7,1), DateTime(2019,1,1), DateTime(2019,7,1)]
            heights = [1000, 1500, 2000, 2500, 3000]

            da = DimArray(
                data,
                (
                    DD.Dim{:geotile}(geotiles),
                    DD.Dim{:date}(dates),
                    DD.Dim{:height}(heights)
                )
            )

            # Save
            JLD2.jldsave(filepath; data=da)

            # Load
            loaded = JLD2.load(filepath, "data")

            # Verify structure preserved
            @test size(loaded) == size(da)
            @test dims(loaded, :geotile) isa DD.Dim{:geotile}
            @test dims(loaded, :date) isa DD.Dim{:date}
            @test dims(loaded, :height) isa DD.Dim{:height}

            # Verify values
            @test all(loaded .== da)
        end
    end

    @testset "Nested dictionary save and load" begin
        mktempdir() do tmpdir
            filepath = joinpath(tmpdir, "test_nested.jld2")

            # Create nested structure
            da1 = DimArray(ones(5, 5), (X=1:5, Y=1:5))
            da2 = DimArray(ones(5, 5) .* 2, (X=1:5, Y=1:5))

            nested = Dict(
                "runs" => Dict(
                    "run1" => da1,
                    "run2" => da2
                )
            )

            # Save
            JLD2.jldsave(filepath; data=nested)

            # Load
            loaded = JLD2.load(filepath, "data")

            # Verify structure
            @test haskey(loaded, "runs")
            @test haskey(loaded["runs"], "run1")
            @test haskey(loaded["runs"], "run2")

            # Verify values
            @test all(loaded["runs"]["run1"] .== da1)
            @test all(loaded["runs"]["run2"] .== da2)
        end
    end
end

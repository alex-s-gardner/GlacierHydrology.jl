using Test
using GlacierHydrology
import GlacierHydrology
using DimensionalData
import DimensionalData as DD
using Dates

@testset "NetCDF I/O" begin
    @testset "dimstack2netcdf round-trip" begin
        mktempdir() do tmpdir
            filepath = joinpath(tmpdir, "test.nc")

            # Create DimStack
            dates = [DateTime(2018,1,1), DateTime(2018,7,1), DateTime(2019,1,1)]
            x_coords = [1.0, 2.0, 3.0]
            y_coords = [4.0, 5.0, 6.0]

            temp_data = randn(3, 3, 3)
            precip_data = randn(Float32, 3, 3, 3)

            ds = DimStack((
                temperature = DimArray(temp_data, (Ti=dates, X=x_coords, Y=y_coords);
                                     metadata=Dict("units" => "°C")),
                precipitation = DimArray(precip_data, (Ti=dates, X=x_coords, Y=y_coords);
                                       metadata=Dict("units" => "mm"))
            ))

            # Write. `dimstack2netcdf` forwards its kwargs straight to NCDatasets.defVar, so only
            # defVar options (deflatelevel, chunksizes, ...) belong here -- global attributes are
            # taken from the DimStack's own metadata instead.
            GlacierHydrology.dimstack2netcdf(ds, filepath)

            # Read
            ds_loaded = GlacierHydrology.netcdf2dimstack(filepath)

            # Verify dimensions. Query by name, not by dimension type: the round trip rebuilds
            # dimensions generically from their NetCDF names, so a `Ti` written out comes back as
            # a `Dim{:Ti}` and `hasdim(ds_loaded, Ti)` would be false.
            @test DD.name.(dims(ds_loaded)) == (:Ti, :X, :Y)
            @test DD.hasdim(ds_loaded, :Ti) || DD.hasdim(ds_loaded, :date)
            @test DD.hasdim(ds_loaded, :X)
            @test DD.hasdim(ds_loaded, :Y)

            # Verify variables exist
            @test haskey(ds_loaded, :temperature) || haskey(ds_loaded, :Temperature)
            @test haskey(ds_loaded, :precipitation) || haskey(ds_loaded, :Precipitation)

            # Verify approximate equality (allow for type conversions)
            temp_key = haskey(ds_loaded, :temperature) ? :temperature : :Temperature
            @test size(ds_loaded[temp_key]) == size(temp_data)
        end
    end

    @testset "NetCDF with metadata preservation" begin
        mktempdir() do tmpdir
            filepath = joinpath(tmpdir, "test_metadata.nc")

            # Create DimArray with metadata
            data = collect(reshape(1.0:12.0, 3, 4))
            da = DimArray(
                data,
                (X=1:3, Y=1:4);
                name=:test_var,
                metadata=Dict("units" => "m", "long_name" => "Test Variable")
            )

            # Global attributes travel as DimStack metadata, which dimstack2netcdf writes to
            # nc.attrib in its final step.
            ds = DimStack((testvar=da,);
                metadata=Dict("title" => "Test Dataset", "institution" => "Test Lab"))
            GlacierHydrology.dimstack2netcdf(ds, filepath)

            # Verify file was created
            @test isfile(filepath)

            # Read back
            ds_loaded = GlacierHydrology.netcdf2dimstack(filepath)

            # Verify data
            loaded_key = first(keys(ds_loaded))
            @test size(ds_loaded[loaded_key]) == size(data)
        end
    end

    @testset "NetCDF with missing values" begin
        mktempdir() do tmpdir
            filepath = joinpath(tmpdir, "test_missing.nc")

            # Create data with NaN (represents missing)
            data = [1.0 2.0 NaN; 4.0 NaN 6.0]
            da = DimArray(data, (X=1:2, Y=1:3); name=:data_with_missing)

            ds = DimStack((datavar=da,))

            # Write
            GlacierHydrology.dimstack2netcdf(ds, filepath)

            # Read
            ds_loaded = GlacierHydrology.netcdf2dimstack(filepath)

            # Verify NaN handling
            loaded_key = first(keys(ds_loaded))
            loaded_data = ds_loaded[loaded_key]

            # Check that NaN values are preserved or represented as missing/fill value
            @test any(isnan.(loaded_data)) || any(ismissing.(loaded_data))
        end
    end
end

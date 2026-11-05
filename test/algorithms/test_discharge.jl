"""
Tests for glacier discharge in utilities_synthesis.jl: per-geotile aggregation, the Gt to km³
ice-equivalent conversion, how the geotile discharge column accumulates, and the Antarctic
SMB-equilibrium discharge.
"""

using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using Dates
using DimensionalData
import DimensionalData as DD
using DataFrames

@testset "Discharge" begin
    @testset "discharge2geotile sums by geotile and converts Gt to km³ ice" begin
        dgeotile = Dim{:geotile}(["lat[+60+62]lon[-142-140]", "lat[+62+64]lon[-142-140]"])
        discharge = DataFrame(longitude=[-141.0, -141.5, -140.5, 10.0], latitude=[61.0, 61.5, 63.0, 0.0],
                              discharge_gtyr=[2.0, 3.0, 7.0, 100.0])

        gt = GGA.discharge2geotile(discharge, dgeotile)
        @test DD.metadata(gt)["units"] == "Gt/yr"
        @test gt[:discharge][geotile=At(dgeotile[1])] == 5.0
        @test gt[:discharge][geotile=At(dgeotile[2])] == 7.0

        km3 = GGA.discharge2geotile(discharge, dgeotile; mass2volume=true)
        @test DD.metadata(km3)["units"] == "km3/yr"
        @test km3[:discharge][geotile=At(dgeotile[1])] ≈ 5.0 * 1000 / GGA.δice
        @test km3[:discharge][geotile=At(dgeotile[2])] ≈ 7.0 * 1000 / GGA.δice
    end

    @testset "geotile discharge column is D·1000/δice km³ per year, accumulated" begin
        D = 4.0 # Gt/yr
        id = "lat[+60+62]lon[-142-140]"
        dgeotile = Dim{:geotile}([id])
        ddate = Dim{:date}(Date(2000, 1, 15):Month(1):Date(2001, 12, 15))
        dheight = Dim{:height}(50.0:100.0:450.0)
        dpscale = Dim{:pscale}([1.0, 2.0])
        dΔT = Dim{:ΔT}([0.0, 1.0])

        mktempdir() do temp_dir
            binned_synthesized_file = joinpath(temp_dir, "synthesized.jld2")
            GGA.FileIO.save(binned_synthesized_file, Dict("dh_hyps" => zeros(dgeotile, ddate, dheight)))

            area_km2 = DimArray(fill(2.0, 1, length(dheight)), (dgeotile, dheight))
            gemb = DimStack((; smb=zeros(dgeotile, ddate, dpscale, dΔT), dv=zeros(dgeotile, ddate, dpscale, dΔT),
                                fac=zeros(dgeotile, ddate, dpscale, dΔT)))
            gemb_fit = DataFrame(id=[id], pscale=[1.0], ΔT=[0.0], group=[1])
            geotiles0 = DataFrame(id=[id])
            for i in 1:19
                geotiles0[!, "rgi$i"] = [i == 1 ? 1.0 : 0.0]
            end

            discharge = GGA.discharge2geotile(DataFrame(longitude=[-141.0], latitude=[61.0], discharge_gtyr=[D]),
                                              dgeotile; mass2volume=true)
            out = GGA.individual_geotile_synthesis_gembfit_dv(binned_synthesized_file, gemb, gemb_fit,
                                                             discharge, geotiles0, area_km2)

            Δdecyear = GGA.decimalyear.(ddate) .- GGA.decimalyear(first(ddate))
            @test out[1, :discharge] ≈ D * 1000 / GGA.δice .* Δdecyear
            @test colmetadata(out, "discharge", "units") == "km3 [i.e.]"

            # a Gt/yr discharge is rejected rather than read as km³
            discharge_gt = GGA.discharge2geotile(DataFrame(longitude=[-141.0], latitude=[61.0], discharge_gtyr=[D]), dgeotile)
            @test_throws "discharge must be in km3/yr" GGA.individual_geotile_synthesis_gembfit_dv(
                binned_synthesized_file, gemb, gemb_fit, discharge_gt, geotiles0, area_km2)
        end
    end
    @testset "discharge2smb: an SMB trend of s m w.e. yr⁻¹ over A km² gives s·A/1000 Gt yr⁻¹" begin
        ddate = Dim{:date}(Date(1975, 1, 15):Month(1):Date(2004, 12, 15))
        Δdecyear = GGA.decimalyear.(ddate) .- GGA.decimalyear(first(ddate))

        s_we = [0.3, 1.2]    # m w.e. yr⁻¹
        A = [150.0, 40.0]    # km²
        # `smb` is cumulative, in metres of ice equivalent
        smb = [DimArray(s * 1000 / GGA.δice .* Δdecyear, ddate) for s in s_we]
        glaciers = DataFrame(CenLat=[-70.0, -65.0, 61.0], CenLon=[10.0, 20.0, -141.0],
                             smb=[smb; [DimArray(5.0 .* Δdecyear, ddate)]],
                             area_km2=[[A[1] / 2, A[1] / 2], [A[2]], [10.0]])

        discharge = GGA.discharge2smb(glaciers)

        @test nrow(discharge) == 2  # the glacier north of 60°S is left to observations
        @test discharge.latitude == [-70.0, -65.0]
        @test discharge.discharge_gtyr ≈ s_we .* A ./ 1000 rtol = 1e-6
    end
end

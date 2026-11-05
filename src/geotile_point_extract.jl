using JLD2
using FileIO
using DataFrames
import GlobalGlacierAnalysis as GGA
import GeometryOps as GO
using CairoMakie
using CSV

location_name = "PlaceGlacier"
x = -122.6360
y = 50.3997

reference_run = "binned/2deg/glacier_dh_best_nmad5_v01_filled_ac_p2_synthesized.jld2"
outfile_suffix = "_dT"  # ASCII stand-in for ΔT: this ends up in a filename

paths = GGA.pathlocal
path2reference = joinpath(paths[:data_dir], reference_run) 
binned_synthesized_dv_file = GGA.gemb_physics_path(replace(path2reference, ".jld2" => "_gembfit_dv$(outfile_suffix).jld2"),
                                                   GGA.gemb_reference_physics_member)

geotiles0 = FileIO.load(binned_synthesized_dv_file, "geotiles")

index = findfirst(GO.intersects.(Ref((x, y)), GGA.extent2rectangle.(geotiles0[:, :extent])))
dates = colmetadata(geotiles0, "dm", "date")

# `dm` is Gt and `mie2cubickm` is area_km2/1000, so Gt / (km²/1000) × 1000 is Gt km⁻² × 10⁶ = kg m⁻²
dm = geotiles0[index, :dm] ./ geotiles0[index, "mie2cubickm"] .* 1000
valid = .!isnan.(dm)

df = DataFrame(
    "datetime" => dates[valid],
    "mass change [kg m⁻²]" => dm[valid]
)

DataFrames.metadata!(df, "location_name", location_name,style=:note)
DataFrames.metadata!(df, "longitude", x, style=:note)
DataFrames.metadata!(df, "latitude", y, style=:note)
DataFrames.metadata!(df, "units", "kg m⁻²", style=:note)
DataFrames.metadata!(df, "source",reference_run, style=:note)

output = joinpath(paths[:project_dir], "Gardner2025_$location_name.csv")
CSV.write(output, df)
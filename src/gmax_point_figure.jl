# ============================================================================
# gmax_point_figure.jl
#
# Calculate and visualize the climatological maximum monthly fractional contribution of glacier runoff to total river flux (gmax).
#
# This script:
# 1. Loads glacier runoff and river flux data from NetCDF files
# 2. Processes time series data for specific geographic locations
# 3. Creates visualizations showing:
#    - River flux components (land, snow, glacier)
#    - Glacier runoff as a fraction of total river flux
#    - The maximum monthly contribution (gmax) with seasonal patterns
# 4. Saves the resulting figures for each location
#
# The script handles data from 2000-2024 and focuses on key river basins near glaciated regions.
# ============================================================================

begin
    using NCDatasets
    using GeoDataFrames
    using DimensionalData
    import GlacierHydrology
    using Dates
    using Statistics
    using CairoMakie
    using Unitful

    dates4plot = (Date(2000,1,1), Date(2024,10,1))
    separate_out_snow = false;

    paths = GlacierHydrology.pathlocal

    glacier_flux_path = joinpath(paths[:project_dir], "Gardner2025_glacier_summary_riverflux.nc")

    surface_mask = GlacierHydrology.reference_surface_mask()
    glacier_rivers_path = GlacierHydrology.glacier_rivers_file(surface_mask)
    glacier_rivers_land_flux_path = replace(glacier_rivers_path, ".arrow" => "_Qs_acc_Qsb_acc.nc")
    glacier_rivers_snow_flux_path = replace(glacier_rivers_path, ".arrow" => "_Qsm_acc.nc")
    land_path = paths[:river]
end

# Dates.unix2datetime(mtime(glacier_flux_path))

begin #[50s]
    # load data from netcdf files
    glacier_flux = NCDataset(glacier_flux_path)
    land_flux = NCDataset(glacier_rivers_land_flux_path)
    separate_out_snow && (snow_flux = NCDataset(glacier_rivers_snow_flux_path))

    # convert to DimensionalData
    glacier_flux = GlacierHydrology.nc2dd(glacier_flux["runoff"]) / (1000 * u"kg/m^3")
    land_flux = GlacierHydrology.nc2dd(land_flux["flux"])
    separate_out_snow && (snow_flux = GlacierHydrology.nc2dd(snow_flux["flux"]))
    # having issues with sortslices
    #@time glacier_flux = sortslices(glacier_flux, dims=:COMID)
    #@time land_flux = sortslices(land_flux, dims=:COMID)

    p = sortperm(dims(glacier_flux,:COMID).val)
    glacier_flux = glacier_flux[:,p]

    p = sortperm(dims(land_flux,:COMID).val)
    land_flux = land_flux[:,p]

    separate_out_snow && (p = sortperm(dims(snow_flux,:COMID).val))
    separate_out_snow && (snow_flux = snow_flux[:,p])

    # do this little hack to get COMID as ForwardOrdered
    land_flux = DimArray(parent(land_flux), (dims(land_flux, :Ti), Dim{:COMID}(val(dims(land_flux, :COMID).val))))
    separate_out_snow && (snow_flux = DimArray(parent(snow_flux), (dims(snow_flux, :Ti), Dim{:COMID}(val(dims(snow_flux, :COMID).val)))))

    #TODO: WHY DOES RUNOFF NOT GO TO ZERO? This needs further investigation... 
    # COMID = 81020092 

    locations = [
        (45047742, "Karachi, Pakistan"),
        (71017956, "Edmonton, Canada"),
        (71031599, "Calgary, Canada"),
        (71039426, "Lethbridge, Canada"),
        (78013077, "Vancouver, Canada"),
        (78016783, "Tacoma, USA"),
        (81014312, "Fairbanks, USA"),
        (67007399, "Lima, Peru"),
        (66010882, "Copiapó, Chile"),
        (66015561, "Temuco, Chile"),
        (65003342, "Mendoza, Argentina"),
        (66013257, "Santiago, Chile"),
        (65001686, "San Juan, Argentina"),
        (67007407, "Lima, Peru"), 
        (45011579, "Lahore, Pakistan"),
        (45008062, "Jammu, India"),
        (45022890, "Delhi, India"),
        (45043858, "Kanpur, India "),
        (45004712, "Peshawar, Pakistan"),
        (46044113, "Dushanbe, Tajikistan"),
        (46047000, "Tashkent, Uzbekistan"),
        (46041153, "Samarkand, Uzbekistan"),
        (46037217, "Tashkent, Uzbekistan"), 
        (47017377, "Ürümqi, China"),
        (47032572, "Jiayuguan, China"),
        (21001926, "Aosta, Italy"),    
        (21001122, "Geneva, Switzerland"),
        (21001840, "Lyon, France"),
        (21002710, "Grenoble, France"),
        (22027761, "Innsbruck, Austria"),
        (24003651, "Luleå, Sweden"),
        (25005999, "Mo i Rana, Norway")
    ]

    date_interval = dates4plot[1]..dates4plot[2]
end

for location in locations

    COMID = location[1]
    name = location[2]

    date_range = (DateTime(2000,3,1), DateTime(2024,10,1))
    fig = GlacierHydrology.plot_point_location_river_flux(land_flux, glacier_flux; date_range, COMID, name, show_title=false)
    display(fig)

    out_name = replace(split(name, ",")[1], " " => "_")
    
    output_path = joinpath(paths[:figures], "riverflux_$out_name.png")
    save(output_path, fig)

    println(output_path)
end

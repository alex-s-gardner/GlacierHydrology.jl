# land_surface_model_routing.jl
#
# Process and route land surface model runoff through river networks.
#
# This script:
# 1. Loads GLDAS land surface model data (CLSM, VIC, NOAH)
# 2. Processes surface runoff, subsurface runoff, and snowmelt
# 3. Routes water through river networks using a linear reservoir model
# 4. Calculates river discharge for glacier-connected river reaches
# 5. Validates results against independent river discharge estimates
# 6. Saves processed data as NetCDF files for downstream analysis
#
# The routing accounts for time delays in subsurface runoff using a 45-day 
# concentration time parameter following Getirana et al. (2012).

begin
    import GlobalGlacierAnalysis as GGA
    import GlacierVolumeChange as GVC
    using NCDatasets
    using Rasters
    using DataFrames
    import GeoInterface as GI
    using CairoMakie
    import GeometryOps as GO
    using GeoDataFrames
    using SortTileRecursiveTree
    using ProgressMeter
    using Rasters.Lookups
    using Dates
    using Statistics
    using NCDatasets
    using FileIO

    # YOU NEED YOU NEED TO MAKE SURE YOU ARE LOGGED INTO EARTHDATA
    # AND HAVE THE NETRC FILE SET UP

    # GLDAS LSM folder
    gldas_folder = "/mnt/devon-r2/data/GLDAS/"
    download_lsm_files = false
    lsm_names = ["GLDAS_CLSM10", "GLDAS_VIC10", "GLDAS_NOAH10"]

    #grid resolution
    grid_res = 1.

    paths = GGA.pathlocal
    rivers_path = joinpath(paths[:river])
    rivers_paths = GGA.allfiles(rivers_path; fn_startswith="riv_pfaf", fn_endswith= "MERIT_Hydro_v07_Basins_v01.shp")
    basins_paths = GGA.allfiles(rivers_path; fn_startswith="cat_pfaf", fn_endswith=".shp")
    # The glacier-fed reaches of the reference member's glacier inventory, from `GGA.glacier_routing`.
    surface_mask = GGA.reference_surface_mask()
    glacier_rivers_path = GGA.glacier_rivers_file(surface_mask)

    glacier_rivers_runoff_path = replace(glacier_rivers_path, ".arrow" => "_runoff.arrow")
    glacier_rivers_runoff_path_jld2 = replace(glacier_rivers_runoff_path, ".arrow" => "_runoff.jld2")

    glacier_routing_path = GGA.glacier_routing_file(surface_mask)
   
    # Use uncorrected ensemble simulations only
    river_Q_path_ens = joinpath(paths[:data_dir], "rivers/Collins2024/Qout_pfaf_ii_GLDAS_ENS_M_1980-01_2009-12_utc/")
    river_Q_files = readdir(river_Q_path_ens; join=true)

    # Total runoff is the sum of subsurface runoff "Qsb_tavg" and surface runoff "Qs_tavg".
    runoff_vars = [[:Qs_acc, :Qsb_acc],  [:Qsm_acc]]
    glacier_rivers_land_flux_paths = [
        replace(glacier_rivers_path, ".arrow" => "_Qs_acc_Qsb_acc.nc"),
        replace(glacier_rivers_path, ".arrow" => "_Qsm_acc.nc")]
end

# download GLDAS LSMs
if download_lsm_files
    # find list of urls that were downloaded from https://disc.gsfc.nasa.gov/datasets?keywords=GLDAS
    gldas_file_urls = GGA.allfiles(gldas_folder; fn_endswith = "_.txt")
    downloadstreams = 6;

    # download the files
    for file_urls in gldas_file_urls
    #file_urls = last(gldas_file_urls)

        # `-i` must stay last: it takes a filename argument, so any flag placed after it is
        # swallowed and that flag's value is parsed as a URI.
        # `GVC.aria2c_cmd` wraps the Aria2_jll binary rather than whatever `aria2c` is on PATH, and
        # works around the OpenSSL `legacy` provider failure it documents.
        cmd = GVC.aria2c_cmd(`--max-tries=10 --retry-wait=1 -c -k 1M -j 16 -x $downloadstreams -s $downloadstreams -d $gldas_folder -i $file_urls`)

        println(cmd)
        run(cmd)
    end
end

# load basin and river data
begin #[~2 min]
    # load in basins and find the centroids and area of each basin
    basins = fill(DataFrame(),length(basins_paths))
    Threads.@threads for i in eachindex(basins_paths) # [~1 min] 
        basins[i] = GeoDataFrames.read(basins_paths[i])
    end
    basins = reduce(vcat, basins)
    sort!(basins, [:COMID])
    basins.COMID = Int64.(basins.COMID) # this is needed as the COMID Type is Vector{Union{Missing, Int64}}:

    # load in river reaches
    rivers = GGA.river_reaches(rivers_paths; col_names=["geometry", "COMID", "NextDownID", "maxup"])
    
    if nrow(rivers) == nrow(basins)
        basins[!,:NextDownID] = rivers[:,:NextDownID]
        basins[!,:HeadBasin] = rivers[:, :maxup] .== 0
    else
        # add next downstream id and head basin flag (head basin = no river inputs)
        basins[!,:NextDownID] = zeros(eltype(rivers[1, :NextDownID]), length(basins.COMID))
        basins[!,:HeadBasin] = falses(length(basins.COMID))
        
        Threads.@threads for r in eachrow(basins)
            ind = searchsortedfirst(rivers.COMID,r.COMID)
            r.NextDownID = rivers[ind, :NextDownID]
            r.HeadBasin = rivers[ind, :maxup] .== 0
        end
    end
    
    # map cells to each downstream rivers (i.e. get list of inputs to each river)
    basins[!, :basin02] = floor.(Int16, basins.COMID./1000000)
    basins[!, :centroid] = GO.centroid.(basins.geometry)
end


# load gridded runoff data
#begin # [~6 min]
for i in eachindex(runoff_vars)

    # read in fluxes from each lsm
    files = GGA.allfiles(gldas_folder; fn_endswith=".nc4", fn_startswith = lsm_names[1])
    Q = Dict()

    for runoff_var in runoff_vars[i]
 
        Q[runoff_var] = cat(Raster.(files; name=runoff_vars[i][1])...; dims=Ti)
    
        for lsm_name in lsm_names
            files = GGA.allfiles(gldas_folder; fn_endswith=".nc4", fn_startswith = lsm_name)
        
            if (lsm_name == lsm_names[1]) && (runoff_var == runoff_vars[i][1])
                continue
            else
                Q[runoff_var] .+= cat(Raster.(files; name=runoff_var)...; dims=Ti)
            end
        end

        # take average accross lsms 
        Q[runoff_var] ./= length(lsm_names)

        # devide by 1000 to go from kg m-2 to m m-2
        Q[runoff_var] ./= 1000

        # For accumulated variables such as Qs_acc, the monthly mean surface runoff is the
        # average 3-hour accumulation over all 3-hour intervals in April 1979. To compute
        # monthly accumulation, use this formula:
        # Qs_acc (April){kg/m2} = Qs_acc (April){kg/m2/3hr} * 8{3hr/day} * 30{days}
        # so we need to devide by 8*30 to get m3 s^-1 per m2

        # to get m3 s^-1 per m2 we simply devide by the number of seconds in a 3hr interval
        dt = 60*60*3
        Q[runoff_var] ./= dt

        # replace missing values with 0
        Q[runoff_var] = coalesce.(Q[runoff_var], 0)
    end
    

    # apply linear reservoir model for subsruface runoff with a time delay factor of 45 days [Getirana et al. 2012]
    if haskey(Q, :Qsb_acc) # [~2 min]
        Tb = 45 # concentration time [days] - [Getirana et al. 2012]
        impulse_resonse = GGA.linear_reservoir_impulse_response_monthly(Tb)

        #p = lines(copy(Q[:Qsb_acc][180, 80, :])) # sanity check
        #xlims!(p.axis, [DateTime(2010, 1, 1), DateTime(2015, 1, 1)])
        GGA.apply_vector_impulse_resonse!(Q[:Qsb_acc], impulse_resonse)
        #lines!(dims(Q[:Qsb_acc][180, 80, :], :Ti).val, Q[:Qsb_acc][180, 80, :].data);p  # sanity check
    end

    # scale snow runoff for testing of impact on gmax
 

    # sum runoff variables
    Q = reduce(+, Q[k] for k in keys(Q))

    # deteremine runoff entering each basin 
    # [does not yet account for entering stream flow, that is the next step]
    begin #[~4.5 min for 2.3M basins]
        dcomid = Dim{:COMID}(basins.COMID)
        dTi = dims(Q, :Ti)

        # make time column major 
        river_inputs = zeros(eltype(Q), dTi, dcomid)

        @showprogress dt = 5 desc = "Adding surface + subsurface runoff to each basin ..." Threads.@threads for r in eachrow(basins) #[2 min]
            # multiply by area in m2 to get m3 s^-1 [raw area is in km2]
            river_inputs[:, At(r.COMID)] = Q[Near(r.centroid[1]), Near(r.centroid[2]), :] .* (r.unitarea .* 1e6)
        end
    end


    river_flux = GGA.flux_accumulate!(river_inputs, basins.COMID, basins.NextDownID, basins.HeadBasin, basins.basin02) #[2 min]

    # subset to glacier river reaches for saving 
    glacier_rivers = GeoDataFrames.read(glacier_routing_path)

    glacier_rivers = unique(vcat(glacier_rivers.RiverIDTrace...))
    glacier_river_land_flux = river_flux[:,At(glacier_rivers)]

    # save data as netcdf
    begin #[4 seconds]
        # data needs to be Float32 (not Real) for saving as netcdf
        data = Float32.(glacier_river_land_flux)

        if isfile(glacier_rivers_land_flux_paths[i])
            rm(glacier_rivers_land_flux_paths[i])
        end
        
        # save as netcdf
        NCDataset(glacier_rivers_land_flux_paths[i],"c") do ds
            data_dims = dims(data)

            # First all dimensions need to be defined
            for dim in data_dims
                dname = string(DimensionalData.name(dim));
                defDim(ds, dname, length(dim))
            end

            # now add the variables
            for dim in data_dims
                dname = string(DimensionalData.name(dim));
                defVar(ds, dname, val(val(dim)), (dname,))
            end

            v = defVar(ds, "flux", parent(data), string.(DimensionalData.name.(data_dims)))
            v.attrib["units"] = "m3 s^-1"
            v.attrib["long_name"] = "river flux"

            # add global attributes
            ds.attrib["title"] = "river flux from the average of 3 land surface models: CLSM, VIC, and NOAH"
        end
    end
end

# sanity check
validate_output = false
if validate_output # [5 min]
    # load in river discharge data and map to the glacier river reaches
    # get dimensions
    fn = river_Q_files[1]
    ds = Dataset(fn, "r")
    t = collect(ds[:time])

    #for fn in river_Q_files [m3 s^-1]
    dti = Dim{:Ti}(t)
    dcomid= Dim{:COMID}(basins.COMID)
    riverQ = zeros(Float32, (dti, dcomid))

    # using Treads causes crash
    # 5 min on single thread
    @showprogress dt = 5 desc = "Load river discharge from Collin's et al. 2024..." for fn in river_Q_files
    #fn = river_Q_files[1]
        
        Dataset(fn, "c") do ds
        #ds = Dataset(fn, "r")
        
            # need to read into memory first to make this efficient
            Qout = collect(ds["Qout"].var)
            rivids = collect(ds["rivid"].var)
            matching_rivids = intersect(dcomid, rivids)

            if !isempty(matching_rivids)
                for rid in matching_rivids
                    riverQ[:, At(rid)] = Qout[findfirst(rivids .== rid), :]
                end
            end
        end
    end

    # sanity check
    a = extrema(dims(riverQ, :Ti))
    b = extrema(dims(river_flux, :Ti))
    c = (max(a[1], b[1]), min(a[2], b[2]))

    Qindx = dims(riverQ, :Ti) .> c[1] .&& dims(riverQ, :Ti) .< c[2]
    fluxindx = dims(river_flux, :Ti) .> c[1] .&& dims(river_flux, :Ti) .< c[2]

    plot(collect(dropdims(mean(riverQ[Qindx, :], dims=:Ti), dims=:Ti)), collect(dropdims(mean(river_flux[fluxindx, :], dims=:Ti), dims=:Ti)))
end


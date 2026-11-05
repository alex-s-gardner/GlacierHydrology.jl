"""
    glacier_discharge(; datadir=pathlocal[:data_dir]) -> DataFrame

Load and combine glacier discharge data from multiple sources.

# Arguments
- `datadir`: Base directory containing glacier data files (default: pathlocal[:data_dir])

# Returns
- DataFrame containing glacier discharge information with columns:
  - `latitude`: Glacier center latitude
  - `longitude`: Glacier center longitude
  - `discharge_gtyr`: Discharge rate in gigatons per year
  - `discharge_err_gtyr`: Uncertainty in discharge rate
  - `frontal_ablation_gtyr`: Frontal ablation rate in gigatons per year

# Description
Combines discharge data from Kochtitzky (2022) for the Northern Hemisphere and 
Fuerst (2023) for Patagonia. Matches Patagonian glaciers with RGI database entries
to obtain coordinates, and standardizes data format across sources.
"""
function glacier_discharge()
    
    # Kochtitzky NH discharge and terminus retreate 
    nothern_hemisphere = CSV.read(pathlocal[:discharge_nh], DataFrame; header=14, skipto=16)

    df = DataFrame(Shapefile.Table(pathlocal[:rgi6_southern_andes]))

    patagonia = CSV.read(pathlocal[:discharge_npi], DataFrame; header=25, skipto=27)
    patagonia = vcat(patagonia, CSV.read(pathlocal[:discharge_spi], DataFrame; header=25, skipto=27))

    # find matching glaciers in RGI 
    patagonia[!, :"Latitude"] .= NaN
    patagonia[!, :"Longitude"] .= NaN
    df.Name[ismissing.(df.Name)] .= "junk"
    patagonia.Names = replace.(patagonia.Names, "Témpano" => "Tempano")
    patagonia.Names[Base.contains.(patagonia.Names, "Grey")] .= "Grey + Dickson"
    patagonia.Names[Base.contains.(patagonia.Names, "Upsala")] .= "Upsala + Cono"
    patagonia.Names[Base.contains.(patagonia.Names, "O'Higgins")] .= "OHiggins"

    for r in eachrow(patagonia)
        index = findfirst((r.Names .== df.Name) .| (r.Names .== df.RGIId))
        if isnothing(index)
            println("cold not find match for: $(r.Names)")

            continue
        end
        r.Latitude = df.CenLat[index]
        r.Longitude = df.CenLon[index]
    end

    df = DataFrame()
    df[!, :latitude] = vcat(nothern_hemisphere.lat, patagonia.Latitude)
    df[!, :longitude] = vcat(nothern_hemisphere.lon, patagonia.Longitude)
    df[!, :discharge_gtyr] = vcat(nothern_hemisphere[:, "2010_2020_mean_discharge_gt_per_year"], (patagonia[:, "Calving FG"] .* δice / 1000))
    df[!, :discharge_err_gtyr] = vcat(nothern_hemisphere[:, "2010_2020_mean_flux_err_gt"], (patagonia[:, "Unc. Calving FG"] .* δice / 1000))
    df[!, :frontal_ablation_gtyr] = vcat(nothern_hemisphere[:, "Frontal_ablation_2010_to_2020_gt_per_yr_mean"], (patagonia[:, "Frontal ablation Minowa (2000-2019)"] .* δice / 1000))
    df[!, :frontal_ablation_gtyr] = vcat(nothern_hemisphere[:, "Frontal_ablation_2010_to_2020_gt_per_yr_mean"], (patagonia[:, "Unc. Frontal ablation Minowa (2000-2019)"] .* δice / 1000))

    return df
end


"""
    read_grace_rgi(; path2file=setpaths()[:grace_rgi])

Read GRACE (Gravity Recovery and Climate Experiment) data for RGI (Randolph Glacier Inventory) regions.

This function loads GRACE data from a MATLAB file and converts region codes from letter-based 
abbreviations to standard RGI numerical identifiers (e.g., "ALA" to "rgi1").

# Arguments
- `path2file`: Path to the GRACE RGI data file (MATLAB format)

# Returns
- A dictionary mapping RGI region identifiers to their corresponding GRACE data

# Note
The function converts region codes from the original letter-based abbreviations (e.g., "ALA" for Alaska)
to the standard RGI numerical format (e.g., "rgi1").
"""
function read_grace_rgi(;path2file=setpaths()[:grace_rgi])

    grace0 = matread(path2file)

    # rgi leter to digit mapping
    old2new = Dict(
        "PAT" => "rgi18", 
        "GRE" => "rgi5", 
        "NAS" => "rgi10", 
        "ICE" => "rgi6", 
        "TRP" => "rgi16",
        "SAW" => "rgi14",
        "SAE" => "rgi15",
        "CDN" => "rgi3",
        "CAS" => "rgi13",
        "AND" => "rgi17",
        "CDS" => "rgi4",
        "ANT" => "rgi19",
        "CEU" => "rgi11",
        "ALA" => "rgi1",
        "SVB" => "rgi7",
        "WNA" => "rgi2", 
        "NEZ" => "rgi18",
        "RAI" => "rgi9",
        "SCA" => "rgi8"
    )


    grace = Dict(any(keys(old2new).== key) ? (old2new[key]) => val : (key) => val for (key, val) in grace0)
    return grace
end

"""
    grace_masschange(; path2file=setpaths()[:grace_rgi])

Load and organize GRACE (Gravity Recovery and Climate Experiment) mass change data for glacier regions.

This function reads GRACE data for RGI (Randolph Glacier Inventory) regions and organizes it into a
dimensional array structure with dimensions for region, date, and error bounds.

# Arguments
- `path2file`: Path to the GRACE RGI data file (MATLAB format). Defaults to the path from `setpaths()[:grace_rgi]`.

# Returns
- A DimensionalData array with dimensions for RGI regions (1-19, plus 98 for HMA and 99 for Global),
  dates, and error bounds (false=value, true=error).

# Note
The error values are multiplied by 2 to represent 95% confidence intervals.
"""
function grace_masschange(; path2file=setpaths()[:grace_rgi])
    grace_raw = read_grace_rgi(; path2file)

    # organize GRACE grace data into a DimArray
    drgi = Dim{:rgi}(collect([1:19..., 98:99...]))
    ddate = Dim{:date}(vec(datenum2date.(grace_raw["rgi1"]["dM_gt_mdl_fill_date"])))
    derror = Dim{:error}([false, true])

    grace = fill(NaN, drgi, ddate, derror)

    for rgi in drgi
        #rgi = drgi[1]
        if rgi == 98
            rgi_id = "HMA"
        elseif rgi == 99
            rgi_id = "Global"
        else
            rgi_id = "rgi$(rgi)"
        end

        if haskey(grace_raw, rgi_id)
            grace[At(rgi), :, At(false)] = vec(grace_raw[rgi_id]["dM_gt_mdl_fill"])

            ## multiply by 2 to get 95% confidence interval
            grace[At(rgi), :, At(true)] = vec(grace_raw[rgi_id]["dM_sigma_gt_mdl_fill"]) * 2
        end
    end

    return grace
end

"""
    glambie2024(; path2glambie=paths = setpaths()[:glambie_2024]) -> DimArray

Load and process the GlaMBIE 2024 glacier mass balance dataset.

# Keywords
- `path2glambie`: Path to the GlaMBIE 2024 CSV file (default: from setpaths())

# Returns
- A DimensionalData.DimArray containing cumulative glacier mass balance in Gt
  with dimensions for RGI region (1-19, 98, 99), date (2000-2024), and error flag

# Notes
- Dates are shifted by up to 6 months from original data to align across hemispheres
- Missing values are linearly interpolated
- Region 98 (HMA) is calculated as the sum of regions 13-15
"""
function glambie2024(; path2glambie=paths = setpaths()[:glambie_2024])
    df = CSV.read(path2glambie, DataFrame)

    drgi = Dim{:rgi}(vcat(collect(1:19), 98, 99))
    derror = Dim{:error}([false, true])
    @warn "GlaMBIE dates are shifted by up to 6 months from the original data to align dates across hemispheres"
    ddate = Dim{:date}(collect(decimalyear2datetime.(2000:2024)); name="GlaMBIE 2024 cumulative glacier mass balance [GLAMBIE]")

    glambie_Gt = fill(NaN, drgi, ddate, derror)

    for (i, rgi) in enumerate(drgi)
        if rgi == 98
            continue
        elseif rgi == 99
            i = i - 1
        end
        scol = (i - 1) * 3 + 1
        date1 = round.(df[:, scol])
        valid = .!ismissing.(date1)
        date1 = collect(decimalyear2datetime.(date1[valid]))

        glambie_Gt[At(rgi), At(date1), At(false)] = coalesce.(df[valid, scol+1], NaN)
        glambie_Gt[At(rgi), At(date1), At(true)] = coalesce.(df[valid, scol+2], NaN)
    end


    glambie_Gt[At(98), :, :] .= 0
    for rgi in 13:15
        glambie_Gt[At(98), :, At(false)] .+= glambie_Gt[At(rgi), :, At(false)]
        glambie_Gt[At(98), :, At(true)] .+= glambie_Gt[At(rgi), :, At(true)] .^ 2
    end

    glambie_Gt[At(98), :, At(true)] = sqrt.(glambie_Gt[At(98), :, At(true)])

    @warn "GlaMBIE data is linearly interpolated over missing date values"
    for rgi in drgi
        for err_flag in derror
            valid = .!isnan.(glambie_Gt[At(rgi), :, At(err_flag)])
            if !all(valid)
                fit = curve_fit(model1_trend, collect(decimalyear.(ddate[valid])), collect(glambie_Gt[At(rgi), :, At(err_flag)][valid]), p_trend)
                glambie_Gt[At(rgi), .!valid, At(err_flag)] = model1_trend(decimalyear.(collect(ddate[.!valid])), fit.param)
            end
        end
    end

    return (glambie_Gt)
end

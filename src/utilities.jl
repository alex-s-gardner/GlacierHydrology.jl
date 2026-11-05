

const MATLAB_EPOCH = Dates.DateTime(-1, 12, 31)

"""
    datenum2date(n) -> DateTime

Convert MATLAB datenum format to Julia DateTime.

# Arguments
- `n`: MATLAB datenum value (days since 0000-01-00)

# Returns
- `DateTime`: Equivalent Julia DateTime object
"""
function datenum2date(n) 
    d = MATLAB_EPOCH + Dates.Millisecond(round(Int64, n * 1000 * 60 * 60 * 24))
    return d
end

rgi2label = Dict(
    "rgi1" => "Alaska",
    "rgi2" => "Western Canada",
    "rgi3" => "Arctic Canada N.",
    "rgi4" => "Arctic Canada S.",
    "rgi5" => "Greenland",
    "rgi6" => "Iceland",
    "rgi7" => "Svalbard",
    "rgi8" => "Scandinavia",
    "rgi9" => "Russian Arctic",
    "rgi10" => "North Asia",
    "rgi11" => "Central Europe",
    "rgi12" => "Caucasus",
    "rgi13" => "Central Asia",
    "rgi14" => "South Asia W.",
    "rgi15" => "South Asia E.",
    "rgi16" => "Tropics",
    "rgi17" => "Southern Andes",
    "rgi18" => "New Zealand",
    "rgi19" => "Antarctic",
    "hma" => "High Mtn. Asia",
    "global" => "Global",
    "global_ep" => "Global [excluding periphery]",
    "ghll" => "All Missions: rgi 2, 10, 11, 12, 13, 14, 15, 17, 18",
    "hll" => "ASTER, ICESat/-2: rgi 1, 3, 4, 5, 6, 7, ,8, 9, 19",
    "hll_ep" => "ASTER, ICESat/-2 no periphery: rgi 1, 3, 4, 6, 7, ,8, 9",
)


"""
    rginum2label(rginum)

Converts an RGI (Randolph Glacier Inventory) region number to its corresponding label.
Takes a numeric RGI region identifier and returns the human-readable label by first
converting the number to a text key using `rginum2txt` and then looking up the label
in the `rgi2label` dictionary.
"""
function rginum2label(rginum) 
    return getindex(rgi2label,getindex(rginum2txt, rginum))
end

rginum2txt = Dict(
    1 => "rgi1" ,
    2 => "rgi2" ,
    3 => "rgi3" ,
    4 => "rgi4" ,
    5 => "rgi5",
    6 => "rgi6" ,
    7 => "rgi7" ,
    8 => "rgi8" ,
    9 => "rgi9" ,
    10 => "rgi10" ,
    11 => "rgi11" ,
    12 => "rgi12" ,
    13 => "rgi13" ,
    14 => "rgi14" ,
    15 => "rgi15" ,
    16 => "rgi16" ,
    17 => "rgi17" ,
    18 => "rgi18" ,
    19 => "rgi19" ,
    98 => "hma" ,
    99 => "global" ,
)

rginum2enclosed_alphanumerics = Dict(
    1 => "①" ,
    2 => "②" ,
    3 => "③" ,
    4 => "④" ,
    5 => "⑤",
    6 => "⑥" ,
    7 => "⑦" ,
    8 => "⑧" ,
    9 => "⑨" ,
    10 => "⑩" ,
    11 => "⑪" ,
    12 => "⑫" ,
    13 => "⑬" ,
    14 => "⑭" ,
    15 => "⑮" ,
    16 => "⑯" ,
    17 => "⑰" ,
    18 => "⑱" ,
    19 => "⑲" ,
    98 => "Ⓗ",
    99 => "Ⓖ",
)


"""
    df_tsfit!(df, tsvars; progress=true, datelimits=nothing) -> DataFrame

Fit a least squares model to time series variables in a DataFrame, calculating offset, trend, amplitude, and phase.

# Arguments
- `df`: DataFrame containing time series data
- `tsvars`: Column names of time series variables to fit

# Keywords
- `progress=true`: Whether to display a progress bar
- `datelimits=nothing`: Optional tuple of (start_date, end_date) to limit the date range for fitting

# Returns
- Modified DataFrame with added columns for offset, trend, amplitude, and phase for each time series variable
"""
function df_tsfit!(df, tsvars; progress=true, datelimits = nothing)
    prog = progress ? ProgressMeter.Progress(length(tsvars); desc="Fitting LSQ model to timeseries...") : nothing

    for tsvar in tsvars

        ddate = colmetadata(df, tsvar, "date")
        x = decimalyear.(ddate)
        x = x .- ceil(mean(x))

        out_var_offset = string(tsvar)*"_offset"
        out_var_trend = string(tsvar)*"_trend"
        out_var_amp = string(tsvar)*"_amplitude"
        out_var_phase = string(tsvar)*"_phase"
        df[!, out_var_offset] = zeros(nrow(df))

        df[!, out_var_trend] = zeros(nrow(df))
        df[!, out_var_amp] = zeros(nrow(df))
        df[!, out_var_phase] = zeros(nrow(df))

        if !isnothing(datelimits)
            date_index = (ddate .> datelimits[1]) .& (ddate .< datelimits[2])
        else
            date_index = trues(length(ddate))
        end

        Threads.@threads for g in eachrow(df)
        
            y = g[tsvar];
            valid = .!isnan.(y) .& date_index

            fit = curve_fit(
                offset_trend_seasonal2, 
                x[valid], 
                y[valid], 
                p_offset_trend_seasonal
                )

            g[out_var_offset] = fit.param[1]
            g[out_var_trend] = fit.param[2]
            g[out_var_amp] = hypot(fit.param[3], fit.param[4])
            g[out_var_phase] = 365.25 * seasonal_peak_fraction(fit.param[3], fit.param[4])
        end

        # Update the progress meter
        progress && ProgressMeter.update!(prog, prog.counter .+ 1)
    end
    return df
end


"""
    connected_groups(connectivity) -> Vector{Int64}

Group elements into connected components based on their connectivity relationships.

# Arguments
- `connectivity`: Vector of vectors where each element contains indices of connected elements

# Returns
- Vector of integers assigning each element to a group ID

# Details
Identifies connected components in a graph represented by the connectivity list.
Elements with empty connectivity lists are assigned to their own individual groups.
"""
function connected_groups(connectivity)

    # Track which tiles have been assigned to groups
    assigned = falses(length(connectivity))
    group_id = zeros(Int64, length(connectivity))

    # Initialize group counter
    group_id0 = 1

    while !all(assigned)
        # Find first unassigned tile with intersections
        start_idx = findfirst(.!assigned)

        if isempty(connectivity[start_idx])
            assigned[start_idx] = true
            group_id[start_idx] = group_id0
            group_id0 += 1
            continue
        end

        # Build group starting from this tile
        current_group = Int64[start_idx]
        assigned[start_idx] = true
        group_id[start_idx] = group_id0

        # Keep expanding group until no more intersections found
        expanded = true
        while expanded
            expanded = false

            # Check all tiles in current group
            for tile_idx in current_group
                # Look at all intersecting tiles
                for intersect_idx in connectivity[tile_idx]
                    if !assigned[intersect_idx]
                        push!(current_group, intersect_idx)
                        assigned[intersect_idx] = true
                        group_id[intersect_idx] = group_id0
                        expanded = true
                    end
                end
            end
        end

        group_id0 += 1
    end

    return group_id
end


"""
    nc2dd(ds, varname)

Convert a NetCDF variable to a DimensionalData.DimArray.

# Arguments
- `ds`: NCDataset containing the variable
- `varname`: String name of the variable to convert

# Returns
- DimArray with dimensions and data from the NetCDF variable

# Notes
- Preserves dimension names and coordinates from the NetCDF dataset
- Converts dimension names to symbols for DimensionalData compatibility
"""
function nc2dd(cfv)
    dnames = NCDatasets.dimnames(cfv)

    if in("units", keys(cfv.attrib)) 
        unit = cfv.attrib["units"]
        try
            units = Unitful.uparse(unit; unit_context=Unitful.unitmodules)
        catch
            units = unit_convert[unit]
        end
    else
        units = NoUnits
    end

    da = DimArray(collect(cfv[:,:]) * units, tuple([Dim{Symbol(dname)}(cfv[dname]) for dname in dnames]...))
    return da
end

"""
    check_for_all_nans(d) -> Nothing

Check if any key in a dictionary contains only NaN values.

# Arguments
- `d`: Dictionary to check

# Returns
- Nothing

# Throws
- Error if any key in the dictionary contains only NaN values
"""
function check_for_all_nans(d)
    for k in keys(d)
        if all(isnan.(d[k]))
            error("all NaNs for $k")
        end
    end
end


"""
    resample(colorscheme::Symbol, n) -> ColorSchemes.ColorScheme

Create a new color scheme by resampling an existing one to a specified number of colors.

# Arguments
- `colorscheme::Symbol`: Name of the source color scheme
- `n`: Number of colors in the new scheme

# Returns
- A new ColorScheme with n evenly spaced colors from the original scheme
"""
function resample(colorscheme::Symbol, n)
    newscheme = ColorSchemes.ColorScheme(
        get(ColorSchemes.colorschemes[colorscheme], (0:n-1) / (n-1))
    )
    return newscheme
end


"""
    dimarray2netcdf(da, filename; name="var", units=nothing, global_attributes=nothing) -> String

Convert a DimensionalData.DimArray -or- DimensionalData.DimStack to a NetCDF file.

# Arguments
- `da`: DimensionalData.DimArray to convert
- `filename`: Path where the NetCDF file will be saved

# Keywords
- `name="var"`: Name for the main variable in the NetCDF file
- `units=nothing`: Units for the main variable (uses DimArray units if not specified)
- `global_attributes=nothing`: Dictionary of global attributes to add to the NetCDF file

# Returns
- Path to the created NetCDF file
"""
function dimarray2netcdf(da::DimArray, filename; name = "var", units = nothing, global_attributes = nothing)
    
    if isfile(filename)
        rm(filename)
    end
    
    NCDataset(filename, "c") do ds
        
        # get dimensions
        da_dims = dims(da)

        # step 1: define dimensions
        for dim in da_dims
            dname = string(DimensionalData.name(dim));
            defDim(ds, dname, length(dim))
        end

        # step 2: add dim variables & metadata
        for dim in da_dims
            dname = string(DimensionalData.name(dim));
            d = defVar(ds, dname, val(val(dim)), (dname,))

            if eltype(dim) <: Number
                d.attrib["units"] = string(Unitful.unit(dim[1]))
            end

            # add metadata
            for (k,v) in DD.metadata(dim)
                d.attrib[k] = v
            end
        end

        # step 3: add variable
        v = defVar(ds, name, ustrip.(parent(da)), string.(DimensionalData.name.(da_dims)))
        
        # step 4: add variable metadata
        if isnothing(units)
            if eltype(da) <: Number
                v.attrib["units"] = string(Unitful.unit(da[1]))
            end
        else
            v.attrib["units"] = string(units)
        end

        for (k,v) in DD.metadata(da)
            d.attrib[k] = v
        end

        # step 5: add global attributes
        if !isnothing(global_attributes)
            for (k,v) in global_attributes
                ds.attrib[k] = v
            end
        end

        return filename
    end;
end



"""
    dimstack2netcdf(ds::DimStack, filename; kwargs...) -> String

Write a `DimStack` to a NetCDF file, preserving all dimensions, data, and metadata.

# Arguments
- `ds::DimStack`: The DimStack to be saved. Each entry should be a `DimArray` with associated dimensions and metadata.
- `filename`: Path (as a `String`) to the output NetCDF file.

# Keywords
- `kwargs...`: Additional keyword arguments passed to `defVar` for each variable (e.g., `attrib`, `fillvalue`).

# Returns
- `String`: The filename of the written NetCDF file.

# Details
- All dimensions in the `DimStack` are defined as NetCDF dimensions.
- Each dimension variable is written, including its metadata.
- All variables in the `DimStack` are written to the file, with their metadata and any optional fill values or attributes.
- Global attributes from the `DimStack` metadata are added to the NetCDF file.
- The function overwrites any existing file at `filename`.
"""
function dimstack2netcdf(ds::DimStack, filename; kwargs...)

    NCDataset(filename, "c"; format=:netcdf4) do nc
        # get dimensions
        ds_dims = dims(ds)

        # step 1: define dimensions
        for dim in ds_dims
            dname = string(DimensionalData.name(dim))
            defDim(nc, dname, length(dim))
        end

        # step 2: add dim variables & metadata
        for dim in ds_dims
            dname = string(DimensionalData.name(dim))
            dim_val = val(val(dim))
            if eltype(dim_val) == Date
                dim_val = DateTime.(dim_val)
            end

            d = defVar(nc, dname, dim_val, (dname,))
       
            # add dimension attributes from dimension metadata
            for (k, v) in DD.metadata(dim)
                d.attrib[k] = v
            end
        end

        # step 3: add variable
        varnames = keys(ds)
        for name in varnames

            var_meta = DD.metadata(ds[name]);
            v = defVar(nc, string(name), ustrip.(parent(ds[name])), string.(DimensionalData.name.(dims(ds[name]))); kwargs...)
     
            # step 4: add variable attributes from variable metadata
            for (k, var0) in var_meta
                v.attrib[k] = var0
            end
        end

        # step 5: add global attributes from global metadata
        for (k, val) in DD.metadata(ds)
            nc.attrib[string(k)] = val
        end
    end

    return filename
end


"""
    netcdf2dimstack(filename) -> DimStack

Load a NetCDF file and convert it to a `DimStack`, preserving dimensions, data, and metadata.

# Arguments
- `filename`: Path to the NetCDF file.

# Returns
- `DimStack`: A DimStack containing all non-dimension variables from the NetCDF file, with associated dimensions and metadata.
"""
function netcdf2dimstack(filename)
    ds0 = NCDataset(filename)
    ds = ncdataset2dimstack(ds0)
    return ds
end


"""
    ncdataset2dimstack(ds0) -> DimStack

Convert an open NCDataset (from NCDatasets.jl) to a DimStack, preserving dimensions, data, and metadata.

# Arguments
- `ds0`: An open NCDataset object.

# Returns
- `DimStack`: A DimStack containing all non-dimension variables from the dataset, with associated dimensions and metadata.

# Notes
- Dimension variables are detected and used to construct Dim objects.
- All global and variable-level attributes are preserved in the resulting DimStack and DimArrays.
"""
function ncdataset2dimstack(ds0)
    # get the dimension names
    dnames = NCDatasets.dimnames(ds0)

    dim0 = Dict()
    for dim in dnames
        if haskey(ds0, dim)
            dim0[dim] = Dim{Symbol(dim)}(ds0[dim])
        else
            dim0[dim] = Dim{Symbol(dim)}(1:ds0.dim[dim])
        end
    end

    vnames = tuple(setdiff(keys(ds0), dnames)...)
    ds = DimStack((; zip(Symbol.(vnames), [DimArray(ds0[vname], getindex.(Ref(dim0), dimnames(ds0[vname])); metadata=Dict(zip(keys(ds0[vname].attrib), getindex.(Ref(ds0[vname].attrib), keys(ds0[vname].attrib))))) for vname in vnames])...); metadata=Dict(zip(keys(ds0.attrib), getindex.(Ref(ds0.attrib), keys(ds0.attrib)))))
    
    return ds
end


"""
    netcdf2dimarray(filename; varname=nothing) -> DimensionalData.DimArray

Convert a NetCDF file to a DimensionalData.DimArray.

# Arguments
- `filename`: Path to the NetCDF file

# Keywords
- `varname=nothing`: Name of the variable to convert. If not specified, attempts to
  automatically detect the single non-dimension variable in the file.

# Returns
- DimArray with dimensions and data from the NetCDF variable, preserving units and metadata

# Throws
- Error if varname is not specified and there isn't exactly one variable in the file
"""
function netcdf2dimarray(filename; varname=nothing)
    # read the netcdf file
    ds = NCDatasets.Dataset(filename) 

    # get the dimension names
    dnames = NCDatasets.dimnames(ds)

    # get the variable name
    if varname === nothing
        varname = setdiff(keys(ds), dnames)
        if length(varname) != 1
            error("varname must be specified if there is not exactly one variable in the netcdf file")
        end
        varname = first(varname)
    end

    # get the units of variable
    da = _var_from_netcdf(ds, varname, dnames)
    
    return da
end


"""
    _var_from_netcdf(ds, varname, dnames) -> DimensionalData.DimArray

Create a DimArray from a NetCDF variable.

# Arguments
- `ds`: NCDataset containing the variable
- `varname`: Name of the variable to convert
- `dnames`: Names of dimensions in the dataset

# Returns
- DimArray with appropriate dimensions and units (if available)
"""
function _var_from_netcdf(ds, varname, dnames)
    units = _units_from_netcdf(ds, varname)
    if !isnothing(units)
        da = DimArray(ds[varname] * units, tuple(_dim_from_netcdf.(Ref(ds), dnames)...); metadata=ds[varname].attrib)
    else
        da = DimArray(ds[varname], tuple(_dim_from_netcdf.(Ref(ds), dnames)...); metadata=ds[varname].attrib)
    end
    return da
end

"""
    _dim_from_netcdf(ds, dname) -> DimensionalData.Dimension

Create a dimension from a NetCDF dimension variable.

# Arguments
- `ds`: NCDataset containing the dimension
- `dname`: Name of the dimension to convert

# Returns
- Dimension with values and units (if available) and original metadata
"""
function _dim_from_netcdf(ds, dname)
    units = _units_from_netcdf(ds, dname)
    if !isnothing(units)
        dim = Dim{Symbol(dname)}(ds[dname][:]*units; metadata=ds[dname].attrib)
    else
        dim = Dim{Symbol(dname)}(ds[dname][:]; metadata=ds[dname].attrib)
    end
    return dim
end

"""
    _units_from_netcdf(ds, varname) -> Union{Unitful.Units, Nothing}

Extract and parse units from a NetCDF variable.

# Arguments
- `ds`: NCDataset containing the variable
- `varname`: Name of the variable to extract units from

# Returns
- Parsed Unitful units if available and parsable
- `NoUnits` if units attribute exists but is empty
- `nothing` if units are not available or cannot be parsed
"""
function _units_from_netcdf(ds, varname)

    # units can not be applied to non-numeric data
    if !(eltype(ds[varname]) <: Number)
        units = nothing
    else
        if haskey(ds[varname].attrib, "units")
            if isempty(ds[varname].attrib["units"])
                units = NoUnits
            else
                units = ds[varname].attrib["units"];
                if isempty(units)
                    units = nothing
                else
                    try
                        units = uparse(units; unit_context=Unitful.unitmodules)
                    catch
                        @warn "could not parse units for $varname"
                        units = nothing
                    end
                end
            end
        else 
            units = nothing
        end
    end
    return units
end

unit_convert = Dict()
unit_convert["m3 s^-1"] = u"m^3/s"
unit_convert["kg s^-1"] = u"kg/s"


"""
    _geotile_load_align(; surface_mask="glacier", geotile_width=2, geotile_order=nothing)

Load and align geotiles for a given surface mask.

# Arguments
- `surface_mask`: Surface mask type (default: "glacier")
- `geotile_width`: Width of geotiles in degrees (default: 2)
- `geotile_order`: Optional geotile order for alignment (default: nothing)

# Returns
Tuple of (geotiles0, reg) where geotiles0 is the processed geotiles DataFrame and reg contains region information
"""
function _geotile_load_align(; surface_mask="glacier", geotile_width=2, geotile_order=nothing, only_geotiles_w_area_gt_0=true)

    # Generate geotiles for this mask
    geotiles0 = geotiles_mask_hyps(surface_mask, geotile_width)

    # Make geotiles mutually exclusive by RGI region
    geotiles0, reg = geotiles_mutually_exclusive_rgi!(copy(geotiles0))

    
    #rename!(geotiles0,"$(surface_mask)_area_km2" => "area_km2")

    if !isnothing(geotile_order)
        # Align geotile indices with elevation change dimensions
        gt_index, dgt_index = intersectindices(geotiles0.id, geotile_order)

        geotiles0 = geotiles0[gt_index, :]
    end

    # Rename area column to include surface mask
    rename!(geotiles0, "$(surface_mask)_area_km2" => "area_km2")

    
    if only_geotiles_w_area_gt_0
        geotiles0 = geotiles0[sum.(geotiles0.area_km2) .> 0, :]
    end

    return geotiles0
end

"""
    ts_seasonal_model(ts, interval)

Fit a simple seasonal model to a time series and extract the amplitude and phase of the seasonal cycle.

# Arguments
- `ts`: A DimArray or similar time series object with a :date or :time dimension.
- `interval`: A tuple of two DateTime objects specifying the interval over which to fit the model.

# Returns
- `amplitude`: The amplitude of the seasonal (annual) cycle.
- `phase_peak_month`: The month (as an integer) when the seasonal peak occurs.

# Description
Fits a model of the form:
    y = a + b*t + c*cos(2πt) + d*sin(2πt)
to the time series over the specified interval, where t is in decimal years (centered).
Returns the amplitude and the month of the seasonal peak.
"""
function ts_seasonal_model(ts; interval=nothing)
    if hasdim(ts, :date)
        ddate = dims(ts, :date)
    elseif hasdim(ts, :time)
        ddate = dims(ts, :time)
    else
        error("ts has no :date or :Ti dimension")
    end

    if isnothing(interval)
        interval = (minimum(ddate), maximum(ddate))
    end

    decyear = decimalyear.(parent(values(ddate[interval[1]..interval[2]])))
    decyear = decyear .- round(mean(decyear))

    # construct design matrix
    n = length(decyear)
    M = hcat(ones(n), decyear, cos.(2π .* decyear), sin.(2π .* decyear))

    # fit seasonal model
    p = M \ ts[interval[1]..interval[2]]

    amplitude = sqrt(p[3]^2 + p[4]^2)

    phase_peak = seasonal_peak_fraction(p[3], p[4])
    phase_peak_month = month(decimalyear2datetime(phase_peak))

    trend = p[2]
    return (;trend, amplitude, phase_peak_month)
end


"""
    calculate_slope(data::DimArray)

Calculate the linear slope of a line from a DimensionalData vector with time dimension.

# Arguments
- `data`: DimArray with time dimension (Ti)

# Returns
- `slope`: The calculated slope value

# Examples
```julia
using DimensionalData

# Create time dimension
ti = Ti(1:10)
data = DimArray(collect(1:10), ti)

# Calculate slope
slope = calculate_slope(data)
```
"""
function calculate_slope(data::DimArray)
    # Extract time dimension and data
    ti_dim = dims(data, :Ti)
    time_values = decimalyear.(collect(val(ti_dim)))
    data_values = collect(data)

    # Remove any NaN values
    valid_mask = .!isnan.(data_values)
    if sum(valid_mask) < 2
        error("Need at least 2 valid data points to calculate slope")
    end

    time_clean = time_values[valid_mask]
    data_clean = data_values[valid_mask]

    # Calculate slope using least squares formula
    n = length(time_clean)
    x̄ = mean(time_clean)
    ȳ = mean(data_clean)

    numerator = sum((time_clean .- x̄) .* (data_clean .- ȳ))
    denominator = sum((time_clean .- x̄) .^ 2)

    if denominator == 0
        error("Cannot calculate slope: all time values are identical")
    end

    slope = numerator / denominator
    return slope
end

function add_single_rgi_column!(df)

    rgis = parse.(Int8, replace.(names(df)[occursin.("rgi", names(df))], "rgi" => ""))
    
    df[!, :rgi] = zeros(Int8, nrow(df))

    for rgi in rgis
        col = "rgi$rgi"
        df[df[:, col].>0, :rgi] .= rgi
    end

    return df
end

"""
    scale2linear(scaled)

Convert a scaled (asymmetric) value back to its corresponding linear value, as used in the scale2dist transformation.

# Arguments
- `scaled::Real`: Scaled value output from `scale2dist`.

# Returns
- `linear::Real`: The corresponding linear value.

# Example
```julia
scale2linear(-2.0)  # returns 0.5
scale2linear(2.0)   # returns 2.0
scale2linear(1.0)   # returns 1.0
```
"""
function scale2linear(scaled)
    if scaled < 1
        linear = -1 / scaled
    else
        linear = scaled
    end

    if linear < -1
        linear += 1
    else
        linear -= 1
    end

    return linear
end

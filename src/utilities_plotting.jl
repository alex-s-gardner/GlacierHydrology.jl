

function plot_model_fit(gemb, discharge, dv_altim, cost, geotiles_in_group; colormap=:thermal)


    dv_altim0 = dropdims(sum(dv_altim[geotile = At(geotiles_in_group)]; dims = :geotile); dims =:geotile)
    gemb0 = dropdims(sum(gemb[geotile=At(geotiles_in_group)], dims=:geotile), dims=:geotile)
    discharge0 = sum(discharge[:discharge][geotile=At(geotiles_in_group)])

    # align dv_altim and dv_gemb
    ddate_gemb = dims(gemb0, :date)
    ex = extrema(dims(dv_altim, :date))
    index_gemb = (ddate_gemb .>= ex[1]) .& (ddate_gemb .<= ex[2])
    gemb0 = gemb0[date=At(ddate_gemb[index_gemb].val)]

    discharge0 = ((decimalyear.(dims(gemb0, :date)) .- decimalyear.(dims(gemb0, :date))[1]) .+ (decimalyear.(dims(gemb0, :date))[2] .- decimalyear.(dims(gemb0, :date))[1]))  .* discharge0

    index_minimum = argmin(cost)
    pscale_best = DimPoints(cost)[index_minimum][1]
    ΔT_best = DimPoints(cost)[index_minimum][2]

    best_smb = gemb_dv_sample(pscale_best, ΔT_best, gemb0[:smb]);
    best_fac = gemb_dv_sample(pscale_best, ΔT_best, gemb0[:fac]);
    best_fit = gemb_dv_sample(pscale_best, ΔT_best, gemb0[:dv]);

    
    fit_gemb  = ts_seasonal_model(best_fit; interval=nothing)
    fit_altim = ts_seasonal_model(dv_altim0; interval=nothing)
    fit_gemb = fit_gemb.amplitude
    fit_altim = fit_altim.amplitude
    fit_diff = round(Int, ((fit_gemb - fit_altim) / fit_altim)*100)
    fit_gemb = round(fit_gemb, digits=2)
    fit_altim = round(fit_altim, digits=2)


    println("dv amplitude: observed: $(fit_altim) km³, modeled: $(fit_gemb) km³: difference: $(fit_diff)%")

    f = plot_model_fit(best_smb, best_fac, discharge0, best_fit, dv_altim0, cost; colormap)
    return f
end



"""
    plot_model_fit(best_smb, best_fac, best_discharge, best_fit, dv0, cost_metric; xlims, xtickspacing, colormap=:thermal)

Create a two-panel figure visualizing the model fit and cost metric for glacier mass balance modeling.

# Arguments
- `best_smb`: Array of best-fit Surface Mass Balance (SMB) values over time.
- `best_fac`: Array of best-fit Firn Air Content (FAC) values over time.
- `best_discharge`: Array of best-fit discharge values over time.
- `best_fit`: Array of best-fit modeled height anomaly values over time.
- `dv0`: Array of observed height anomaly values over time.
- `cost_metric`: 2D array of cost metric values (e.g., RMSE or MAD) as a function of model parameters.
- `xlims`: Tuple specifying the x-axis limits for the time series plot.
- `xtickspacing`: Spacing between x-axis ticks for the time series plot.
- `colormap`: (Optional) Colormap to use for plotting (default: `:thermal`).

# Returns
- `f`: A CairoMakie Figure object with two subplots:
    - Left: Time series of SMB, FAC, discharge, observed, and modeled height anomalies.
    - Right: Contour plot of the cost metric as a function of model parameters.

"""
function plot_model_fit(best_smb, best_fac, best_discharge, best_fit, dv0, cost_metric; colormap=:thermal)
    f = _publication_figure(columns=2, rows=1)
    ax1 = CairoMakie.Axis(f[1, 1])
    ax2 = CairoMakie.Axis(f[1, 2])


    plot_best_fit!(ax1, best_smb, best_fac, best_discharge, best_fit, dv0; colormap)
    ax2, crange = plot_cost_metric!(ax2, cost_metric; colormap)
    Colorbar(f[1, 3]; limits=crange, colormap)#label="cost"

    return f
end

"""
    plot_best_fit!(ax, best_smb, best_fac, best_discharge, best_fit, dv0; xlims, xtickspacing, colormap=:thermal)

Plot the best-fit glacier mass balance model components and observations on a given axis.

# Arguments
- `ax`: A CairoMakie Axis object to plot on.
- `best_smb`: Array of best-fit Surface Mass Balance (SMB) values over time.
- `best_fac`: Array of best-fit Firn Air Content (FAC) values over time.
- `best_discharge`: Array of best-fit discharge values over time.
- `best_fit`: Array of best-fit modeled height anomaly values over time.
- `dv0`: Array of observed height anomaly values over time.

# Keyword Arguments
- `xlims`: Tuple specifying the x-axis limits for the time series plot.
- `xtickspacing`: Spacing between x-axis ticks for the time series plot.
- `colormap`: (Optional) Colormap to use for plotting (default: `:thermal`).

# Returns
- The modified axis object with the plotted lines and legend.

This function normalizes all time series to a common reference period (the first year), then plots
the SMB, FAC, discharge, observed, and modeled height anomalies as lines on the provided axis.
A legend is added for clarity.
"""
function plot_best_fit!(ax, best_smb, best_fac, best_discharge, best_fit, dv0; colormap=:thermal)

    ddate = dims(best_smb, :date)
    decyear = decimalyear.(collect(val(ddate)))
    clrs = Makie.resample_cmap(colormap, 6)

    best_smb = best_smb .- mean(best_smb[ddate[1]..(ddate[1]+Year(1))])
    best_fac = best_fac .- mean(best_fac[ddate[1]..(ddate[1]+Year(1))])
    best_discharge = best_discharge .- mean(best_discharge[ddate[1]..(ddate[1]+Year(1))])
    best_fit = best_fit .- mean(best_fit[ddate[1]..(ddate[1]+Year(1))])
    dv0 = dv0 .- mean(dv0[ddate[1]..(ddate[1]+Year(1))])

    ax.ylabel="height anomaly"
    ax.ytickformat=values -> ["$(round(Int,value))m" for value in values]
    #ax.xticks = xlims[1]:xtickspacing:xlims[2]
    #xlims!(ax, xlims)
    #title="best fit for $single_geotile_test [pscale = $pscale0, ΔT = $ΔT0, mad = $(round(cost_metric_minimum; digits=2))]"

    CairoMakie.lines!(ax, decyear, parent(dv0); label="observed", color=clrs[1])
    CairoMakie.lines!(ax, decyear, parent(best_fac); label="FAC", color=clrs[2])
    CairoMakie.lines!(ax, decyear, parent(best_discharge); label="discharge", color=clrs[3])
    CairoMakie.lines!(ax,decyear, parent(best_smb); label="SMB", color=clrs[4])
    CairoMakie.lines!(ax,decyear, parent(best_fit); label="modeled", color=clrs[5])

    axislegend(ax, position=:lb, patchsize=(20.0f0, 1.0f0), padding=(5.0f0, 5.0f0, 5.0f0, 5.0f0), labelsize=12, rowgap=1) # orientation=:horizontal, framevisible=false
    return ax
end

"""
    plot_cost_metric!(ax, cost_metric; colormap=:thermal)

Plot a contour map of the cost metric as a function of precipitation scale (pscale) and ΔT.

# Arguments
- `ax`: A CairoMakie Axis object to plot on.
- `cost_metric`: 2D array (or NamedDimsArray) of cost metric values, with dimensions :pscale and :ΔT.

# Keyword Arguments
- `colormap`: Colormap to use for the contour plot (default: `:thermal`).

# Returns
- `(ax, crange)`: The modified axis object and the color range tuple used for the plot.

This function creates a contour plot of the cost metric over the parameter space of precipitation scale and ΔT.
It highlights the minimum cost location with a marker and label, and adds a colorbar for reference.
"""
function plot_cost_metric!(ax, cost_metric; colormap=:thermal, step=1/35)

    dpscale = dims(cost_metric, :pscale)
    dΔT = dims(cost_metric, :ΔT)

    ax.xlabel = "precipitation scaling"
    ax.ylabel = "temperature offset [K]"

    # `pscale` is multiplicative, so it is drawn on the symmetric 1/n..n axis `scale_linear_ticks`
    # builds. `ΔT` is an additive offset in kelvin and is already linear: transforming it would put
    # 0 K at the axis value that means "no scaling" and compress the negative half.
    (x, xticks, xticklabels) = scale_linear_ticks(dpscale.val)
    y = collect(dΔT.val)

    # normalize cost metric to 0-1
    cost_metric .-= minimum(cost_metric[.!isinf.(cost_metric)])
    cost_metric ./= maximum(cost_metric[.!isinf.(cost_metric)])

    cost_metric[isinf.(cost_metric)] .= 9999
    crange = (0, 1)
   

    contour!(ax, x, y, cost_metric.data; colorrange=crange, levels=0:step:maximum(cost_metric), colormap)

    ax.xticks = (xticks, xticklabels)

    index_minimum = argmin(cost_metric)
    pscale_best = DimPoints(cost_metric)[index_minimum][1]
    ΔT_best = DimPoints(cost_metric)[index_minimum][2]

    x_best = scale2linear(pscale_best)
    y_best = ΔT_best
    scatter!(ax, x_best, y_best; color=:black, markersize=15, marker=:xcross)

    return ax, crange
end

"""
    plot_hist_gemb_altim_trend_amplitude(geotiles0)

Plot histograms comparing model and altimetry-derived trends and amplitudes.

# Arguments
- `geotiles0`: A DataFrame or table containing columns for model and altimetry trends and amplitudes.
    Expected columns: `"dv_trend"`, `"dv_altim_trend"`, `"dv_amplitude"`, `"dv_altim_amplitude"`.

# Returns
- `f`: A CairoMakie Figure object with two panels:
    - Left: Histogram of trends (model vs. altimetry)
    - Right: Histogram of amplitudes (model vs. altimetry)

# Description
Creates a two-panel figure. The first panel shows step histograms of the model and altimetry-derived trends.
The second panel shows step histograms of the model and altimetry-derived amplitudes. Legends are included for clarity.
"""
function plot_hist_gemb_altim_trend_amplitude(geotiles0)
    f = _publication_figure(columns=2, rows=1)
    var1 = "dh_trend"
    var2 = "dh_synthesis_trend"

    ax1 = CairoMakie.Axis(f[1, 1]; xlabel="trend", ylabel="count")
    ax1.xtickformat = values -> ["$(round(Int,value))m yr⁻¹" for value in values]
    CairoMakie.stephist!(geotiles0[:, var1], bins=-5:0.25:5; label="model")
    CairoMakie.stephist!(geotiles0[:, var2], bins=-5:0.25:5; label="observed")
    axislegend(ax1, framevisible=false)

    var1 = "dh_amplitude"
    var2 = "dh_synthesis_amplitude"
    ax2 = CairoMakie.Axis(f[1, 2]; xlabel="amplitude")
    ax2.xtickformat = values -> ["$(round(Int,value))m" for value in values]

    CairoMakie.stephist!(geotiles0[:, var1], bins=0:0.25:5; label="model")
    CairoMakie.stephist!(geotiles0[:, var2], bins=0:0.25:5; label="observed")
    axislegend(ax2, framevisible=false)

    return f
end



"""
    plot_multiregion_dvdm(regions; kwargs...) -> Figure, Vector, Vector

Create a multi-region plot of glacier mass change time series with stacked regions.

This function visualizes mass change time series for multiple glacier regions, stacking them
vertically with appropriate offsets to avoid overlap. It supports plotting multiple variables
with error bounds and customizable styling.

# Arguments
- `regions`: Dictionary containing mass change data for different regions and variables

# Keyword Arguments
- `variables=["dm", "dm_altim"]`: Variables to plot (last variable is plotted on top)
- `units="Gt"`: Units for mass change (e.g., "Gt", "m", "m w.e.")
- `rgi_regions=setdiff(collect(dims(runs_rgi["dm_altim"], :rgi)), [13, 14, 15, 99])`: RGI regions to include
- `showlines=false`: Whether to show grid lines
- `fontsize=15`: Font size for plot text
- `colormap=:Dark2_4`: Color map for regions
- `region_order=nothing`: Custom ordering of regions (default: sorted by total mass change)
- `ylims=nothing`: Custom y-axis limits
- `title=nothing`: Plot title
- `palette=nothing`: Custom color palette
- `delta_offset=nothing`: Vertical offset between regions (default depends on units)
- `all_error_bounds=false`: Whether to show error bounds for all variables
- `daterange=DateTime(2000,1,1):Month(1):DateTime(2025,1,1)`: Date range to plot
- `numbersinylabel=false`: Whether to include numeric values in y-axis labels

# Returns
- `f`: The Makie Figure object
- `region_order`: Vector of region IDs in the order they were plotted
- `ylims`: The y-axis limits used in the plot
"""
function plot_multiregion_dvdm(
    regions=regions;
    variables=["dm", "dm_altim"], # last variable is plotted last
    units="Gt",
    rgi_regions=setdiff(collect(dims(runs_rgi["dm_altim"], :rgi)), [13, 14, 15, 99]),
    showlines=false,
    fontsize=15,
    colormap=:Dark2_4,
    region_order=nothing,
    ylims=nothing,
    title=nothing,
    palette=nothing,
    delta_offset=nothing,
    all_error_bounds=false,
    daterange=DateTime(2000, 1, 1):Month(1):DateTime(2025, 1, 1),
    numbersinylabel=false,)

    drgi = Dim{:rgi}(rgi_regions)
    dvarname = Dim{:varname}(variables)

    if isnothing(delta_offset)
        if (units .== "m") || (units .== "m w.e.")
            delta_offset = -10
        else
            delta_offset = -100
        end
    end

    total_mass_change = zeros(drgi)

    #total mass change of feature variable (i.e. last variable in variables)
    first_valid = findfirst(.!isnan.(regions[dvarname[end]][At(drgi[1]), minimum(daterange)..maximum(daterange), At(false)]))
    last_valid = findlast(.!isnan.(regions[dvarname[end]][At(drgi[1]), minimum(daterange)..maximum(daterange), At(false)]))

    if isempty(first_valid)
        error("No valid data for $dvarname[end] in $daterange")
    end

    for rgi in drgi
        a = regions[dvarname[end]][At(rgi), minimum(daterange)..maximum(daterange), At(false)][first_valid]
        b = regions[dvarname[end]][At(rgi), minimum(daterange)..maximum(daterange), At(false)][last_valid]
        total_mass_change[At(rgi)] = b - a
    end

    if isnothing(region_order)
        region_order = drgi[sortperm(collect(total_mass_change[:]); rev=true)]
    end

    total_mass_change = total_mass_change[At(rgi_regions)]

    # make plots
    xtick_delta = 5
    years = unique(year.(daterange))

    yticks = zeros(length(rgi_regions))
    xticks = (floor(Int64, minimum(years) ./ xtick_delta)*xtick_delta):xtick_delta:(ceil(Int64, maximum(years) ./ xtick_delta)*xtick_delta)

    xlims = (xticks.start, xticks.stop)

    # backgroundcolor=RGBf(0.98, 0.98, 0.98)
    f = CairoMakie.Figure(; size=(500, 750), fontsize)

    # initialize
    n = length(rgi_regions)
    (n == 1) ? (n = 2) : (n = n)

    if isnothing(palette)
        palette = (; color=Makie.resample_cmap(colormap, n))
    end

    y2ticks = zeros(n)
    y2ticklabels = fill("", n)
    y2ticklabelcolor = palette.color

    yticklabels = rginum2label.(region_order)

    if numbersinylabel
        yticklabels = ["$(rginum2label(id)) $(rginum2enclosed_alphanumerics[id])" for id in region_order]
    end

    # this is a hack... axes need to be defined early or things break
    ax1 = CairoMakie.Axis(f[1, 1];
        #title,
        palette,
        xticks=WilkinsonTicks(5),
        yticklabelrotation=0 / (360 / (2 * pi)),
        ygridcolor=:black,
        ygridwidth=0.5,
        yticks=((1:n) / n, string.((1:n) / n)),
        yticklabelcolor=y2ticklabelcolor,
    )

    if .!isnothing(title)
        ax1.title = title
        ax1.titlecolor = RGBf(0.85, 0.85, 0.85) / 2
    end

    ax2 = CairoMakie.Axis(f[1, 1];
        yaxisposition=:right,
        yticks=((1:n) / n, string.((1:n) / n)),
        yticklabelcolor=y2ticklabelcolor,
    )

    CairoMakie.xlims!(ax1, xlims...)

    ymin = Inf
    ymax = -Inf

    region_offsets = zeros(region_order)

    # calculate offsets 
    varname = dvarname[end]
    valid_index = .!isnan.(regions[varname][At(region_order[1]), minimum(daterange)..maximum(daterange), At(false)])
    last = zeros(sum(valid_index))

    for rgi in region_order
        mid0 = regions[varname][At(rgi), minimum(daterange)..maximum(daterange), At(false)][valid_index]
        region_offsets[At(rgi)] = floor(min(minimum((last .- (mid0 .- mid0[1]))), 0) / delta_offset) * delta_offset .+ delta_offset .- mid0[1]

        last = mid0 .+ region_offsets[At(rgi)]
    end

    # plot error bounds for
    endi = length(dvarname)
    if all_error_bounds
        starti = 1
    else
        starti = length(dvarname)
    end

    for varname in dvarname[starti:endi]
        valid_index = .!isnan.(regions[varname][At(region_order[1]), minimum(daterange)..maximum(daterange), At(false)])
        last = zeros(sum(valid_index))

        for rgi in region_order

            mid0 = regions[varname][At(rgi), minimum(daterange)..maximum(daterange), At(false)][valid_index]

            if all(isnan.(mid0))
                continue
            end

            mid = mid0 .+ region_offsets[At(rgi)]
            low = mid .- regions[varname][At(rgi), minimum(daterange)..maximum(daterange), At(true)][valid_index]
            high = mid .+ regions[varname][At(rgi), minimum(daterange)..maximum(daterange), At(true)][valid_index]

            CairoMakie.band!(ax1, collect(decimalyear.(dims(mid0, :date))), collect(low), collect(high); color=(:grey, 0.4))

            ymin = min(ymin, minimum(low))
            ymax = max(ymax, maximum(high))

        end
    end

    # plot lines over error bounds
    for varname in dvarname
        #varname = "dm_altim"
        valid_index = .!isnan.(regions[varname][At(region_order[1]), minimum(daterange)..maximum(daterange), At(false)])
        last = zeros(sum(valid_index))

        c = 0
        for (i, rgi) in enumerate(region_order)

            #rgi = region_order[1]
            mid0 = regions[varname][At(rgi), minimum(daterange)..maximum(daterange), At(false)][valid_index] .+ region_offsets[At(rgi)]

            if all(isnan.(mid0))
                continue
            end

            c += 1
            if varname == dvarname[end]
                ln = CairoMakie.lines!(ax1, collect(decimalyear.(dims(mid0, :date))), collect(mid0))

                if !occursin("yr⁻¹", units)
                    y2ticklabels[i] = string(round(Int64, (mid0[end] - mid0[1]))) * " $(units)"
                else
                    fit = ts_seasonal_model(mid0)
                    y2ticklabels[i] = string(round(Int64, fit.trend)) * " $(units)"
                end

                yticks[i] = mid0[1]
                y2ticks[i] = mid0[end]

                y2ticklabelcolor[i] = palette.color[c]
            else
                ln = CairoMakie.lines!(ax1, collect(decimalyear.(dims(mid0, :date))), collect(mid0); color=(:black, 0.3))
            end

        end
    end

    if isnothing(ylims)
        ylims = [ymin, ymax]
    end

    CairoMakie.ylims!(ax1, minimum(ylims), maximum(ylims))
    CairoMakie.reset_limits!(ax1) # this is needed to be able to retrive limits

    minorgridspacing = -delta_offset / 5
    ylim = round.(Int, ax1.yaxis.attributes.limits.val ./ minorgridspacing) .* minorgridspacing
    #yticks = collect(region_offsets[At(collect(region_order))])

    ax1.yticks = (yticks, yticklabels)
    ax1.yticklabelcolor = y2ticklabelcolor
    ax1.ygridvisible = showlines

    #!showlines && CairoMakie.hidexdecorations!(ax1,)
    try
        ax2.yticks = (y2ticks, y2ticklabels)
    catch e
        printstyled("NOTE: if you end up here there is a bug with CairoMakie in which a cryptic error is thrown when the ytick positions exceed the yaxis limits\n"; color=:red)
        rethrow(e)
    end

    ax2.yticklabelcolor = y2ticklabelcolor

    foo = ax1.yaxis.attributes.limits.val
    CairoMakie.ylims!(ax2, ylims)
    CairoMakie.ylims!(ax2, ylims)

    minorgridspacing = -(delta_offset / 5)

    ax2.yminorticks = ylim[1]:minorgridspacing:ylim[2]

    !showlines && CairoMakie.hidespines!(ax1, :t, :r, :l)
    CairoMakie.hidespines!(ax2)
    ax2.yminorticksvisible = false
    ax2.yminorgridvisible = showlines
    ax2.yticksvisible = false
    ax2.ygridvisible = false

    CairoMakie.hidexdecorations!(ax2)

    xtickcolor = RGBf(0.85, 0.85, 0.85)
    ax1.xaxisposition = :bottom
    ax1.xgridvisible = true
    ax1.xgridcolor = xtickcolor
    ax1.xticklabelcolor = xtickcolor ./ 2
    ax1.xticks = (Float64.(xticks), string.(xticks))
    ax1.xtickcolor = xtickcolor
    ax1.xlabelvisible = true
    ax1.yticksvisible = false
    ax1.xticksvisible = true
    ax1.bottomspinecolor = xtickcolor

    CairoMakie.ylims!(f.content[1], ylims)

    return f, region_order, ylims
end


"""
    show_error_bar_table(regional_results; cols2display = 3)

Display a formatted table of regional results in the console, organized in sections.

Prints the regional results DataFrame in multiple sections, with each section showing a subset
of columns. The first two columns (typically region identifiers) are always displayed, followed
by `cols2display` data columns in each section. Section headers are printed in yellow and
contain the name of the first column in that section.

# Arguments
- `regional_results`: DataFrame containing regional data to be displayed
- `cols2display`: Number of data columns to show in each section (default: 3)
"""
function show_error_bar_table(regional_results; cols2display=3)
    # dump table to console
    colnames = names(regional_results)
    last_set = ceil(Int, (length(colnames) - 2) / (cols2display + 1) - 1)
    for i = 0:last_set
        printstyled("\n------------------------------------------------------- $(colnames[3 + (4 * i)])-------------------------------------------------------\n", color=:yellow)
        if (i == last_set) && (length(colnames) > (last_set + 4))
            show(regional_results[:, [1:2; (3 .+ (4*i)):end]], allrows=true, allcols=true)
        else
            show(regional_results[:, [1:2; (3:(2+cols2display)) .+ (4 * i)]], allrows=true, allcols=true)
        end
    end
end




"""
    plot_point_location_river_flux(land_flux, glacier_flux; dates4plot, COMID, name)

Plot river flux components (land and glacier) for a specific location, without snow flux.

This is a convenience wrapper for `plot_point_location_river_flux` that omits the snow flux argument.

# Arguments
- `land_flux`: DimArray of land river flux (time × COMID).
- `glacier_flux`: DimArray of glacier river flux (time × COMID).
- `dates4plot`: Tuple of (start_date, end_date) for plotting.
- `COMID`: Integer river reach identifier.
- `name`: String for plot title.

# Returns
- `f`: Makie Figure object.
"""
function plot_point_location_river_flux(land_flux, glacier_flux; date_range, COMID, name, show_title=true)
    f = plot_point_location_river_flux(land_flux, glacier_flux, nothing; date_range, COMID, name, show_title)
    return f
end

"""
    plot_point_location_river_flux(land_flux, glacier_flux, snow_flux; dates4plot, COMID, name)

Plot river flux components and glacier contribution for a specific location.

Creates a figure with:
- Stacked area plots of land, glacier, and (optionally) snow river fluxes over time for a given COMID.
- Glacier fraction of total river flux, with seasonal maximum (gmax) highlighted.

# Arguments
- `land_flux`: DimArray of land river flux (time × COMID).
- `glacier_flux`: DimArray of glacier river flux (time × COMID).
- `snow_flux`: DimArray of snow river flux (time × COMID), or `nothing` if not available.
- `dates4plot`: Tuple of (start_date, end_date) for plotting.
- `COMID`: Integer river reach identifier.
- `name`: String for plot title.

# Returns
- `f`: Makie Figure object.
"""
function plot_point_location_river_flux(land_flux, glacier_flux, snow_flux; date_range, COMID, name, show_title=true)
   

    f = _publication_figure(columns=1, rows=1)

    date_interval = date_range[1]..date_range[2]

    land_flux0 = land_flux[Ti=date_interval, COMID=At(COMID)]
    glacier_flux0 = glacier_flux[Ti=date_interval, COMID=At(COMID)]

    if !isnothing(snow_flux)
         separate_out_snow = true
         snow_flux0 = snow_flux[Ti = date_interval, COMID=At(COMID)]
    else
        separate_out_snow = false
    end

    x = decimalyear.(val(dims(land_flux0, :Ti)))
    xticks = collect(floor(minimum(x)/5)*5:5:ceil(maximum(x)/5)*5)

    separate_out_snow && (y_snow = parent(snow_flux0) ./ 1000)

    ax1 = CairoMakie.Axis(f[1:2, 1]; xticklabelsvisible=false, xticksvisible=false, xticks=xticks)
    if show_title
        ax1.title = name
    end

    # plot land flux
    y = parent(land_flux0) ./ 1000
    separate_out_snow && (y .-= y_snow)

    x0 = vcat(x[1], x, x[end])
    yunits = Unitful.unit(land_flux[1])
    y = vcat(0 * yunits, y, 0 * yunits)
    max_y = maximum(y)
    poly!(ax1, GI.Point.(ustrip(x0), ustrip(y)); color=(:peru, 1), label="land")

    # plot snow flux
    if separate_out_snow
        y = y_snow
        yunits = Unitful.unit(y[1])
        y = vcat(0.0 * yunits, y, 0.0 * yunits)
        max_y = max(maximum(y), max_y)
        poly!(ax1, GI.Point.(x0, ustrip(y)); color=(:gray, 0.5), label="snow")
    end

    y = parent(glacier_flux0) ./ 1000
    yunits = Unitful.unit(y[1])
    y = vcat(0 * yunits, y, 0 * yunits)
    poly!(ax1, GI.Point.(ustrip(x0), ustrip(y)); color=(:blue, 0.5), label="glacier")

    axislegend(ax1; merge=true, backgroundcolor=(:white, 0), framevisible=false, position=:lt)
    ax1.ylabel = "river flux (m³/s) x 10⁻³"
    xlims!(ax1, minimum(xticks), maximum(xticks))
    ylims!(ax1, 0, ustrip(max_y))

    ax2 = CairoMakie.Axis(f[3, 1], xticks=xticks)
    glacier_frac = glacier_flux0 ./ (land_flux0 .+ glacier_flux0) * 100

    gmax = groupby(glacier_frac, Ti => Bins(month, 1:12))
    gmax = mean.(gmax; dims=:Ti)
    gmax = cat(gmax...; dims=dims(gmax, :Ti))
    dTi = dims(gmax, :Ti)
    gmax_month = dTi[argmax(gmax)]
    gmax = round(Int8, maximum(gmax))

    foo = glacier_frac[month.(dims(glacier_frac, :Ti)).==gmax_month]
    gmax_min = round(Int8, minimum(foo))
    gmax_max = round(Int8, maximum(foo))

    lines!(ax2, x, parent(glacier_frac); color=(:blue, 0.5))
    lines!(ax2, [minimum(x), maximum(x)], [gmax_min, gmax_min]; color=(:gray, 0.5), linestyle=:dash)
    lines!(ax2, [minimum(x), maximum(x)], [gmax_max, gmax_max]; color=(:gray, 0.5), linestyle=:dash)
    lines!(ax2, [minimum(x), maximum(x)], [gmax, gmax]; color=(:black, 1))
    ax2.ylabel = "glacier fraction [%]"

    xlims!(ax2, minimum(xticks), maximum(xticks))
    text!(maximum(x), mean((gmax_max, gmax)); text="gmax = $gmax [$(Dates.monthname(gmax_month))] ", align=(:right, :center), color=(:black, 0.8), overdraw=true)

    return f
end


function plot_ref_pscale_ΔT_summary(path2runs_synthesized, binned_synthesized_dv_file_ref; rgi2plot=[99], seasonality_weight=seasonality_weight, distance_from_origin_penalty=distance_from_origin_penalty, ΔT_to_pscale_weight=ΔT_to_pscale_weight, show_title=true, bin_width=0.25)

    # The reference run's calibrated synthesis is in its physics member's folder; the altimetry
    # member it was synthesized from is not.
    dv_ref, physics_member = gemb_physics_source(binned_synthesized_dv_file_ref)
    ensemble_reference_file = replace(dv_ref, "_gembfit_dv.jld2" => ".jld2")

    region_fits = ensemble_summary(path2runs_synthesized, ensemble_reference_file; physics_member)
    gemb_fit = GeoDataFrames.read(gembfit_path(ensemble_reference_file, physics_member))

    title1 = "W_s: $(seasonality_weight), W_d: $(distance_from_origin_penalty), W_m2p: $(ΔT_to_pscale_weight)"

    # `pscale` is multiplicative and is binned on the symmetric 1/n..n axis. `ΔT` is an additive offset
    # in kelvin, so it is binned in kelvin directly -- transforming it would place 0 K where the scale
    # axis means "no scaling" and compress the negative half against the positive one.
    (pscale_linear, pscale_ticks, pscale_ticklabels) = scale_linear_ticks(gemb_fit[:, :pscale])

    e0 = ceil(Int, maximum(abs.(pscale_linear)))
    pscale_linear_bins = (-e0-bin_width):bin_width:(e0+bin_width)
    pscale_linear_bin_centers = (pscale_linear_bins[1:end-1] .+ pscale_linear_bins[2:end]) ./ 2

    e0 = ceil(Int, maximum(abs.(gemb_fit[:, :ΔT])))
    ΔT_bins = (-e0-bin_width):bin_width:(e0+bin_width)
    ΔT_bin_centers = (ΔT_bins[1:end-1] .+ ΔT_bins[2:end]) ./ 2


    f = Vector{Any}(undef, length(rgi2plot))
    for i in eachindex(rgi2plot)
        rgi = rgi2plot[i]
        runoff = round.(Int, region_fits[varname=At("runoff"), parameter=At("trend"), rgi=At(rgi), error=At(false)])
        acc = round.(Int, region_fits[varname=At("acc"), parameter=At("trend"), rgi=At(rgi), error=At(false)])
        fac = round.(Int, region_fits[varname=At("fac"), parameter=At("trend"), rgi=At(rgi), error=At(false)])

        title2 = "runoff: $(runoff), acc: $(acc), fac:$(fac)"

        title = "$(rginum2label(rgi)) [$(title2)] Gt/yr: $(title1)"

        if rgi == 99
            index = gemb_fit.rgi .> 0
        elseif rgi == 98
            index = (gemb_fit.rgi .> 12) .& (gemb_fit.rgi .< 16)
        else
            index = gemb_fit.rgi .== rgi
        end
        
        h = StatsBase.fit(Histogram, scale2linear.(gemb_fit[index, :pscale]), pscale_linear_bins)
        f[i] = _publication_figure(; columns=2, rows=1)
        if show_title
            f[i][0, 1:2] = Label(f[i], title)
        end

        ax = Makie.Axis(f[i][1, 1]; xlabel="precipitation scaling", ylabel="count", xticks=(pscale_ticks, pscale_ticklabels))
        barplot!(ax, pscale_linear_bin_centers, h.weights;)

        h = StatsBase.fit(Histogram, gemb_fit[index, :ΔT], ΔT_bins)
        ax = Makie.Axis(f[i][1, 2]; xlabel="temperature offset [K]", ylabel="count")
        barplot!(ax, ΔT_bin_centers, h.weights;)
    end

    return f
end


function plot_ensemble_rgi_timeseries(vars2plot, path2runs_synthesized; center_on_dates=DateTime(2002, 5, 1):Month(1):DateTime(2005, 12, 15), rgi2plot=nothing, physics_member=gemb_reference_physics_member)

    # load results for all runs for each RGI [18s]
    runs_rgi = runs2rgi(path2runs_synthesized; physics_member)

    dvarname = Dim{:varname}(vars2plot)

    # center data
    runs_rgi = runs_center!(runs_rgi, center_on_dates)
    check_for_all_nans(runs_rgi)

    clrs = DimArray(Makie.resample_cmap(:thermal, length(dvarname)), dvarname; name="colors")

    if isnothing(rgi2plot)
        rgi2plot = dims(runs_rgi[vars2plot[1]], :rgi)
    end

    f = Vector{Any}(undef, length(rgi2plot))
    for i in eachindex(rgi2plot)
        rgi = rgi2plot[i]
        f[i] = _publication_figure(; columns=1, rows=1)
        ax = Makie.Axis(f[i][1, 1]; title="$(rginum2label(rgi))", ylabel="mass/volume anomaly [Gt/km³]")
        for i = eachindex(path2runs_synthesized)
            for var2plot in vars2plot
                decyear = decimalyear.(collect(dims(runs_rgi[var2plot][run=At(path2runs_synthesized[i]), rgi=At(rgi)], :date).val))
                lines!(ax, decyear, runs_rgi[var2plot][run=At(path2runs_synthesized[i]), rgi=At(rgi)].data; label=var2plot, color=clrs[varname=At(var2plot)])
            end
        end
        axislegend(ax, unique=true)
    end

    return f
end


function scale_linear_ticks(scale0)
   
    xlinear = scale2linear.(scale0)

    xticks = round.([-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6], digits=4)
    xticklabels = ["1/7", "1/6", "1/5", "1/4", "1/3", "1/2", "1", "2", "3", "4", "5", "6", "7"]

    xlims = extrema(xlinear)
    index = (xticks .<= ceil(Int, maximum(xlims))) .& (xticks .>= floor(Int,    minimum(xlims)))
    xticks = xticks[index]
    xticklabels = xticklabels[index]

    return (xlinear, xticks, xticklabels)
end


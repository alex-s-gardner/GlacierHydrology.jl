"""
Plot calibrated mass balance components for a given geotile

This script loads GEMB data, discharge data, adjusts GEMB for discharge,
and plots the mass balance components for a specified geotile.
"""

import GlobalGlacierAnalysis as GGA
import GlacierVolumeChange as GVC
using CairoMakie
using FileIO
using DataFrames
using Dates
using Unitful
Unitful.register(GGA.MyUnits)
using Statistics

# Configuration
geotile2plot = "lat[-44-42]lon[+170+172]"
plots_show = true
plots_save = true
gemb_run_id = 9

# Load GEMB ensemble data
println("Loading GEMB ensemble data...")
surface_mask = GGA.reference_surface_mask()
gemb_physics_member = GGA.gemb_reference_physics_member
gemb = GGA.gemb_ensemble_dv(; gemb_run_id, surface_mask, physics_member=gemb_physics_member)

# Load discharge data
println("Loading discharge data...")
discharge = GGA.global_discharge_filled(;
    surface_mask,
    discharge_global_fn=GGA.discharge_global_file(surface_mask),
    gemb,
    discharge2smb_max_latitude=-60,
    discharge2smb_equilibrium_period=(Date(1979), Date(2000)),
    pscale=1,
    ΔT=0,
    geotile_width=2,
    force_remake_before=GGA.gemb_ensemble_built(GGA.gemb_info(; gemb_run_id, physics_member=gemb_physics_member), [surface_mask]),
    force_remake_before_hypsometry=nothing
)

# Adjust GEMB for discharge
println("Adjusting GEMB for discharge...")
gemb = GGA.dv_adjust4discharge!(gemb, discharge)

# Get reference ensemble file and construct path to synthesized data
reference_ensemble_file = GVC.reference_ensemble_file
binned_synthesized_file = replace(reference_ensemble_file, "aligned.jld2" => "synthesized.jld2")
binned_synthesized_dv_file = GGA.gemb_physics_path(replace(binned_synthesized_file, ".jld2" => "_gembfit_dv.jld2"),
                                                   gemb_physics_member)

println("Loading mass balance components from: $binned_synthesized_dv_file")

# Load the DataFrame
geotiles_df = load(binned_synthesized_dv_file, "geotiles")

# Check if the geotile exists in the data
geotile_ids = geotiles_df[!, :id]
if !(geotile2plot in geotile_ids)
    println("Geotile $geotile2plot not found in data.")
    println("Total available geotiles: $(length(geotile_ids))")
    println("First 20 available geotiles:")
    for (i, gt) in enumerate(geotile_ids)
        if i > 20
            break
        end
        println("  $gt")
    end
    error("Geotile not found")
end

# Extract row for the specified geotile
geotile_row = geotiles_df[geotiles_df[!, :id] .== geotile2plot, :][1, :]

# Define mass balance components to plot
component_names = ["smb", "acc", "melt", "runoff", "refreeze", "rain", "fac", "ec"]

# Get dates from column metadata for the first variable
dates = colmetadata(geotiles_df, :smb, "date")
time_years = year.(dates) .+ (month.(dates) .- 1) ./ 12

# Define colors for different components
component_colors = Dict(
    "smb" => :blue,
    "acc" => :cyan,
    "melt" => :red,
    "runoff" => :orange,
    "refreeze" => :purple,
    "rain" => :green,
    "fac" => :brown,
    "ec" => :magenta,
    "dv" => :black,
    "dm" => :gray
)

component_labels = Dict(
    "smb" => "Surface Mass Balance",
    "acc" => "Accumulation",
    "melt" => "Melt",
    "runoff" => "Runoff",
    "refreeze" => "Refreeze",
    "rain" => "Rain",
    "fac" => "Firn Air Content",
    "ec" => "Evaporation/Condensation",
    "dv" => "Volume Change",
    "dm" => "Mass Change"
)

# Get conversion factor from km³ ice equivalent to m w.e.
# Values are stored in km³ i.e. per geotile
# To convert to m w.e. rate:
# 1. diff gives km³/month
# 2. Divide by area (km²) to get km/month = 1000 m/month
# 3. Convert ice to water equivalent: multiply by ρ_ice/ρ_water = δice/1000
mie2cubickm = geotile_row[:mie2cubickm]  # Conversion factor: area_km2 / 1000
area_km2 = mie2cubickm * 1000  # Geotile glacier area in km²

println("Geotile area: $(area_km2) km²")
println("Conversion factor (mie2cubickm): $(mie2cubickm)")

# Collect valid components and their data
plot_data = []
for varname in component_names
    # Check if variable exists in DataFrame
    if !(Symbol(varname) in propertynames(geotiles_df))
        println("Skipping $varname - not found in data")
        continue
    end

    # Extract values for this geotile
    values = geotile_row[Symbol(varname)]

    # Skip if missing
    if ismissing(values) || isnothing(values)
        println("Skipping $varname - no data")
        continue
    end

    # Calculate rate:
    # values are in km³ i.e. (cumulative)
    # diff gives km³ i.e./month
    # Divide by mie2cubickm to get m i.e./month
    # Multiply by δice/1000 to convert ice to water equivalent
    rates_km3_per_month = diff(values)  # km³ i.e./month
    rates = rates_km3_per_month ./ mie2cubickm .* (GGA.δice / 1000)  # m w.e./month

    # Use midpoint times for the rates
    time_years_mid = (time_years[1:end-1] .+ time_years[2:end]) ./ 2

    # Get color and label
    color = get(component_colors, varname, :black)
    label = get(component_labels, varname, varname)

    push!(plot_data, (varname=varname, label=label, color=color, times=time_years_mid, rates=rates))
end

# Create figure with subplots
n_plots = length(plot_data)
n_cols = 2
n_rows = ceil(Int, n_plots / n_cols)

fig = Figure(size=(1400, 300 * n_rows))

# Add overall title
Label(fig[0, :], "Mass Balance Component Rates for Geotile: $geotile2plot",
      fontsize=20, font=:bold)

# Create subplots
for (idx, data) in enumerate(plot_data)
    row = div(idx - 1, n_cols) + 1
    col = mod(idx - 1, n_cols) + 1

    ax = Axis(fig[row, col],
        xlabel="Year",
        ylabel="Rate (m/month w.e.)",
        title=data.label
    )

    # Plot line
    lines!(ax, data.times, data.rates,
        color=data.color,
        linewidth=1.5)

    # Add horizontal line at zero
    hlines!(ax, [0], color=:gray, linestyle=:dash, linewidth=1)
end

# Display and/or save
if plots_show
    display(fig)
end

if plots_save
    output_dir = GGA.pathlocal[:figures]
    mkpath(output_dir)
    # Replace special characters in filename
    filename_safe = replace(geotile2plot, "[" => "_", "]" => "_", "+" => "p", "-" => "m")
    output_file = joinpath(output_dir, "geotile_mass_balance_components_$(filename_safe).png")
    save(output_file, fig, px_per_unit=2)
    println("Figure saved to: $output_file")
end

# Create second figure: dv comparison
println("\nCreating dv comparison plot...")

# Check if all required variables exist
required_vars = [:dv, :acc, :runoff, :ec, :fac, :dv_altim]
if all(v -> Symbol(v) in propertynames(geotiles_df), required_vars)

    # Extract values
    dv_values = geotile_row[:dv]
    dv_altim_values = geotile_row[:dv_altim]
    acc_values = geotile_row[:acc]
    runoff_values = geotile_row[:runoff]
    ec_values = geotile_row[:ec]
    fac_values = geotile_row[:fac]

    # Get dates for dv_altim (may be different from GEMB dates)
    dates_altim = colmetadata(geotiles_df, :dv_altim, "date")
    time_years_altim = year.(dates_altim) .+ (month.(dates_altim) .- 1) ./ 12

    # Calculate rates for dv
    dv_rates_km3 = diff(dv_values)
    dv_rates = dv_rates_km3 ./ mie2cubickm .* (GGA.δice / 1000)

    # Calculate rates for dv_altim
    dv_altim_rates_km3 = diff(dv_altim_values)
    dv_altim_rates = dv_altim_rates_km3 ./ mie2cubickm .* (GGA.δice / 1000)
    time_years_altim_mid = (time_years_altim[1:end-1] .+ time_years_altim[2:end]) ./ 2

    # Calculate rates for components
    acc_rates = diff(acc_values) ./ mie2cubickm .* (GGA.δice / 1000)
    runoff_rates = diff(runoff_values) ./ mie2cubickm .* (GGA.δice / 1000)
    ec_rates = diff(ec_values) ./ mie2cubickm .* (GGA.δice / 1000)
    fac_rates = diff(fac_values) ./ mie2cubickm .* (GGA.δice / 1000)

    # Calculate the sum: acc - runoff - ec + fac
    calculated_dv = acc_rates .- runoff_rates .- ec_rates .+ fac_rates

    # Time axis
    time_years_mid = (time_years[1:end-1] .+ time_years[2:end]) ./ 2

    # Create figure
    fig2 = Figure(size=(1200, 600))

    Label(fig2[0, :], "Volume Change Comparison for Geotile: $geotile2plot",
          fontsize=20, font=:bold)

    # Top panel: time series
    ax1 = Axis(fig2[1, 1],
        xlabel="Year",
        ylabel="Rate (m/month w.e.)",
        title="dv vs Calculated (acc - runoff - ec + fac)"
    )

    lines!(ax1, time_years_mid, dv_rates,
        label="dv (GEMB model)",
        color=:blue,
        linewidth=2)

    lines!(ax1, time_years_mid, calculated_dv,
        label="acc - runoff - ec + fac",
        color=:red,
        linestyle=:dash,
        linewidth=2)

    lines!(ax1, time_years_altim_mid, dv_altim_rates,
        label="dv_altim (synthesis)",
        color=:green,
        linewidth=2,
        alpha=0.7)

    hlines!(ax1, [0], color=:gray, linestyle=:dash, linewidth=1)
    axislegend(ax1, position=:lt)

    # Bottom panel: residual
    ax2 = Axis(fig2[2, 1],
        xlabel="Year",
        ylabel="Residual (m/month w.e.)",
        title="Residual (dv - calculated)"
    )

    residual = dv_rates .- calculated_dv
    lines!(ax2, time_years_mid, residual,
        color=:black,
        linewidth=1.5)

    hlines!(ax2, [0], color=:gray, linestyle=:dash, linewidth=1)

    # Calculate statistics
    rmse = sqrt(mean(residual.^2))
    mean_residual = mean(residual)
    println("Residual statistics:")
    println("  RMSE: $(round(rmse, digits=6)) m/month w.e.")
    println("  Mean: $(round(mean_residual, digits=6)) m/month w.e.")

    # Add text with statistics
    text!(ax2, 0.02, 0.95,
        text="RMSE: $(round(rmse, digits=6)) m/month\nMean: $(round(mean_residual, digits=6)) m/month",
        align=(:left, :top),
        space=:relative,
        fontsize=12)

    # Display and/or save
    if plots_show
        display(fig2)
    end

    if plots_save
        output_file2 = joinpath(output_dir, "geotile_dv_comparison_$(filename_safe).png")
        save(output_file2, fig2, px_per_unit=2)
        println("DV comparison figure saved to: $output_file2")
    end
else
    println("Skipping dv comparison - missing required variables")
end

# Create third figure: acc, rain, and acc - rain
println("\nCreating acc/rain comparison plot...")

if all(v -> Symbol(v) in propertynames(geotiles_df), [:acc, :rain])

    acc_values = geotile_row[:acc]
    rain_values = geotile_row[:rain]

    # Calculate rates in m w.e./month
    acc_rates = diff(acc_values) ./ mie2cubickm .* (GGA.δice / 1000)
    rain_rates = diff(rain_values) ./ mie2cubickm .* (GGA.δice / 1000)
    acc_minus_rain_rates = acc_rates .- rain_rates

    time_years_mid = (time_years[1:end-1] .+ time_years[2:end]) ./ 2

    fig3 = Figure(size=(1200, 900))

    Label(fig3[0, :], "Accumulation and Rain for Geotile: $geotile2plot",
          fontsize=20, font=:bold)

    # Subplot 1: acc
    ax1 = Axis(fig3[1, 1],
        ylabel="Rate (m w.e./month)",
        title="Accumulation (acc)"
    )
    lines!(ax1, time_years_mid, acc_rates, color=:cyan, linewidth=1.5)
    hlines!(ax1, [0], color=:gray, linestyle=:dash, linewidth=1)

    # Subplot 2: rain
    ax2 = Axis(fig3[2, 1],
        ylabel="Rate (m w.e./month)",
        title="Rain"
    )
    lines!(ax2, time_years_mid, rain_rates, color=:green, linewidth=1.5)
    hlines!(ax2, [0], color=:gray, linestyle=:dash, linewidth=1)

    # Subplot 3: acc - rain
    ax3 = Axis(fig3[3, 1],
        xlabel="Year",
        ylabel="Rate (m w.e./month)",
        title="Accumulation - Rain (snow only)"
    )
    lines!(ax3, time_years_mid, acc_minus_rain_rates, color=:blue, linewidth=1.5)
    hlines!(ax3, [0], color=:gray, linestyle=:dash, linewidth=1)

    # Link x-axes
    linkxaxes!(ax1, ax2, ax3)
    hidexdecorations!(ax1, grid=false)
    hidexdecorations!(ax2, grid=false)

    if plots_show
        display(fig3)
    end

    if plots_save
        output_file3 = joinpath(output_dir, "geotile_acc_rain_$(filename_safe).png")
        save(output_file3, fig3, px_per_unit=2)
        println("Acc/rain figure saved to: $output_file3")
    end
else
    println("Skipping acc/rain plot - missing required variables")
end

println("\nDone!")

# Reproduce three GEMB-pipeline audit findings from the production data, read-only.
#
#   1. The geotile `discharge` column was divided by δice/1000 a second time, after
#      `discharge2geotile(...; mass2volume=true)` had already converted Gt to km³ ice equivalent.
#   2. Antarctic discharge from `discharge2smb` (SMB trend × area) was km³ ice equivalent per year
#      stored as `discharge_gtyr`, so the Antarctic total read ~9 % high.
#   3. The GEMB ensemble was area-weighted by RGI6 hypsometry and calibrated against altimetry
#      weighted by RGI7 hypsometry for `glacier_rgi7` members.
#
# Run from the repository root (needs the JPL data directories in `src/local_paths.jl`):
#
#     ~/.juliaup/bin/julia +release --project=. scripts/audit/discharge_units_and_surface_masks.jl
#
# Nothing is written. Output with the discharge file of 2026-10-05 and the hypsometry of 2025-12:
#
#     Antarctic discharge (south of 60°S, 2182 glaciers)
#       stored as Gt/yr:              110.6
#       in Gt/yr (× δice/1000):       101.4
#       80–115 Gt/yr gate passes:    stored true, corrected true
#     Geotile discharge column, sum over geotiles, per year (km³ i.e./yr)
#       before the fixes:              199.5
#       double conversion only:        188.6   (× 1.0905 = 1000/δice)
#       corrected:                     172.9
#     RGI7 / RGI6 glacier area per geotile (877 geotiles with ice in both)
#       per-geotile quantiles    5 %: 0.67   50 %: 1.00   95 %: 1.22
#       area-weighted quantiles  5 %: 0.93   50 %: 1.00   95 %: 1.08
#       ice area in geotiles where the ratio differs from 1 by more than 10 %: 7.3 %
#       geotiles with ice in only one inventory: RGI6 27, RGI7 0
#
# The discharge file is whatever `global_discharge_filled` last wrote to
# `pathlocal[:discharge_global]`, so the values move when it is regenerated.

import GlobalGlacierAnalysis as GGA
using DataFrames
using DimensionalData
using FileIO
using Printf
using Statistics

discharge_file = GGA.pathlocal[:discharge_global]
δice = GGA.δice
volume2mass = δice / 1000

# --- 1 and 2: discharge units -------------------------------------------------------------------
discharge = FileIO.load(discharge_file, "discharge")
antarctic = discharge.latitude .< -60

# South of 60°S every row comes from `discharge2smb`, whose stored value was km³ i.e./yr.
antarctic_stored = sum(discharge.discharge_gtyr[antarctic])
antarctic_gt = antarctic_stored * volume2mass

println("discharge file: $(discharge_file)")
println("Antarctic discharge (south of 60°S, $(count(antarctic)) glaciers)")
@printf("  stored as Gt/yr:             %6.1f\n", antarctic_stored)
@printf("  in Gt/yr (× δice/1000):      %6.1f\n", antarctic_gt)
@printf("  80–115 Gt/yr gate passes:    stored %s, corrected %s\n",
        80 <= antarctic_stored <= 115, 80 <= antarctic_gt <= 115)

# Corrected per-glacier table: Antarctic rows converted to Gt.
discharge_fixed = copy(discharge)
discharge_fixed[antarctic, :discharge_gtyr] .*= volume2mass

dgeotile = dims(GGA._geotile_area_km2(; surface_mask="glacier", geotile_width=2), :geotile)
rate_old = GGA.discharge2geotile(discharge, dgeotile; mass2volume=true)[:discharge]
rate_new = GGA.discharge2geotile(discharge_fixed, dgeotile; mass2volume=true)[:discharge]

# The column per year of elapsed time. Before the fixes the km³ rate was divided by volume2mass a
# second time, on top of the Antarctic rows already being km³ stored as Gt.
column_both = sum(rate_old) / volume2mass
column_double_only = sum(rate_new) / volume2mass
column_new = sum(rate_new)
println("Geotile discharge column, sum over geotiles, per year (km³ i.e./yr)")
@printf("  before the fixes:            %7.1f\n", column_both)
@printf("  double conversion only:      %7.1f   (× %.4f = 1000/δice)\n", column_double_only, column_double_only / column_new)
@printf("  corrected:                   %7.1f\n", column_new)

# --- 3: RGI6 vs RGI7 area per geotile -----------------------------------------------------------
area6 = dropdims(sum(GGA._geotile_area_km2(; surface_mask="glacier", geotile_width=2); dims=:height); dims=:height)
area7 = dropdims(sum(GGA._geotile_area_km2(; surface_mask="glacier_rgi7", geotile_width=2); dims=:height); dims=:height)

ids = intersect(val(dims(area6, :geotile)), val(dims(area7, :geotile)))
a6 = [area6[geotile=At(id)] for id in ids]
a7 = [area7[geotile=At(id)] for id in ids]
both = (a6 .> 0) .& (a7 .> 0)
ratio = a7[both] ./ a6[both]
weight = a6[both]

function weighted_quantile(x, w, p)
    order = sortperm(x)
    cw = cumsum(w[order]) ./ sum(w)
    return x[order][searchsortedfirst(cw, p)]
end

println("RGI7 / RGI6 glacier area per geotile ($(count(both)) geotiles with ice in both)")
@printf("  per-geotile quantiles    5 %%: %.2f   50 %%: %.2f   95 %%: %.2f\n",
        (quantile(ratio, p) for p in (0.05, 0.5, 0.95))...)
@printf("  area-weighted quantiles  5 %%: %.2f   50 %%: %.2f   95 %%: %.2f\n",
        (weighted_quantile(ratio, weight, p) for p in (0.05, 0.5, 0.95))...)
@printf("  ice area in geotiles where the ratio differs from 1 by more than 10 %%: %.1f %%\n",
        100 * sum(weight[abs.(ratio .- 1) .> 0.1]) / sum(weight))
@printf("  geotiles with ice in only one inventory: RGI6 %d, RGI7 %d\n",
        count((a6 .> 0) .& (a7 .== 0)), count((a7 .> 0) .& (a6 .== 0)))

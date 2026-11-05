############################################################################################
#                           Gardner2025_geotiles_rates_myr.gpkg                          #
#                   m.i.e to m.w.e conversion = 910 kg/m3 / 1000 kg/m3 = 0.91              #
############################################################################################
Variables:
- id: geotile id
- area_km2: glacier area (RGI7)
- group: geotile grouping id
- pscale: optimal precipitation scaling [unitless]
- mscale: optimal melt scaling [unitless]
- rgiid: RGI region id (RGI6)
- discharge_trend: solid ice discharge through flux gate located near grounding line [m.i.e./yr]
- dh_synthesis_trend: trend in height anomaly from synthesis of satellite data [m/yr]
- dh_synthesis_amplitude: seasonal amplitude in height anomaly from synthesis of satellite data [m]
- dh_synthesis_phase: day of seasonal maximum in height anomaly from synthesis of satellite data [day of year]
- runoff_trend: trend in glacier meltwater runoff (excludes rain) [m.i.e./yr]
- runoff_amplitude: seasonal amplitude in glacier meltwater runoff (excludes rain) [m.i.e.]
- runoff_phase: day of seasonal maximum in glacier meltwater runoff (excludes rain) [day of year]
- fac_trend: trend in firn air content [m/yr]
- fac_amplitude: seasonal amplitude in firn air content [m]
- fac_phase: day of seasonal maximum in firn air content [day of year]
- smb_trend: trend in surface mass balance [m.i.e./yr]
- smb_amplitude: seasonal amplitude in surface mass balance [m.i.e.]
- smb_phase: day of seasonal maximum in surface mass balance [day of year]
- acc_trend: trend in accumulation [m.i.e./yr]
- acc_amplitude: seasonal amplitude in accumulation [m.i.e.]
- acc_phase: day of seasonal maximum in accumulation [day of year]
- melt_trend: trend in surface melt [m.i.e./yr]
- melt_amplitude: seasonal amplitude in surface melt [m.i.e.]
- melt_phase: day of seasonal maximum in surface melt [day of year]
- ec_trend: trend in evaporation/condensation [m.i.e./yr]
- ec_amplitude: seasonal amplitude in evaporation/condensation [m.i.e.]
- ec_phase: day of seasonal maximum in evaporation/condensation [day of year]
- refreeze_trend: trend in refreezing [m.i.e./yr]
- refreeze_amplitude: seasonal amplitude in refreezing [m.i.e.]
- refreeze_phase: day of seasonal maximum in refreezing [day of year]
- dh_trend: trend in height change [m/yr]
- dh_amplitude: seasonal amplitude in height change [m]
- dh_phase: day of seasonal maximum in height change [day of year]
- dm_trend: trend in mass change [m.w.e./yr]
- dm_amplitude: seasonal amplitude in mass change [m.w.e.]
- dm_phase: day of seasonal maximum in mass change [day of year]
- dm_synthesis_trend: trend in mass change from synthesis of satellite data corrected for changes in fac [m.w.e./yr]
- dm_synthesis_amplitude: seasonal amplitude in mass change from synthesis of satellite data corrected for changes in fac [m.w.e.]
- dm_synthesis_phase: day of seasonal maximum in synthesis of satellite data corrected for changes in fac [day of year]
############################################################################################


############################################################################################
#                         Gardner2025_regional_timseries.nc                                #
############################################################################################
Dimensions:
- rgi: RGI region id (RGI6)
- date: center date of monthly data

Variables (each variable has a corresponding _error variable representing the
95% confidence interval of the formal error, an additional minimum fractional error of 20% 
should be applied to all component rates):
- fac: firn air content anomaly [km^3]
- dv_synthesis: volume anomaly from synthesis of satellite data [km^3]
- dv: volume anomaly from optimal model fit [km^3]
- dm_synthesis: mass anomaly from firn corrected synthesis of satellite data [Gt]
- dm: mass anomaly from optimal model fit [Gt]
- dm_grace: mass anomaly from GRACE [Gt]
- dm_glambie: mass anomaly reported by GlaMBIE Team, 2025 [Gt]
- smb: surface mass balance anomaly from optimal model fit [Gt]
- acc: accumulation anomaly [Gt]
- melt: melt anomaly [Gt]
- refreeze: refreeze anomaly [Gt]
- runoff: meltwater runoff anomaly (excludes rain) [Gt]
- ec: evaporation/condensation/sublimation/deposition anomaly from optimal model fit [Gt]
- discharge: solid ice discharge anomaly [Gt]

Note: Rain variables have been excluded because snow is scaled in GEMB runs to account for
avalanching and wind redistribution, which unrealistically scales rain by the same factor.
############################################################################################


############################################################################################
#                          Gardner2025_endorheic_fraction.csv                              #
############################################################################################
Variables:
- rgi: RGI region id (RGI6)
- region_name: RGI region name (RGI6)
- dm_endorheic_fraction: fraction of mass change that terminates in a closed basin (i.e. endorheic)
- runoff_endorheic_fraction: fraction of meltwater runoff that terminates in a closed basin (i.e. endorheic, excludes rain)
############################################################################################


############################################################################################
#                          Gardner2025_regional_results.csv                                #
############################################################################################
All trends and accelerations are computed over the period March 2000 to December 2024.
Values are reported as best estimate ± 95% confidence interval. A minimum error of 20% is 
applied on all component rates excluding rates of volume and mass change that are well 
represented by the formal errors.

Variables:
- rgi: RGI region id (RGI7)
- region_name: RGI region name (RGI7)
- dv_synthesis_[km³/yr]: trend in volume anomaly from synthesis of satellite data ± 95% confidence interval [km³/yr]
- dv_synthesis_[km³/yr²]: acceleration in volume anomaly from synthesis of satellite data ± 95% confidence interval [km³/yr²]
- dv_synthesis_amplitude_[km³]: seasonal amplitude in volume anomaly from synthesis of satellite data ± 95% confidence interval [km³]
- dv_[km³/yr]: trend in volume anomaly from optimal model fit ± 95% confidence interval [km³/yr]
- dv_[km³/yr²]: acceleration in volume anomaly from optimal model fit ± 95% confidence interval [km³/yr²]
- dv_amplitude_[km³]: seasonal amplitude in volume anomaly from optimal model fit ± 95% confidence interval [km³]
- dm_[Gt/yr]: trend in mass anomaly ± 95% confidence interval [Gt/yr]
- dm_[Gt/yr²]: acceleration in mass anomaly ± 95% confidence interval [Gt/yr²]
- acc_[Gt/yr]: trend in accumulation anomaly ± 95% confidence interval [Gt/yr]
- acc_[Gt/yr²]: acceleration in accumulation anomaly ± 95% confidence interval [Gt/yr²]
- runoff_[Gt/yr]: trend in meltwater runoff anomaly (excludes rain) ± 95% confidence interval [Gt/yr]
- runoff_[Gt/yr²]: acceleration in meltwater runoff anomaly (excludes rain) ± 95% confidence interval [Gt/yr²]
- melt_[Gt/yr]: trend in melt anomaly ± 95% confidence interval [Gt/yr]
- melt_[Gt/yr²]: acceleration in melt anomaly ± 95% confidence interval [Gt/yr²]
- refreeze_[Gt/yr]: trend in refreeze anomaly ± 95% confidence interval [Gt/yr]
- refreeze_[Gt/yr²]: acceleration in refreeze anomaly ± 95% confidence interval [Gt/yr²]
- ec_[Gt/yr]: trend in evaporation/condensation/sublimation/deposition anomaly ± 95% confidence interval [Gt/yr]
- ec_[Gt/yr²]: acceleration in evaporation/condensation/sublimation/deposition anomaly ± 95% confidence interval [Gt/yr²]
- discharge_[Gt/yr]: trend in solid ice discharge through flux gate located near grounding line ± 95% confidence interval [Gt/yr]
- smb_[Gt/yr]: trend in surface mass balance anomaly ± 95% confidence interval [Gt/yr]
- smb_[Gt/yr²]: acceleration in surface mass balance anomaly ± 95% confidence interval [Gt/yr²]
- fac_[km³/yr]: trend in firn air content anomaly ± 95% confidence interval [km³/yr]
- fac_[km³/yr²]: acceleration in firn air content anomaly ± 95% confidence interval [km³/yr²]
- gsi: glacier sustainability index [unitless]
- area_km2: glacier area (RGI7) [km²]
############################################################################################


############################################################################################
#                          Gardner2025_gmax_buffer_population.nc                           #
############################################################################################
Dimensions:
- country: country name
- gmax: minimum gmax threshold [%] (values: 25, 50)
- runoff: minimum climatological maximum monthly glacier meltwater runoff threshold (excludes rain) [m³/s] (values: 0, 1, 2.5, 5, 10, 100, 1000)
- buffer: maximum buffer distance from river [m] (values: 1000, 5000, 10000, 15000, 20000, 25000, 30000, 40000, 50000)

Variables:
- population: population living within a given buffer distance of a river with a gmax and runoff threshold
############################################################################################


############################################################################################
#                          Gardner2025_glacier_summary_gmax.gpkg                           #
############################################################################################
Variables:
- COMID: MERIT Hydro river reach unique ID
- ocean_terminating: river reach flows into ocean [true/false]
- continent: continent of river reach
- country: country of river reach
- lengthkm: length of river reach [km]
- runoff_max_avg: climatological maximum monthly glacier meltwater runoff (excludes rain) [m³/s]
- runoff_max_month: month of climatological maximum monthly glacier meltwater runoff (excludes rain)
- gmax_max: maximum gmax between March 1, 2000 and October 1, 2024
- gmax_min: minimum gmax between March 1, 2000 and October 1, 2024
- gmax_range: range in gmax between March 1, 2000 and October 1, 2024
- gmax_avg: average gmax between March 1, 2000 and October 1, 2024
- gmax_month: month of climatological maximum gmax
############################################################################################

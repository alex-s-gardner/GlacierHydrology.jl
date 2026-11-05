module GlobalGlacierAnalysis

   module MyUnits
      using Unitful; @unit Gt "Gt" M 1u"Pg" false
   end

   using Unitful
   Unitful.register(MyUnits)
   using GlobalGlacierAnalysis.MyUnits
   
   # import geographic packages
   using Proj
   using Rasters
   using Shapefile
   using LibGEOS
   using GeoInterface

   # file IO
   using FileIO
   using MAT
   using CSV
   using NCDatasets
   using JLD2
   

   # standard packages
   using DataFrames
   using Dates
   using Statistics

   # plotting packages
   using ColorSchemes
   using CairoMakie
   
   # other packages
   using Extents
   using Interpolations
   using Printf
   using StatsBase
   import TOML
   using LsqFit
   using Distances
   using DataInterpolations
   using FlexiJoins
   using LinearAlgebra
   using Rotations
   using CoordinateTransformations
   using SortTileRecursiveTree
   using AbstractTrees
   using ProgressMeter
   using Random
   using GeometryBasics

   using BinStatistics

   # import packages
   import DimensionalData as DD
   import GeometryOps as GO
   import GeoInterface as GI
   import GeoFormatTypes as GFT
   import GeometryOpsCore as GOC
   import SortTileRecursiveTree

   # Glacier volume change comes from GlacierVolumeChange.jl: GGA reads each ensemble member's volume
   # product (`GVC.read_volume`) and uses its geotile, time and fitting helpers listed here.
   import GlacierVolumeChange as GVC
   using GlacierVolumeChange:
      _geotile_area_km2, _publication_figure, allfiles, binned_filled_fileparts, centroid,
      decimalyear, decimalyear2datetime, extent2rectangle, file_is_current, geotile_extent,
      geotile_hypsometry, geotiles_golden_test, geotiles_mask_hyps,
      geotiles_mutually_exclusive_rgi!, intersectindices, meters2lonlat_distance, model1_trend,
      offset_trend, offset_trend_acceleration_seasonal2, offset_trend_p, offset_trend_seasonal2,
      p1, p3, p_offset_trend_seasonal, p_trend, project_height_bins, seasonal_peak_fraction, slope,
      validrange, within

   # this is where the local paths are defined
   include("local_paths.jl")

   # add utilities
   include("utilities_project.jl")
   include("utilities_gemb.jl")
   include("utilities_gemb_tiles.jl")
   include("utilities_wgms.jl")
   include("utilities_gemb_physics.jl")
   include("utilities_routing.jl")
   include("utilities_synthesis.jl")
   include("utilities_postprocessing.jl")
   include("utilities.jl")
   include("utilities_plotting.jl")
   include("utilities_readers.jl")

end

# config template with config defaults
# not a separate module, only sourced from palm_dynamic_config

# PALM case, domain, and configuration parameters
# (only an example, needs to be rewriten in user config)
domain = ''
resolution = ''
scenario = ''
nested_domain = False

# file name of output dynamic driver ("" means the standard name)
dynamic_driver_file = ""

# import grid parameters for dynamic driver from static driver
grid_from_static = True
# file name of static driver ("" means the standard name)
static_driver_file = ""
# reference coordinate system of PALM simulation32633
proj_palm = "EPSG:32633"

# projection lon-lat4326
proj_wgs84 = 'EPSG:4326'

# vertical grid
# layer height (dz = 0.0 means dz is assigned from dx)
dz = 0.0
# we need description of the PALM vertical structure for dynamic driver
nz = 80  # z in grids
dz_stretch_level = 1200.0 # in meters
dz_stretch_factor = 1.02
dz_max = 60.0

# time origin and extent of the simulation (format YYYY-MM-DD hh:mm:ss)
origin_time = ""
simulation_hours = 2

# WRF related configurations
wrf_hybrid_levs = True
# Smoothing of PALM terrain for WRF vertical interpolation to avoid sharp
# horizontal gradients. None = off
vinterp_terrain_smoothing = None

# wrf-chem input files path and default file mask
wrf_dir_name = ""  # "" means that standard path will be calculated in the init
wrf_file_mask =  "wrfout_d01_*"
# process radiation from wrf, true to include
radiation_from_wrf = False
wrf_rad_file_mask = "auxhist6_*"
# smoothing distance for radiation
radiation_smoothing_distance = 10000.0

# PALM species names - version for phstatp2 mechanism
species_names = ['NO', 'NO2', 'O3', 'PM10', 'PM25']

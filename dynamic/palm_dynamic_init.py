#!/usr/bin/python3
# -*- coding: utf-8 -*-

#------------------------------------------------------------------------------#
#
# Scripts for processing of WRF-CHEM files to PALM dynamic driver.
#
#------------------------------------------------------------------------------#
'''
This file provides initialization of the basic variables and structures base on
the default and user config values. It can be adjusted to the needs of particular
user system or structure of the data storage.
'''

import os
from pathlib import Path
import numpy as np

# paths of directories
dir_base = os.path.abspath(Path(dir_scripts).parent.parent)
dir_in = os.path.join(dir_base, domain)
if scenario == '':
    dir_in_scen = dir_in
else:
    dir_in_scen = os.path.join(dir_in,scenario)
dir_out = os.path.join(dir_base, domain)

# file names of PALM PIDS drivers
# extension of file name for scenario
scenario_ext = ("_"+scenario if scenario != "" else "")
# static driver netcdf file name
if static_driver_file == "":
    static_driver_file = os.path.join(dir_out, domain+"_static_driver"+"_d"+resolution+scenario_ext+".nc")
# dynamic driver netcdf file name
if dynamic_driver_file == "":
    dynamic_driver_file = os.path.join(dir_out, domain+"_dynamic_driver"+"_d"+resolution+scenario_ext+".nc")

# parameters of dynamic driver
# minimal number of free surface canopy layers above top of terrain with building and plant canopy
nscl_free = 3
# default path of wrf files in case it is not set in user config
if wrf_dir_name == '':
    wrf_dir_name = os.path.join(dir_in, 'wrf')
    interp_dir_name = wrf_dir_name

# Settings for geostrophic wind
gw_gfs_margin_deg = 5. #smoothing area in degrees lat/lon
gw_wrf_margin_km = 10. #smoothing area in km
#gw_alpha = .143 #GW vertical interpolation by power law
gw_alpha = 1. #ignore wind power law, interpolate linearly

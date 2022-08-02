**Scripts for processing of WRF-CHEM files to PALM dynamic driver.**

Version: v.1.0

The scripts are based on the [wrf-CAMx interface](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/iofiles/wrf_interface)

**Usage:**
- `palm_dynamic -c <config_name> [-w]`
- The optional parameter -w allows to skip horizontal and vertical interpolation in case it has already been done.
- Example: `python3 palm_dynamic.py -c augsburg_validation_summer_10`


**Script requirement:**
- The name of the case configuration on the command line.
- The corresponding configuration file `config_name.conf` in the `configuration` subdirectory. 
- The values which agree with defaults do not need to be set. `palm_dynamic_init.py` contains setting and calculation of standard  initialization values for particular system and can be adjusted.

**Required modules:**
- [numpy](https://pypi.org/project/numpy)
- [scipy](https://pypi.org/project/scipy)
- [pyproj](https://pypi.org/project/pyproj)
- [netCDF4](https://pypi.org/project/netCDF4)
- [metpy](https://unidata.github.io/MetPy)

In the current version, the only supported projection in WRF-CHEM is Lambertiam conformal conic, which is WRF default and recommended projection for
mid-latitudes.

The scripts support both variants of WRF-CHEM vertical levels - the sigma levels (default until WRF version 3.*) and the hybrid levels (default since WRF 4.*). However, it is necessary to correctly configure this option via the setting `wrf_hybrid_levs = True/False`.

# Configuration parameters
Description of the particular configuration options are (defaults are in parenthesis):

## Domain and case related configurations
    - domain              Name of the simulation case ("").
    - resolution          Name of the particular domain resolution scenario ("").
    - scenario            Name of the individual scenario in the case ("").
    - nested_domain       False indicates parent and True nested domain. (False).

    - dynamic_driver_file File name of output dynamic driver ("").
    - grid_from_static    True - import from static driver, False - given in config (True).
    - static_driver_file  File name of the static driver in case of grid_from_static ("").
    - proj_palm           Reference coordinate system of PALM simulation ("EPSG:32633").
    - proj_wgs84          Reference coordinate system of lon-lat projection ("EPSG:4326").
    
    - dz                  Height of the PALM vertical grid layer (0.0). The default value dz = 0.0 means dz is assigned from dx.
    - nz                  Number of vertical layers of PALM domain (200).
    - dz_stretch_level    Height in meters from which stretching of vertical levels starts in PALM (5000.0).
    - dz_stretch_factor   Coefficient of the stretching of the vertical layers in PALM (1.0).
    - dz_max              Max height of the stretched vertical layers (100.0).
    
    - origin_time         Origin time of the PALM simulation in the format YYYY-MM-DD hh:mm:ss (""). Default reads from static driver.
    - simulation_hours    Extent of the simulation in hours.

## WRF-CHEM related configurations
    - wrf_dir_name        File path of the wrf-chem input files (""). 
    - wrf_file_mask       File mask of the wrf-chem input files  ("wrfout_*.e000")
    - wrf_hybrid_levs     True means hybrid levels in WRF files, False means sigma levels (True).
    - vinterp_terrain_smoothing
        - the standard deviation for Gaussian kernel of smoothing method of the PALM terrain for WRF vertical interpolation to avoid sharp horizontal gradients. Value None disables the smoothing. (None)
    - interp_dir_name     File path to interpolated files
    - wrfchem_spec        Wrf-chem chemical species to be included in dynamic driver, `list of gaseous species: no, no2, no3, pm10, PM2_5_DRY, o3, co, hno3, ho, h2o2, nh3`

    - aerosol_wrfchem     True means aerosols are included, (False)
    - wrfchem_bin_limits  Wrf-chem aerosol siz bins `[3.9e-8, 1.56e-7, 6.25e-7, 2.5e-6, 1.0e-5]`
    - listspec            PALM aerosol species, only aerosols, options listspec = ['SO4', 'OC', 'BC', 'DU', 'SS', 'NH', 'NO']
    - nbin                SALSA parameter, # size bins in subrange `[1,7]`
    - reglim              SALSA parameter, subrange limits `[3.9e-8, 5.0e-8, 2.5e-6]`
    - nf2a                SALSA parameter, insoluble fraction, currently only soluble supported in PALM (1.0)
                 
    - radiation_from_wrf  enable or disable processing of radiation from WRF files (True).
    - wrf_rad_file_mask   file mask of the wrf radiation input files ("auxhist6_*"). The default setting reads radiation from WRF auxiliary history files. This setting allows to use finer time step for WRF radiation outputs than for other values.
    - radiation_smoothing_distance smoothing distance for radiation values in m (10000.0).

## Horizontal parameters of the PALM domain which have to be set in case of grid_from_static = False
* `nx, ny`              Number of horizontal grids of the domain in x and y directions.
* `dx, dy`              Grid cell size of the domain in x and y directions.
* `origin_x, origin_y`  Origin x and y of the domain.
* `origin_z`            Origin of the domain in the vertical direction.

###  Major Changes from wrf + CAMx scripts made to wrf_chem_for_palm
- Proj future warning resolved
- Chemical species are read from wrf-chem files and interpolated at the same time as the dynamic variables
- Interpolated files are saved to a different directory than the wrf-chem data files
- A variety of chemical species can be inlcuded.
- Aerosols can be included and are weighted based on specified aerosol size bins
- Base state temperature corrected to 290. for WRFChem
- Changes to palm_dynamic_output were made to include the various aerosol and chemical species

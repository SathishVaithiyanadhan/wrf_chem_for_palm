#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------#
#
# Scripts for processing of WRF-CHEM files to PALM dynamic driver.
#
#------------------------------------------------------------------------------#
'''
This file creates and writes the dynamic driver netcdf file based on preprepared
transformed and interpolated wrf-chem files.
'''
import os
import time
import numpy as np
import netCDF4
from palm_wrf_utils import palm_wrf_gw
from palm_dynamic_config import *
import palm_dynamic_aerosol

def palm_dynamic_output(wrf_files, interp_files, dynamic_driver_file, times_sec,
                        dimensions, z_levels, z_levels_stag, ztop,
                        z_soil_levels, dx, dy, lon_center, lat_center,
                        rad_times_proc, rad_values_proc, sma, nested_domain):
    print('\nProcessing interpolated files to dynamic driver')
    # dimension of the time coordinate
    dimtimes = len(times_sec)
    # other coordinates
    dimnames = ['z', 'zw', 'zsoil', 'x','xu', 'y', 'yv']                # z  height agl in m, zw staggered, zsoil 4 lev from wrf 
    dimsize_names = ['zdim' , 'zwdim', 'zsoildim', 'xdim', 'xudim', 'ydim', 'yvdim']    # palm: zw = z - 1
    x = np.arange(dx, dimensions['xdim']*dx+dx, dx)
    y = np.arange(dy, dimensions['ydim']*dy+dy, dy)
    # fill values
    fillvalue_float = float(-9999.0)
    try:
        os.remove(dynamic_driver_file)
    except:
        pass
    # boundaries
    boundary = ['left','right','south', 'north','top']
    # atmospheric variables
    atmos_var = ['pt','qv','u','v','w','soil_m','soil_t']
    if len(wrfchem_spec) > 0:
        dynam_chem_variables = atmos_var + wrfchem_spec
    else:
        dynam_chem_variables = atmos_var

    # prepare influx/outflux area sizes
    zstag_all = np.r_[0., z_levels_stag, ztop]
    zwidths = zstag_all[1:] - zstag_all[:-1]
    areas_xb = np.zeros((len(z_levels), 1))
    areas_xb[:,0] = zwidths * dy
    areas_yb = np.zeros((len(z_levels), 1))
    areas_yb[:,0] = zwidths * dx
    areas_zb = dx*dy
    area_boundaries = (areas_xb.sum()*dimensions['ydim']*2
            + areas_yb.sum()*dimensions['xdim']*2
            + areas_zb*dimensions['xdim']*dimensions['ydim'])

    # create netcdf output file
    print('Creating Dynamic Driver:', dynamic_driver_file)
    outfile = netCDF4.Dataset(dynamic_driver_file, "w", format="NETCDF4" )
    # create dimensions (time and coordinates)
    outfile.createDimension('time', dimtimes)
    for _dim in zip(dimnames, dimsize_names):
        outfile.createDimension(_dim[0], dimensions[_dim[1]])
    # create variables (time and coordinates)
    _val_times            = outfile.createVariable('time',"f4", ("time"))
    _val_times[:]         = times_sec[:]

    _val_z_levels         = outfile.createVariable('z',"f4", ("z"))
    _val_z_levels[:]      = z_levels[:]

    _val_z_levels_stag    = outfile.createVariable('zw',"f4", ("zw"))
    _val_z_levels_stag[:] = z_levels_stag[:]

    _val_z_soil_levels    = outfile.createVariable('zsoil',"f4", ("zsoil"))
    _val_z_soil_levels[:] = z_soil_levels[:]

    _val_y             = outfile.createVariable('y',"f4", ("y"))
    _val_y[:]          = y[:]

    _val_x             = outfile.createVariable('x',"f4", ("x"))
    _val_x[:]          = x[:]
    #---------------------------------------------------------------------------
    # include aerosols, add dimensions and variables
    if aerosol_wrfchem:
        print('\tCreating aerosol variables in dynamic driver')
        # create dimensions
        outfile.createDimension('Dmid', sum(nbin))

        # include for gaseous aerosols from chemistry
        idx_len = len(listspec)
        outfile.createDimension('composition_index', idx_len)
        outfile.createDimension('max_string_length', 25)

        # variables
        val_dmid = outfile.createVariable('Dmid', "f4",('Dmid',))
        val_dmid.setncattr('units', 'm')
        _dmid, _bim = palm_dynamic_aerosol.define_bins(nbin, reglim)
        _dmid = np.array(_dmid)
        val_dmid[:] = _dmid[:]

        val_idx = outfile.createVariable('composition_index', "i4", ('composition_index',))
        for n in range(0, idx_len):
            val_idx[n] = n+1

        # 2D vertical profile variables
        _val_composition_name    = outfile.createVariable('composition_name',"S1", ("composition_index", "max_string_length"))
        _val_composition_name.setncattr('long_name', "aerosol composition name")
        spec_name = listspec
        composition_name         = np.array([spec_name],dtype = 'S25')
        _val_composition_name[:] = netCDF4.stringtochar(composition_name)

        # upwind values for aerosol concen and mass frac
        aero_massfrac_a          = palm_dynamic_aerosol.aerosolProfileMass(interp_files, listspec)
        _val_aerosol_mass_a      = outfile.createVariable('init_atmosphere_mass_fracs_a', "f4", ("z", "composition_index"),fill_value=fillvalue_float)
        _val_aerosol_mass_a.setncattr('long_name',"initial mass fraction profile: a bins")
        _val_aerosol_mass_a[:,:] = aero_massfrac_a

        _val_aerosol_mass_b      = outfile.createVariable('init_atmosphere_mass_fracs_b', "f4", ("z", "composition_index"),fill_value=fillvalue_float)
        _val_aerosol_mass_b.setncattr('long_name',"initial mass fraction profile: b bins")
        if nf2a == 1.0:
            _val_aerosol_mass_b[:,:] = 0
        elif nf2a <= 1.0:
            _val_aerosol_mass_b[:,:] = 1 - aero_massfrac_a

        # aerosol concen.#
        aero_concen = palm_dynamic_aerosol.aerosolProfileConc(interp_files, nbin, reglim, wrfchem_bin_limits)
        _val_aerosol_con = outfile.createVariable('init_atmosphere_aerosol', "f4", ("z", "Dmid"), fill_value = fillvalue_float)
        _val_aerosol_con.setncattr('lod', 1)
        _val_aerosol_con.setncattr('long_name',"initial vertical profile of aerosol concentration")
        _val_aerosol_con.setncattr('units',"#/m3")
        _val_aerosol_con[:,:] = aero_concen

        # time dependent variables - aerosols
        for side in boundary:
            if (side == 'left' or side == 'right'):
                _dim1 = "z"
                _dim2 = "y"
            elif (side == 'south' or side == 'north'):
                _dim1 = "z"
                _dim2 = "x"
            elif side =='top':
                _dim1 = "y"
                _dim2 = "x"
            # mass fractions
            _var_a = outfile.createVariable('ls_forcing_'+ side + '_mass_fracs_a', "f4", ("time", _dim1, _dim2, "composition_index"),
                    fill_value=fillvalue_float)
            _var_a.setncattr('lod', 1)
            _var_a.setncattr('long_name', "boundary conditions of mass fraction profile: a bins")
            _var_a[:,:,:,:] = palm_dynamic_aerosol.aerosolMassWrfchemBoundary(dimensions, _dim1, _dim2, interp_files, listspec, side)

            _var_b = outfile.createVariable('ls_forcing_'+ side + '_mass_fracs_b', "f4", ("time", _dim1, _dim2, "composition_index"),
                    fill_value=fillvalue_float)
            _var_b.setncattr('lod', 1)
            _var_b.setncattr('long_name', "boundary conditions of mass fraction profile: b bins")
            if (nf2a >= 1.0):
                _var_b[:,:,:,:] = 0
            elif (nf2a < 1.0):
                _var_b[:,:,:,:] = 1-_var_a

            # aerosol concen.#
            _var_con = outfile.createVariable('ls_forcing_'+ side + '_aerosol', "f4", ("time", _dim1, _dim2, "Dmid"),
                    fill_value=fillvalue_float)
            _var_con.setncattr('lod', 1)
            _var_con.setncattr('long_name',"boundary condition of aerosol concentration")
            _var_con[:,:,:,:] = palm_dynamic_aerosol.aerosolConWrfchemBoundary(dimensions, _dim1, _dim2, interp_files, side, nbin, reglim, wrfchem_bin_limits)

    #---------------------------------------------------------------------------
    # create dynamical & chemical variables in output
    def add_interpDim(dynam_chem_variables):
        print('\tCreating dynamical variables in dynamic driver')
        # surface pressure
        _val_surface_forcing_surface_pressure = outfile.createVariable('surface_forcing_surface_pressure', "f4",("time"))
        # geostrophic wind
        _val_ls_forcing_ug = outfile.createVariable('ls_forcing_ug', "f4", ("time", "z"),fill_value=fillvalue_float)
        _val_ls_forcing_vg = outfile.createVariable('ls_forcing_vg', "f4", ("time", "z"),fill_value=fillvalue_float)
        # create variables other dynamic & chemical variables
        for var in dynam_chem_variables:
            if (var == 'pt' or var == 'qv'):
                _val_init_var = outfile.createVariable('init_atmosphere_'+ var, "f4", ("z", "y", "x"),fill_value=fillvalue_float)
                _val_init_var.setncattr('lod', 2)
            elif (var == 'soil_m' or var == 'soil_t'):
                _val_init_var = outfile.createVariable('init_'+ var,"f4", ("zsoil", "y", "x"),fill_value=fillvalue_float)
                _val_init_var.setncattr('lod', 2)
            elif var == 'u':
                _val_init_var = outfile.createVariable('init_atmosphere_'+ var, "f4", ("z", "y", "xu"),fill_value=fillvalue_float)
                _val_init_var.setncattr('lod', 2)
            elif var == 'v':
                _val_init_var = outfile.createVariable('init_atmosphere_'+ var, "f4", ("z", "yv", "x"),fill_value=fillvalue_float)
                _val_init_var.setncattr('lod', 2)
            elif var == 'w':
                _val_init_var = outfile.createVariable('init_atmosphere_'+ var, "f4", ("zw", "y", "x"),fill_value=fillvalue_float)
                _val_init_var.setncattr('lod', 2)
            elif var in wrfchem_spec:
                _val_init_var = outfile.createVariable('init_atmosphere_'+ var.upper(), "f4", ("z",),fill_value=fillvalue_float)
                if (var == 'PM10' or var == 'PM2_5_DRY'):
                    _val_init_var.setncattr('lod', 1)
                    _val_init_var.setncattr('units','kg/m3')
                else:
                    _val_init_var.setncattr('lod', 1)
                    _val_init_var.setncattr('units','ppm')

            # create time dependent variables        
            if not nested_domain:
                if (var == 'soil_m' or var == 'soil_t'):
                    continue
                else:
                    for side in boundary:
                        if (side == 'left' or side == 'right'):
                            if var == 'u':
                                _val_ls_forcing = outfile.createVariable('ls_forcing_'+ side +'_'+ var, "f4", ("time", "z", "y"),
                                        fill_value=fillvalue_float)
                            elif var == 'v':
                                _val_ls_forcing = outfile.createVariable('ls_forcing_'+ side +'_'+ var, "f4", ("time", "z", "yv"),
                                        fill_value=fillvalue_float)
                            elif var == 'w':
                                _val_ls_forcing = outfile.createVariable('ls_forcing_'+ side +'_'+ var, "f4", ("time", "zw", "y"),
                                        fill_value=fillvalue_float)
                            elif var in wrfchem_spec:
                                _val_ls_forcing = outfile.createVariable('ls_forcing_'+ side +'_'+ var.upper(), "f4", ("time", "z", "y"),
                                        fill_value=fillvalue_float)
                            else:
                                _val_ls_forcing = outfile.createVariable('ls_forcing_'+ side +'_'+ var, "f4", ("time", "z", "y"),
                                        fill_value=fillvalue_float)

                        elif (side == 'south' or side == 'north'):
                            if var == 'u':
                                _val_ls_forcing = outfile.createVariable('ls_forcing_'+ side +'_'+ var, "f4", ("time", "z", "xu"),
                                        fill_value=fillvalue_float)
                            elif var == 'v':
                                _val_ls_forcing = outfile.createVariable('ls_forcing_'+ side +'_'+ var, "f4", ("time", "z", "x"),
                                        fill_value=fillvalue_float)
                            elif var == 'w':
                                _val_ls_forcing = outfile.createVariable('ls_forcing_'+ side +'_'+ var, "f4", ("time", "zw", "x"),
                                        fill_value=fillvalue_float)
                            elif var in wrfchem_spec:
                                _val_ls_forcing = outfile.createVariable('ls_forcing_'+ side +'_'+ var.upper(), "f4", ("time", "z", "x"),
                                        fill_value=fillvalue_float)
                            else:
                                _val_ls_forcing = outfile.createVariable('ls_forcing_'+ side +'_'+ var, "f4", ("time", "z", "x"),
                                        fill_value=fillvalue_float)
                        else:
                            if var == 'u':
                                _val_ls_forcing = outfile.createVariable('ls_forcing_'+ side +'_'+ var, "f4", ("time", "y", "xu"),
                                        fill_value=fillvalue_float)
                            elif var == 'v':
                                _val_ls_forcing = outfile.createVariable('ls_forcing_'+ side +'_'+ var, "f4", ("time", "yv", "x"),
                                        fill_value=fillvalue_float)
                            elif var == 'w':
                                _val_ls_forcing = outfile.createVariable('ls_forcing_'+ side +'_'+ var, "f4", ("time", "y", "x"),
                                        fill_value=fillvalue_float)
                            elif var in wrfchem_spec:
                                _val_ls_forcing = outfile.createVariable('ls_forcing_'+ side +'_'+ var.upper(), "f4", ("time", "y", "x"),
                                        fill_value=fillvalue_float)
                            else:
                                _val_ls_forcing = outfile.createVariable('ls_forcing_'+ side +'_'+ var, "f4", ("time", "y", "x"),
                                        fill_value=fillvalue_float)
                        # atrributes
                        _val_ls_forcing.setncattr('lod', 2)
                        if var not in atmos_var:
                            if (var == 'PM10' or var == 'PM2_5_DRY'):
                                _val_ls_forcing.setncattr('units', 'kg/m3')
                            else:
                                _val_ls_forcing.setncattr('units', 'ppm')
    # function - create dynamical & chemical variables
    add_interpDim(dynam_chem_variables)

    outfile.close()
    
    # read interpolated files and write values for dynamical & chemical variables
    def add_interpValues(dynam_chem_variables):
        print('\tCreating initialisation variables in dynamic driver')
        infile = netCDF4.Dataset(interp_files[0], "r", format="NETCDF4")
        outfile = netCDF4.Dataset(dynamic_driver_file, "r+", format="NETCDF4")
        # initialization variables
        for var in dynam_chem_variables:
            if var == 'soil_m':
                init_var = infile.variables['init_'+var]
                _val_init_var = outfile.variables['init_'+var]
                for k in range(0,_val_init_var.shape[0]):
                    _val_init_var[k, :, :] = init_var[0, k, :, :] * sma[:, :]

            elif var == 'soil_t':
                init_var = infile.variables['init_'+var]
                _val_init_var = outfile.variables['init_'+var]
                _val_init_var[:, :, :] = init_var[0, :, :, :]

            elif (var == 'pt' or var == 'qv' or var == 'w' or var == 'u' or var == 'v'):
                init_var = infile.variables['init_atmosphere_'+ var]
                _val_init_var = outfile.variables['init_atmosphere_'+ var]
                if var == 'u':
                    _val_init_var[:, :, :] = init_var[0, :, :, 1:]
                elif var == 'v':
                    _val_init_var[:, :, :] = init_var[0, :, 1:, :]
                elif (var == 'pt' or var == 'qv' or var == 'w'):
                    _val_init_var[:, :, :] = init_var[0, :, :, :]

            elif var in wrfchem_spec:
                if var == 'PM10':
                    init_var = infile.variables['init_atmosphere_'+ var.lower()]
                    _val_init_var = outfile.variables['init_atmosphere_'+ var]
                    _val_init_var[:] = init_var[0,:,:,:].mean(axis=(1,2))
                else:
                    init_var = infile.variables['init_atmosphere_'+ var]
                    _val_init_var = outfile.variables['init_atmosphere_'+ var.upper()]
                    _val_init_var[:] = init_var[0,:,:,:].mean(axis=(1,2))

        infile.close()
        outfile.close()

        # time dependent variables - dynamical & chemical variables
        if not nested_domain:
            print('\tCreating time dependent variables in dynamic driver')
            outfile = netCDF4.Dataset(dynamic_driver_file, "r+", format="NETCDF4")
            for ts in range(0, len(interp_files)):
                # geostrophic wind
                #print('Open wrf file: '+wrf_files[ts])
                nc_wrf = netCDF4.Dataset(wrf_files[ts], 'r')
                ug, vg = palm_wrf_gw(nc_wrf, lon_center, lat_center, z_levels)
                _val_ls_forcing_ug = outfile.variables['ls_forcing_ug']
                _val_ls_forcing_vg = outfile.variables['ls_forcing_vg']
                _val_ls_forcing_ug = ug
                _val_ls_forcing_vg = vg
                nc_wrf.close()

                print("\tProcessing interpolated file: ",interp_files[ts])
                infile = netCDF4.Dataset(interp_files[ts], "r", format="NETCDF4")
                # surface pressure
                surface_forcing_surface_pressure = infile.variables['surface_forcing_surface_pressure']
                _val_surface_forcing_surface_pressure = outfile.variables['surface_forcing_surface_pressure']
                _val_surface_forcing_surface_pressure[ts] = np.average(surface_forcing_surface_pressure[:,:,:], axis = (1,2))[0]
                # other dynamical & chemical variables except: soil_t, soli_m, u, v, w
                for var in dynam_chem_variables:
                    if (var == 'soil_m' or var == 'soil_t' or var == 'u' or var == 'v' or var == 'w'):
                        continue
                    else:
                        if (var == 'h2so4' or var=='hno3' or var== 'nh3' or var== 'ocnv' or var== 'ocsv' or var == 'no' or var == 'no2' or var == 'no3'):
                            var = var.upper()
                        
                        init_var = infile.variables['init_atmosphere_'+ var.lower()]
                        for side in boundary:
                            _val_ls_forcing_var = outfile.variables['ls_forcing_'+ side +'_'+var]
                            if side == 'left':
                                _val_ls_forcing_var[ts, :, :] = init_var[0, :, :, 0]
                            elif side == 'right':
                                _val_ls_forcing_var[ts, :, :] = init_var[0, :, :, dimensions['xdim'] - 1]
                            elif side == 'south':
                                _val_ls_forcing_var[ts, :, :] = init_var[0, :, 0, :]
                            elif side == 'north':
                                _val_ls_forcing_var[ts, :, :] = init_var[0, :, dimensions['ydim'] - 1, :]
                            elif side == 'top':
                                _val_ls_forcing_var[ts, :, :] = init_var[0, dimensions['zdim'] - 1, :, :]
                
                # dynamical variables:u, v, w - mass balancing
                init_atmosphere_u  = infile.variables['init_atmosphere_u']
                init_atmosphere_v  = infile.variables['init_atmosphere_v']
                init_atmosphere_w  = infile.variables['init_atmosphere_w']

                uxleft  = init_atmosphere_u[0, :, :, 0]
                uxright = init_atmosphere_u[0, :, :, dimensions['xdim'] - 1]
                vysouth = init_atmosphere_v[0, :, 0, :]
                vynorth = init_atmosphere_v[0, :, dimensions['ydim'] - 1, :]
                wztop   = init_atmosphere_w[0, dimensions['zwdim'] - 1, :, :]
                mass_disbalance = ((uxleft * areas_xb).sum()
                        - (uxright * areas_xb).sum()
                        + (vysouth * areas_yb).sum()
                        - (vynorth * areas_yb).sum()
                        - (wztop * areas_zb).sum())
                mass_corr_v = mass_disbalance / area_boundaries
                #print('Mass disbalance: {0:8g} m3/s (avg = {1:8g} m/s)'.format(mass_disbalance, mass_corr_v))
                uxleft  -= mass_corr_v
                uxright += mass_corr_v
                vysouth -= mass_corr_v
                vynorth += mass_corr_v
                wztop   += mass_corr_v
                
                _val_ls_forcing_u_left  = outfile.variables['ls_forcing_left_u']
                _val_ls_forcing_u_right = outfile.variables['ls_forcing_right_u']
                _val_ls_forcing_u_south = outfile.variables['ls_forcing_south_u']
                _val_ls_forcing_u_north = outfile.variables['ls_forcing_north_u']
                _val_ls_forcing_u_top   = outfile.variables['ls_forcing_top_u']

                _val_ls_forcing_u_left[ts, :, :] = uxleft
                _val_ls_forcing_u_right[ts, :, :] = uxright
                _val_ls_forcing_u_south[ts, :, :] = init_atmosphere_u[0, :, 0, 1:]
                _val_ls_forcing_u_north[ts, :, :] = init_atmosphere_u[0, :, dimensions['ydim'] - 1, 1:]
                _val_ls_forcing_u_top[ts, :, :] = init_atmosphere_u[0, dimensions['zdim'] - 1, :, 1:]

                _val_ls_forcing_v_left  = outfile.variables['ls_forcing_left_v']
                _val_ls_forcing_v_right = outfile.variables['ls_forcing_right_v']
                _val_ls_forcing_v_south = outfile.variables['ls_forcing_south_v']
                _val_ls_forcing_v_north = outfile.variables['ls_forcing_north_v']
                _val_ls_forcing_v_top   = outfile.variables['ls_forcing_top_v']

                _val_ls_forcing_v_left[ts, :, :] = init_atmosphere_v[0, :, 1:, 0]
                _val_ls_forcing_v_right[ts, :, :] = init_atmosphere_v[0, :, 1:, dimensions['xdim'] - 1]
                _val_ls_forcing_v_south[ts, :, :] = vysouth
                _val_ls_forcing_v_north[ts, :, :] = vynorth
                _val_ls_forcing_v_top[ts, :, :] = init_atmosphere_v[0, dimensions['zdim'] - 1, 1:, :]

                _val_ls_forcing_w_left  = outfile.variables['ls_forcing_left_w']
                _val_ls_forcing_w_right = outfile.variables['ls_forcing_right_w']
                _val_ls_forcing_w_south = outfile.variables['ls_forcing_south_w']
                _val_ls_forcing_w_north = outfile.variables['ls_forcing_north_w']
                _val_ls_forcing_w_top   = outfile.variables['ls_forcing_top_w']

                _val_ls_forcing_w_left[ts, :, :] = init_atmosphere_w[0, :, :, 0]
                _val_ls_forcing_w_right[ts, :, :] = init_atmosphere_w[0, :, :, dimensions['xdim'] - 1]
                _val_ls_forcing_w_south[ts, :, :] = init_atmosphere_w[0, :, 0, :]
                _val_ls_forcing_w_north[ts, :, :] = init_atmosphere_w[0, :, dimensions['ydim'] - 1, :]
                _val_ls_forcing_w_top[ts, :, :] = wztop

                infile.close()
            outfile.close()
    # read interpolated files and write values to outfile
    add_interpValues(dynam_chem_variables)

    # add radiation input
    if len(rad_times_proc) > 0:
        outfile = netCDF4.Dataset(dynamic_driver_file, "a", format="NETCDF4")
        # radiation time dimension and variable
        outfile.createDimension('time_rad', len(rad_times_proc))
        _val_times =  outfile.createVariable('time_rad',"f4", ("time_rad"))
        _val_times[:] = rad_times_proc[:]
        # radiation variables
        var = outfile.createVariable('rad_sw_in', "f4", ("time_rad"), fill_value=fillvalue_float)
        var.setncattr('lod', 1)
        var.units = 'W/m2'
        var[:] = rad_values_proc[0][:]
        var = outfile.createVariable('rad_lw_in', "f4", ("time_rad"), fill_value=fillvalue_float)
        var.setncattr('lod', 1)
        var.units = 'W/m2'
        var[:] = rad_values_proc[1][:]
        var = outfile.createVariable('rad_sw_in_dif', "f4", ("time_rad"), fill_value=fillvalue_float)
        var.setncattr('lod', 1)
        var.units = 'W/m2'
        var[:] = rad_values_proc[2][:]

        outfile.close()
    # correct chemical species names to PALM output
    outfile = netCDF4.Dataset(dynamic_driver_file, "a", format="NETCDF4")
    if 'PM2_5_DRY' in dynam_chem_variables:
        outfile.renameVariable('init_atmosphere_PM2_5_DRY','init_atmosphere_pm2_5')
        outfile.renameVariable('ls_forcing_left_PM2_5_DRY', 'ls_forcing_left_pm2_5')
        outfile.renameVariable('ls_forcing_right_PM2_5_DRY','ls_forcing_right_pm2_5')
        outfile.renameVariable('ls_forcing_south_PM2_5_DRY','ls_forcing_south_pm2_5')
        outfile.renameVariable('ls_forcing_north_PM2_5_DRY','ls_forcing_north_pm2_5')
        outfile.renameVariable('ls_forcing_top_PM2_5_DRY','ls_forcing_top_pm2_5')
    if 'PM10' in dynam_chem_variables:
        outfile.renameVariable('init_atmosphere_PM10','init_atmosphere_pm10')
        outfile.renameVariable('ls_forcing_left_PM10', 'ls_forcing_left_pm10')
        outfile.renameVariable('ls_forcing_right_PM10','ls_forcing_right_pm10')
        outfile.renameVariable('ls_forcing_south_PM10','ls_forcing_south_pm10')
        outfile.renameVariable('ls_forcing_north_PM10','ls_forcing_north_pm10')
        outfile.renameVariable('ls_forcing_top_PM10','ls_forcing_top_pm10')
    outfile.close()

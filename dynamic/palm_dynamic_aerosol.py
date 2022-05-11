#------------------------------------------------------------------------------
#
# Script for processing of WRF-CHEM files to PALM dynamic driver.
#
#------------------------------------------------------------------------------
import os
import netCDF4
import numpy as np
from palm_dynamic_config import *

def translate_aerosol_species(name):
    testname = name
    translation_table = {
            "SO4": "so4",
            "NO" : "no3",
            "NH" : "nh4",
            "BC" : "bc",
            "OC" : "smpa,smpbb,glysoa_sfc,biog1_c,biog1_o",
            "SS" : "cl,na",
            "DU" : "co3,ca,oin"
            }
    if testname in translation_table:
        return translation_table[testname]
    else:
        return None

def aerosolProfileMass(interp_files, listspec):

    infile  = netCDF4.Dataset(interp_files[0], "r", format="NETCDF4")
    alt     = infile.variables['init_atmosphere_alt'][0]
    u       = infile.variables['init_atmosphere_u']
    v       = infile.variables['init_atmosphere_v']
    z_level = infile.variables['z']

    # upwind location & mass frac
    aero_massfrac_a = np.zeros((z_level.size, len(listspec)))

    for zlev in range(0, z_level.size):
        u_wnd = u[0,zlev,:,:]
        v_wnd = v[0,zlev,:,:]
        wnd_dir = np.mod(180 + np.rad2deg(np.arctan2(u_wnd, v_wnd)),360)
        wnd_avg = np.mean(wnd_dir)
        if  0 < wnd_avg <= 45:
            prf_x = round(wnd_dir.shape[0]/2)
            prf_y = 0
        elif 315 < wnd_avg <=0:
            prf_x = round(wnd_dir.shape[0]/2)
            prf_y = 0
        elif 45 < wnd_avg <= 135:
            prf_x = round(wnd_dir.shape[0]/2)
            prf_y = round(wnd_dir.shape[1]/2)
        elif 135 < wnd_avg <= 225:
            prf_x = round(wnd_dir.shape[0]/2)
            prf_y = wnd_dir.shape[1]
        else:
            prf_x = round(wnd_dir.shape[0]/2)
            prf_y = round(wnd_dir.shape[1]/2)
        # mass at each z-level (z, composition_index)
        aero_mass = np.zeros((len(listspec)))

        for naero in range(0, len(listspec)):
            in_spec = listspec[naero]
            if (in_spec == 'H2SO4' or in_spec=='HNO3' or in_spec== 'NH3' or in_spec== 'OCNV' or in_spec== 'OCSV'):
                spec_mass = infile.variables['init_atmosphere_'+ in_spec.lower()][0]
            else:
                spec_mass = infile.variables['init_atmosphere_'+ in_spec][0]
            aero_mass[naero] = spec_mass[zlev, prf_y, prf_x]
        total_mass = np.sum(aero_mass)

        # mass fraction values
        for naero in range(0, len(listspec)):
            aero_massfrac_a[zlev, naero] = aero_mass[naero]/total_mass

    return aero_massfrac_a

    infile.close()

def aerosolProfileConc(interp_files, nbin, reglim, wrfchem_bin_limits):

    infile  = netCDF4.Dataset(interp_files[0], "r", format="NETCDF4")
    alt     = infile.variables['init_atmosphere_alt'][0]
    u       = infile.variables['init_atmosphere_u']
    v       = infile.variables['init_atmosphere_v']
    z_level = infile.variables['z']

    # define palm bins
    bin_dmid, bin_lims = define_bins(nbin, reglim)
    # define overlap of palm & wrfchem bins
    open_bins, overlap_ratio = aerosol_binoverlap(bin_lims, wrfchem_bin_limits)
    open_bins = sorted(set(open_bins), key=open_bins.index)

    aero_con = np.zeros((z_level.size, bin_dmid.size))

    for zlev in range(0, z_level.size):
        u_wnd = u[0,zlev,:,:]
        v_wnd = v[0,zlev,:,:]
        wnd_dir = np.mod(180 + np.rad2deg(np.arctan2(u_wnd, v_wnd)),360)
        wnd_avg = np.mean(wnd_dir)
        if  0 < wnd_avg <= 45:
            prf_x = round(wnd_dir.shape[0]/2)
            prf_y = 0
        elif 315 < wnd_avg <=0:
            prf_x = round(wnd_dir.shape[0]/2)
            prf_y = 0
        elif 45 < wnd_avg <= 135:
            prf_x = round(wnd_dir.shape[0]/2)
            prf_y = round(wnd_dir.shape[1]/2)
        elif 135 < wnd_avg <= 225:
            prf_x = round(wnd_dir.shape[0]/2)
            prf_y = wnd_dir.shape[1]
        else:
            prf_x = round(wnd_dir.shape[0]/2)
            prf_y = round(wnd_dir.shape[1]/2)

        # aerosol concen#
        # ug/kg-dryair to aerosol# concen.(# m-3):inverse density (m3 kg-1)
        # convert units, weight Dmid by size-bin use overlap ratio
        for n_dmid in range(0, bin_dmid.size):
            outval = 0.0
            for abin in range(0, len(open_bins)):
                inval = infile.variables['init_aerosol'+ open_bins[abin]][0]
                _val = inval[zlev, prf_y, prf_x] * alt[zlev, prf_y, prf_x]
                _val1 = _val * overlap_ratio[n_dmid, abin]
                outval = outval + _val1
            aero_con[zlev, n_dmid] = outval

    return aero_con

    infile.close()

def aerosolMassWrfchemBoundary(dimensions, _dim1, _dim2, interp_files, listspec, side):
    _dim1_size = dimensions[_dim1+'dim']
    _dim2_size = dimensions[_dim2+'dim']

    # mass_frac_a: val_side[time, z, y, composition_index]
    val_side = np.zeros( (len(interp_files), _dim1_size, _dim2_size, len(listspec)) )

    for ts in range(0, len(interp_files)):
        infile = netCDF4.Dataset(interp_files[ts], "r", format="NETCDF4")
        val_spec = np.zeros( (_dim1_size, _dim2_size, len(listspec)) )

        for n_spec in range(0, len(listspec)):
            in_spec = listspec[n_spec]
            if (in_spec == 'H2SO4' or in_spec=='HNO3' or in_spec== 'NH3' or in_spec== 'OCNV' or in_spec== 'OCSV'):
                in_spec = in_spec.lower()

            if (side == 'left'):
                val_spec[:,:,n_spec] = infile.variables['init_atmosphere_'+ in_spec][0, :, :, 0]
            elif (side == 'right'):
                val_spec[:,:,n_spec] = infile.variables['init_atmosphere_'+ in_spec][0, :, :, dimensions['xdim'] - 1]
            elif (side =='south'):
                val_spec[:,:,n_spec] = infile.variables['init_atmosphere_'+ in_spec][0, :, 0, :]
            elif (side == 'north'):
                val_spec[:,:,n_spec] = infile.variables['init_atmosphere_'+ in_spec][0, :, dimensions['ydim'] - 1, :]
            elif (side == 'top'):
                val_spec[:,:,n_spec] = infile.variables['init_atmosphere_'+ in_spec][0, dimensions['zdim']-1, :, :]
        val_sum = np.sum(val_spec,axis = 2)

        # write mass frac values
        for n_spec in range(0, len(listspec)):
            val_side[ts, :, :, n_spec] = val_spec[:,:,n_spec]/val_sum
        infile.close()

    return val_side

def aerosolConWrfchemBoundary(dimensions, _dim1, _dim2, interp_files, side, nbin, reglim, wrfchem_bin_limits):
    # define palm bins
    bin_dmid, bin_lims = define_bins(nbin, reglim)
    # define overlap of palm & wrfchem bins
    open_bins, overlap_ratio = aerosol_binoverlap(bin_lims, wrfchem_bin_limits)
    open_bins = sorted(set(open_bins), key=open_bins.index)
    # adjust dimensions
    _dim1_size = dimensions[_dim1+'dim']
    _dim2_size = dimensions[_dim2+'dim']
    # aerosol concentration number [time, z, y, Dmid]
    aero_conc = np.zeros((len(interp_files), _dim1_size, _dim2_size, bin_dmid.size))
    for ts in range(0, len(interp_files)):
        infile = netCDF4.Dataset(interp_files[ts], "r", format="NETCDF4")
        outval = np.zeros( (_dim1_size, _dim2_size) )
        for n_dmid in range(0, bin_dmid.size):

            for abin in range(0, len(open_bins)):
                if (side == 'left'):
                    alt   = infile.variables['init_atmosphere_alt'][0, :, :, 0]
                    inval = infile.variables['init_aerosol'+ open_bins[abin]][0, :, :, 0]
                elif (side == 'right'):
                    alt   = infile.variables['init_atmosphere_alt'][0, :, :, dimensions['xdim'] - 1]
                    inval = infile.variables['init_aerosol'+ open_bins[abin]][0, :, :, dimensions['xdim'] - 1]
                elif (side =='south'):
                    alt   = infile.variables['init_atmosphere_alt'][0, :, 0, :]
                    inval = infile.variables['init_aerosol'+ open_bins[abin]][0, :, 0, :]
                elif (side == 'north'):
                    alt   = infile.variables['init_atmosphere_alt'][0, :, dimensions['ydim'] - 1, :]
                    inval = infile.variables['init_aerosol'+ open_bins[abin]][0, :, dimensions['ydim'] - 1, :]
                elif (side == 'top'):
                    alt   = infile.variables['init_atmosphere_alt'][0, dimensions['zdim']-1, :, :]
                    inval = infile.variables['init_aerosol'+ open_bins[abin]][0, dimensions['zdim']-1, :, :]
                # convert and factor for size-bin
                inval = inval * alt
                inval = inval * overlap_ratio[n_dmid, abin]
                outval = outval + inval

            aero_conc[ts, :, :, n_dmid] = outval
        infile.close()

    return aero_conc

def define_bins(nbin, reglim):
    # degine Dmid based on nbin and reglim
    nbins = np.sum(nbin) # = subrange 1 + subrange 2
    # Log-normal to sectional
    vlolim = np.zeros(nbins)
    vhilim = np.zeros(nbins)
    dmid   = np.zeros(nbins)
    bin_limits = np.zeros(nbins)
    # Sectional bin limits
    ratio_d = reglim[1] / reglim[0]
    for b in range(nbin[0]):
      vlolim[b] = np.pi / 6.0 * (reglim[0] * ratio_d **(float(b) / nbin[0]))**3
      vhilim[b] = np.pi / 6.0 * (reglim[0] * ratio_d **(float(b+1) / nbin[0]))**3
      dmid[b] = np.sqrt((6.0 * vhilim[b] / np.pi )**0.33333333 * (6.0 * vlolim[b] / np.pi)**0.33333333)
    ratio_d = reglim[2] / reglim[1]

    for b in np.arange(nbin[0], np.sum(nbin),1):
      c = b-nbin[0]
      vlolim[b] = np.pi / 6.0 * (reglim[1] * ratio_d **(float(c) / nbin[1]))**3
      vhilim[b] = np.pi / 6.0 * (reglim[1] * ratio_d **(float(c+1) / nbin[1])) ** 3
      dmid[b] = np.sqrt((6.0 * vhilim[b] / np.pi)**0.33333333 * ( 6.0 * vlolim[b] / np.pi)**0.33333333)
    bin_limits = (6.0 * vlolim / np.pi )**0.33333333
    bin_limits = np.append(bin_limits, reglim[-1])

    return dmid, bin_limits

def range_overlap(range1, range2):
  x1, x2 = range1.start, range1.stop
  y1, y2 = range2.start, range2.stop

  return x1 <= y2 and y1 <= x2

def aerosol_binoverlap(palm_binlim, wrfchem_binlim):
    overlap_ratio = np.zeros( (len(palm_binlim)-1, len(wrfchem_binlim)-1) )
    aerobin_open = []
    for pbin in range(0,len(palm_binlim)-1):
        palm_range = range(int(palm_binlim[pbin]* 1e+9), int(palm_binlim[pbin+1]* 1e+9))
        for wbin in range(0, len(wrfchem_binlim)-1):
            wrfchem_range = range(int(wrfchem_binlim[wbin]* 1e+9)+1, int(wrfchem_binlim[wbin+1]* 1e+9))

            if range_overlap(palm_range, wrfchem_range):
                aerobin_open.append('_a0'+ str(wbin+1))
                overlap = len(set(palm_range) & set(wrfchem_range))
                overlap_ratio[pbin, wbin] = overlap/len(wrfchem_range)

    return aerobin_open, overlap_ratio

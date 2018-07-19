# * Copyright (C) 2017 J.Mittaz University of Reading
# * This code was developed for the EC project Fidelity and Uncertainty in
# * Climate Data Records from Earth Observations (FIDUCEO). 
# * Grant Agreement: 638822
# *
# * This program is free software; you can redistribute it and/or modify it
# * under the terms of the GNU General Public License as published by the Free
# * Software Foundation; either version 3 of the License, or (at your option)
# * any later version.
# * This program is distributed in the hope that it will be useful, but WITHOUT
# * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# * more details.
# * 
# * A copy of the GNU General Public License should have been supplied along
# * with this program; if not, see http://www.gnu.org/licenses/
# * ------------------------------------------------------------------------
# * MT: 18-10-2017: added quality flag fields
# * MT: 30-10-2017: uncertainty variables renamed in line with FCDR-CDR file format spec fv1.1.1
# * MT: 09-11-2017: channel_correlation_matrix added (sensor specific)
# * MT: 10-11-2017: spatial_correlation_scale added (sensor specific)
# * JM: 06-07-2018: Real channel/spatial correlation scales added plus SRF
# * JM: 07-07-2018: Output GBCS L1C option for SST with channel covariance

from fiduceo.fcdr.writer.fcdr_writer import FCDRWriter
from fiduceo.fcdr.writer.templates import avhrr
import sys
import netCDF4
import numpy as np
import datetime
import xarray
from optparse import OptionParser
import FCDR_HIRS.metrology as met
import write_l1c_data as l1c

class read_netcdf(object):

    def add_nan_values(self,values):
        with np.errstate(invalid='ignore'):
            gd = np.isfinite(values) & (values < -1e20)
        if np.sum(gd) > 0:
            values[gd] = np.NaN
        return values

    def scale_values(self,values):
        with np.errstate(invalid='ignore'):
            gd = (values > -1e20) & np.isfinite(values)
        if np.sum(gd) > 0:
            values[gd] = values[gd]*100.
        return values

    def read_data(self,filename):

        ncid = netCDF4.Dataset(filename,'r')

        self.sources = ncid.sources
        self.noaa_string = ncid.noaa_string
        self.version = ncid.version
        self.spatial_correlation_scale = ncid.spatial_correlation_scale
        year = ncid.variables['year'][:]
        month = ncid.variables['month'][:]
        day = ncid.variables['day'][:]
        hours = ncid.variables['hours'][:]
        gd = (year > 1900) & (day > 0) & (month > 0) & (hours >= 0)
        if 0 == np.sum(gd):
            raise Exception('No good data in netcdf')
        yr = year[gd]
        mn = month[gd]
        dy = day[gd]
        temp = hours[gd]
        hour = temp.astype(np.int)
        temp = (temp - hour)*60.
        minute = temp.astype(np.int)
        temp = (temp - minute)*60.
        second = temp.astype(np.int)
        microsec = ((temp-second)*1e6).astype(np.int)
        self.date_time = []
        for i in range(len(hour)):
            self.date_time.append(datetime.datetime(yr[i],mn[i],dy[i],hour[i],minute[i],second[i],microsec[i]))
        self.time = netCDF4.date2num(self.date_time,'seconds since 1970-01-01')
        lat = ncid.variables['latitude'][:,:]
        self.lat = lat[gd,:]
        lon = ncid.variables['longitude'][:,:]
        self.lon = lon[gd,:]
        satza = ncid.variables['satza'][:,:]
        self.satza = satza[gd,:]
        solza = ncid.variables['solza'][:,:]
        self.solza = solza[gd,:]
        relaz = ncid.variables['relaz'][:,:]
        self.relaz = relaz[gd,:]
        ch1 = ncid.variables['ch1'][:,:]
        self.ch1 = ch1[gd,:]
        ch2 = ncid.variables['ch2'][:,:]
        self.ch2 = ch2[gd,:]
        ch3a = ncid.variables['ch3a'][:,:]
        self.ch3a = ch3a[gd,:]
        self.ch3a_there=False
        if np.ma.is_masked(self.ch3a):
            if np.any(~self.ch3a.mask):
                self.ch3a_there=True
        if np.any(np.isfinite(self.ch3a)):
            self.ch3a_there=True            
        ch3b = ncid.variables['ch3b'][:,:]
        self.ch3b = ch3b[gd,:]
        ch4 = ncid.variables['ch4'][:,:]
        self.ch4 = ch4[gd,:]
        try:
            ch5 = ncid.variables['ch5'][:,:]
            self.ch5 = ch5[gd,:]
            self.ch5_there = True
        except:
            self.ch5_there = False
        u_random_ch1 = ncid.variables['ch1_random'][:,:]
        self.u_random_ch1 = u_random_ch1[gd,:]
        u_random_ch2 = ncid.variables['ch2_random'][:,:]
        self.u_random_ch2 = u_random_ch2[gd,:]
        if self.ch3a_there:
            u_random_ch3a = ncid.variables['ch3a_random'][:,:]
            self.u_random_ch3a = u_random_ch3a[gd,:]
        u_random_ch3b = ncid.variables['ch3b_random'][:,:]
        self.u_random_ch3b = u_random_ch3b[gd,:]
        u_random_ch4 = ncid.variables['ch4_random'][:,:]
        self.u_random_ch4 = u_random_ch4[gd,:]
        if self.ch5_there:
            u_random_ch5 = ncid.variables['ch5_random'][:,:]
            self.u_random_ch5 = u_random_ch5[gd,:]
        u_non_random_ch1 = ncid.variables['ch1_non_random'][:,:]
        self.u_non_random_ch1 = u_non_random_ch1[gd,:]
        u_non_random_ch2 = ncid.variables['ch2_non_random'][:,:]
        self.u_non_random_ch2 = u_non_random_ch2[gd,:]
        if self.ch3a_there:
            u_non_random_ch3a = ncid.variables['ch3a_non_random'][:,:]
            self.u_non_random_ch3a = u_non_random_ch3a[gd,:]
        u_non_random_ch3b = ncid.variables['ch3b_non_random'][:,:]
        self.u_non_random_ch3b = u_non_random_ch3b[gd,:]
        u_non_random_ch4 = ncid.variables['ch4_non_random'][:,:]
        self.u_non_random_ch4 = u_non_random_ch4[gd,:]
        if self.ch5_there:
            u_non_random_ch5 = ncid.variables['ch5_non_random'][:,:]
            self.u_non_random_ch5 = u_non_random_ch5[gd,:]
        u_common_ch3b = ncid.variables['ch3b_common'][:,:]
        self.u_common_ch3b = u_common_ch3b[gd,:]
        u_common_ch4 = ncid.variables['ch4_common'][:,:]
        self.u_common_ch4 = u_common_ch4[gd,:]
        if self.ch5_there:
            u_common_ch5 = ncid.variables['ch5_common'][:,:]
            self.u_common_ch5 = u_common_ch5[gd,:]
        u_common_ch1 = ncid.variables['ch1_common'][:,:]
        self.u_common_ch1 = u_common_ch1[gd,:]
        u_common_ch2 = ncid.variables['ch2_common'][:,:]
        self.u_common_ch2 = u_common_ch2[gd,:]
        if self.ch3a_there:
            u_common_ch3a = ncid.variables['ch3a_common'][:,:]
            self.u_common_ch5 = u_common_ch3a[gd,:]
        scan_qual = ncid.variables['quality_scanline_bitmask'][:]
        chan_qual = ncid.variables['quality_channel_bitmask'][:,:]
        self.scan_qual = scan_qual[gd]
        self.chan_qual = chan_qual[gd,:]
        dBT3_over_dT = ncid.variables['dBT3_over_dT'][:,:]
        self.dBT3_over_dT = dBT3_over_dT[gd,:]
        dBT4_over_dT = ncid.variables['dBT4_over_dT'][:,:]
        self.dBT4_over_dT = dBT4_over_dT[gd,:]
        if self.ch5_there:
            dBT5_over_dT = ncid.variables['dBT5_over_dT'][:,:]
            self.dBT5_over_dT = dBT5_over_dT[gd,:]
        dRe1_over_dCS = ncid.variables['dRe1_over_dCS'][:,:]
        self.dRe1_over_dCS = dRe1_over_dCS[gd,:]
        dRe2_over_dCS = ncid.variables['dRe2_over_dCS'][:,:]
        self.dRe2_over_dCS = dRe2_over_dCS[gd,:]
        if self.ch3a_there:
            dRe3a_over_dCS = ncid.variables['dRe3a_over_dCS'][:,:]
            self.dRe3a_over_dCS = dRe3a_over_dCS[gd,:]
        dBT3_over_dCS = ncid.variables['dBT3_over_dCS'][:,:]
        self.dBT3_over_dCS = dBT3_over_dCS[gd,:]
        dBT4_over_dCS = ncid.variables['dBT4_over_dCS'][:,:]
        self.dBT4_over_dCS = dBT4_over_dCS[gd,:]
        if self.ch5_there:
            dBT5_over_dCS = ncid.variables['dBT5_over_dCS'][:,:]
            self.dBT5_over_dCS = dBT5_over_dCS[gd,:]
        dBT3_over_dCICT = ncid.variables['dBT3_over_dCICT'][:,:]
        self.dBT3_over_dCICT = dBT3_over_dCICT[gd,:]
        dBT4_over_dCICT = ncid.variables['dBT4_over_dCICT'][:,:]
        self.dBT4_over_dCICT = dBT4_over_dCICT[gd,:]
        if self.ch5_there:
            dBT5_over_dCICT = ncid.variables['dBT5_over_dCICT'][:,:]
            self.dBT5_over_dCICT = dBT5_over_dCICT[gd,:]
        smoothPRT = ncid.variables['dBT5_over_dCICT'][:]
        self.smoothPRT = smoothPRT[gd]
        self.cal_cnts_noise = ncid.variables['cal_cnts_noise'][:]
        self.cnts_noise = ncid.variables['cnts_noise'][:]
        self.spatial_correlation_scale = ncid.spatial_correlation_scale
        self.ICT_Temperature_Uncertainty = ncid.ICT_Temperature_Uncertainty
        self.noaa_string = ncid.noaa_string
        self.orbital_temperature = ncid.orbital_temperature
        scanline = ncid.variables['scanline'][:]
        self.scanline = scanline[gd]
        orig_scanline = ncid.variables['orig_scanline'][:]
        self.orig_scanline = orig_scanline[gd] 

        badNav = ncid.variables['badNavigation'][:]
        self.badNav = badNav[gd]
        badCal = ncid.variables['badCalibration'][:]
        self.badCal = badCal[gd]
        badTime = ncid.variables['badTime'][:]
        self.badTime = badTime[gd]
        missingLines = ncid.variables['missingLines'][:]
        self.missingLines = missingLines[gd]
        solar3 = ncid.variables['solar_contam_3b'][:]
        self.solar3 = solar3[gd]
        solar4 = ncid.variables['solar_contam_4'][:]
        self.solar4 = solar4[gd]
        solar5 = ncid.variables['solar_contam_5'][:]
        self.solar5 = solar5[gd]

        self.nx = self.lat.shape[1]
        self.ny = self.lat.shape[0]

        ncid.close()

        self.lat = np.ma.filled(self.lat,np.NaN)
        self.lon = np.ma.filled(self.lon,np.NaN)
        self.time = np.ma.filled(self.time,np.NaN)
        self.satza = np.ma.filled(self.satza,np.NaN)
        self.solza = np.ma.filled(self.solza,np.NaN)
        self.relaz = np.ma.filled(self.relaz,np.NaN)
        self.ch1 = np.ma.filled(self.ch1,np.NaN)
        self.ch2 = np.ma.filled(self.ch2,np.NaN)
        if self.ch3a_there:
            self.ch3a = np.ma.filled(self.ch3a,np.NaN)
        self.ch3b = np.ma.filled(self.ch3b,np.NaN)
        self.ch4 = np.ma.filled(self.ch4,np.NaN)
        if self.ch5_there:
            self.ch5 = np.ma.filled(self.ch5,np.NaN)
        self.u_random_ch1 = np.ma.filled(self.u_random_ch1,np.NaN)
        self.u_random_ch2 = np.ma.filled(self.u_random_ch2,np.NaN)
        if self.ch3a_there:
            self.u_random_ch3a = np.ma.filled(self.u_random_ch3a,np.NaN)
        self.u_random_ch3b = np.ma.filled(self.u_random_ch3b,np.NaN)
        self.u_random_ch4 = np.ma.filled(self.u_random_ch4,np.NaN)
        if self.ch5_there:
            self.u_random_ch5 = np.ma.filled(self.u_random_ch5,np.NaN)
        self.u_non_random_ch1 = np.ma.filled(self.u_non_random_ch1,np.NaN)
        self.u_non_random_ch2 = np.ma.filled(self.u_non_random_ch2,np.NaN)
        if self.ch3a_there:
            self.u_non_random_ch3a = np.ma.filled(self.u_non_random_ch3a,np.NaN)
        self.u_non_random_ch3b = np.ma.filled(self.u_non_random_ch3b,np.NaN)
        self.u_non_random_ch4 = np.ma.filled(self.u_non_random_ch4,np.NaN)
        if self.ch5_there:
            self.u_non_random_ch5 = np.ma.filled(self.u_non_random_ch5,np.NaN)
        self.u_common_ch1 = np.ma.filled(self.u_common_ch1,np.NaN)
        self.u_common_ch2 = np.ma.filled(self.u_common_ch2,np.NaN)
        if self.ch3a_there:
            self.u_common_ch3a = np.ma.filled(self.u_common_ch3a,np.NaN)
        self.u_common_ch3b = np.ma.filled(self.u_common_ch3b,np.NaN)
        self.u_common_ch4 = np.ma.filled(self.u_common_ch4,np.NaN)
        if self.ch5_there:
            self.u_common_ch5 = np.ma.filled(self.u_common_ch5,np.NaN)
        self.dBT3_over_dT = np.ma.filled(self.dBT3_over_dT,np.NaN)
        self.dBT4_over_dT = np.ma.filled(self.dBT4_over_dT,np.NaN)
        if self.ch5_there:
            self.dBT5_over_dT = np.ma.filled(self.dBT5_over_dT,np.NaN)
        self.dRe1_over_dCS = np.ma.filled(self.dRe1_over_dCS,np.NaN)
        self.dRe2_over_dCS = np.ma.filled(self.dRe2_over_dCS,np.NaN)
        if self.ch3a_there:
            self.dRe3a_over_dCS = np.ma.filled(self.dRe3a_over_dCS,np.NaN)
        self.dBT3_over_dCS = np.ma.filled(self.dBT3_over_dCS,np.NaN)
        self.dBT4_over_dCS = np.ma.filled(self.dBT4_over_dCS,np.NaN)
        if self.ch5_there:
            self.dBT5_over_dCS = np.ma.filled(self.dBT5_over_dCS,np.NaN)
        self.dBT3_over_dCICT = np.ma.filled(self.dBT3_over_dCICT,np.NaN)
        self.dBT4_over_dCICT = np.ma.filled(self.dBT4_over_dCICT,np.NaN)
        if self.ch5_there:
            self.dBT5_over_dCICT = np.ma.filled(self.dBT5_over_dCICT,np.NaN)
        self.cal_cnts_noise = np.ma.filled(self.cal_cnts_noise,np.NaN)
        self.cnts_noise = np.ma.filled(self.cnts_noise,np.NaN)

        self.lat = self.add_nan_values(self.lat)
        self.lon = self.add_nan_values(self.lon)
        self.time = self.add_nan_values(self.time)
        self.satza = self.add_nan_values(self.satza)
        self.relaz = self.add_nan_values(self.relaz)
        self.ch1 = self.add_nan_values(self.ch1)
        self.ch2 = self.add_nan_values(self.ch2)
        if self.ch3a_there:
            self.ch3a = self.add_nan_values(self.ch3a)
        self.ch3b = self.add_nan_values(self.ch3b)
        self.ch4 = self.add_nan_values(self.ch4)
        if self.ch5_there:
            self.ch5 = self.add_nan_values(self.ch5)
        self.u_random_ch1 = self.add_nan_values(self.u_random_ch1)
        self.u_random_ch2 = self.add_nan_values(self.u_random_ch2)
        if self.ch3a_there:
            self.u_random_ch3a = self.add_nan_values(self.u_random_ch3a)
        self.u_random_ch3b = self.add_nan_values(self.u_random_ch3b)
        self.u_random_ch4 = self.add_nan_values(self.u_random_ch4)
        if self.ch5_there:
            self.u_random_ch5 = self.add_nan_values(self.u_random_ch5)
        self.u_non_random_ch1 = self.add_nan_values(self.u_non_random_ch1)
        self.u_non_random_ch2 = self.add_nan_values(self.u_non_random_ch2)
        if self.ch3a_there:
            self.u_non_random_ch3a = self.add_nan_values(self.u_non_random_ch3a)
        self.u_non_random_ch3b = self.add_nan_values(self.u_non_random_ch3b)
        self.u_non_random_ch4 = self.add_nan_values(self.u_non_random_ch4)
        if self.ch5_there:
            self.u_non_random_ch5 = self.add_nan_values(self.u_non_random_ch5)
        self.u_common_ch1 = self.add_nan_values(self.u_common_ch1)
        self.u_common_ch2 = self.add_nan_values(self.u_common_ch2)
        if self.ch3a_there:
            self.u_common_ch3a = self.add_nan_values(self.u_common_ch3a)
        self.u_common_ch3b = self.add_nan_values(self.u_common_ch3b)
        self.u_common_ch4 = self.add_nan_values(self.u_common_ch4)
        if self.ch5_there:
            self.u_common_ch5 = self.add_nan_values(self.u_common_ch5)
        self.dBT3_over_dT = self.add_nan_values(self.dBT3_over_dT)
        self.dBT4_over_dT = self.add_nan_values(self.dBT4_over_dT)
        if self.ch5_there:
            self.dBT5_over_dT = self.add_nan_values(self.dBT5_over_dT)
        self.dRe1_over_dCS = self.add_nan_values(self.dRe1_over_dCS)
        self.dRe2_over_dCS = self.add_nan_values(self.dRe2_over_dCS)
        if self.ch3a_there:
            self.dRe3a_over_dCS = self.add_nan_values(self.dRe3a_over_dCS)
        self.dBT3_over_dCS = self.add_nan_values(self.dBT3_over_dCS)
        self.dBT4_over_dCS = self.add_nan_values(self.dBT4_over_dCS)
        if self.ch5_there:
            self.dBT5_over_dCS = self.add_nan_values(self.dBT5_over_dCS)
        self.dBT3_over_dCICT = self.add_nan_values(self.dBT3_over_dCICT)
        self.dBT5_over_dCICT = self.add_nan_values(self.dBT4_over_dCICT)
        if self.ch5_there:
            self.dBT5_over_dCICT = self.add_nan_values(self.dBT5_over_dCICT)
        self.cal_cnts_noise = self.add_nan_values(self.cal_cnts_noise)
        self.cnts_noise = self.add_nan_values(self.cnts_noise)
            
#        self.ch1 = self.scale_values(self.ch1)
#        self.ch2 = self.scale_values(self.ch2)
#        if self.ch3a_there:
#            self.ch3a = self.scale_values(self.ch3a)
#        self.u_random_ch1 = self.scale_values(self.u_random_ch1)
#        self.u_random_ch2 = self.scale_values(self.u_random_ch2)
#        if self.ch3a_there:
#            self.u_random_ch3a = self.scale_values(self.u_random_ch3a)
#        self.u_non_random_ch1 = self.scale_values(self.u_non_random_ch1)
#        self.u_non_random_ch2 = self.scale_values(self.u_non_random_ch2)
#        if self.ch3a_there:
#            self.u_non_random_ch3a = self.scale_values(self.u_non_random_ch3a)

    def __init__(self,filename):

        self.read_data(filename)
#
# Run Gerrits CURUC routines
#
#
# Force bad data to be nan's
#
def set_to_nan(TL):

    # Set -1e30 to NaN
    gd = np.isfinite(TL)
    newTL = TL[gd]
    gd2 = (newTL < -1e20)
    if np.sum(gd2) > 0:
        newTL[gd2] = float('nan')
        TL[gd] = newTL
    if np.ma.is_masked(TL):
        gd = ~TL.mask
        TL[gd] = float('nan')

    return TL

#
# Copy over to sensitivity using the coords (xarray) arrays
#
def copy_over_C(inarray,n_l_coord,n_e_coord,inverse=False):

    if inverse:
        outarray = np.zeros((len(n_e_coord),len(n_l_coord)))
        for i in range(len(n_l_coord)):
            outarray[:,i] = inarray[n_l_coord[i],n_e_coord]
    else:
        outarray = np.zeros((len(n_l_coord),len(n_e_coord)))
        for i in range(len(n_l_coord)):
            outarray[i,:] = inarray[n_l_coord[i],n_e_coord]

    return outarray

#
# Do copy in a block and deal with different channel selection
#
def copy_C3(dBT3_over_dX,dBT4_over_dX,dBT5_over_dX,\
                n_l,n_e,chans,inverse=False,xchan=False):

    if xchan:
        if inverse:
            output = np.zeros((len(n_e),len(n_l),len(chans)))
            for i in range(len(chans)):
                if 3 == chans[i]:
                    output[:,:,i] = copy_over_C(dBT3_over_dX,n_l,n_e,inverse=True)
                elif 4 == chans[i]:
                    output[:,:,i] = copy_over_C(dBT4_over_dX,n_l,n_e,inverse=True)
                elif 5 == chans[i]:
                    output[:,:,i] = copy_over_C(dBT5_over_dX,n_l,n_e,inverse=True)
        else:
            output = np.zeros((len(n_l),len(n_e),len(chans)))
            for i in range(len(chans)):
                if 3 == chans[i]:
                    output[:,:,i] = copy_over_C(dBT3_over_dX,n_l,n_e)
                elif 4 == chans[i]:
                    output[:,:,i] = copy_over_C(dBT4_over_dX,n_l,n_e)
                elif 5 == chans[i]:
                    output[:,:,i] = copy_over_C(dBT5_over_dX,n_l,n_e)
    else:
        if inverse:
            output = np.zeros((len(chans),len(n_e),len(n_l)))
            for i in range(len(chans)):
                if 3 == chans[i]:
                    output[i,:,:] = copy_over_C(dBT3_over_dX,n_l,n_e,inverse=True)
                elif 4 == chans[i]:
                    output[i,:,:] = copy_over_C(dBT4_over_dX,n_l,n_e,inverse=True)
                elif 5 == chans[i]:
                    output[i,:,:] = copy_over_C(dBT5_over_dX,n_l,n_e,inverse=True)
        else:
            output = np.zeros((len(chans),len(n_l),len(n_e)))
            for i in range(len(chans)):
                if 3 == chans[i]:
                    output[i,:,:] = copy_over_C(dBT3_over_dX,n_l,n_e)
                elif 4 == chans[i]:
                    output[i,:,:] = copy_over_C(dBT4_over_dX,n_l,n_e)
                elif 5 == chans[i]:
                    output[i,:,:] = copy_over_C(dBT5_over_dX,n_l,n_e)
                    
    return output
#
# No 12 micron channel case
#
def copy_C2(dBT3_over_dX,dBT4_over_dX,\
                n_l,n_e,chans,inverse=False,xchan=False):

    if xchan:
        if inverse:
            output = np.zeros((len(n_e),len(n_l),len(chans)))
            for i in range(len(chans)):
                if 3 == chans[i]:
                    output[:,:,i] = copy_over_C(dBT3_over_dX,n_l,n_e,inverse=True)
                elif 4 == chans[i]:
                    output[:,:,i] = copy_over_C(dBT4_over_dX,n_l,n_e,inverse=True)
                else:
                    raise Exception('chans out of range in copy_C2')
        else:
            output = np.zeros((len(n_l),len(n_e),len(chans)))
            for i in range(len(chans)):
                if 3 == chans[i]:
                    output[:,:,i] = copy_over_C(dBT3_over_dX,n_l,n_e)
                elif 4 == chans[i]:
                    output[:,:,i] = copy_over_C(dBT4_over_dX,n_l,n_e)
                else:
                    raise Exception('chans out of range in copy_C2')
    else:
        if inverse:
            output = np.zeros((len(chans),len(n_e),len(n_l)))
            for i in range(len(chans)):
                if 3 == chans[i]:
                    output[i,:,:] = copy_over_C(dBT3_over_dX,n_l,n_e,inverse=True)
                elif 4 == chans[i]:
                    output[i,:,:] = copy_over_C(dBT4_over_dX,n_l,n_e,inverse=True)
                else:
                    raise Exception('chans out of range in copy_C2')
        else:
            output = np.zeros((len(chans),len(n_l),len(n_e)))
            for i in range(len(chans)):
                if 3 == chans[i]:
                    output[i,:,:] = copy_over_C(dBT3_over_dX,n_l,n_e)
                elif 4 == chans[i]:
                    output[i,:,:] = copy_over_C(dBT4_over_dX,n_l,n_e)
                else:
                    raise Exception('chans out of range in copy_C2')

    return output

#
# Visible channel case - only from CS data
#

#
# Do copy in a block and deal with different channel selection
#
def copy_C3_vis(dRe1_over_dCS,dRe2_over_dCS,dRe3a_over_dCS,\
                n_l,n_e,chans,inverse=False,xchan=False):

    if xchan:
        if inverse:
            output = np.zeros((len(n_e),len(n_l),len(chans)))
            for i in range(len(chans)):
                if 0 == chans[i]:
                    output[:,:,i] = copy_over_C(dRe1_over_dCS,n_l,n_e,inverse=True)
                elif 1 == chans[i]:
                    output[:,:,i] = copy_over_C(dRe2_over_dCS,n_l,n_e,inverse=True)
                elif 2 == chans[i]:
                    output[:,:,i] = copy_over_C(dRe3a_over_dCS,n_l,n_e,inverse=True)
        else:
            output = np.zeros((len(n_l),len(n_e),len(chans)))
            for i in range(len(chans)):
                if 0 == chans[i]:
                    output[:,:,i] = copy_over_C(dRe1_over_dCS,n_l,n_e)
                elif 1 == chans[i]:
                    output[:,:,i] = copy_over_C(dRe2_over_dCS,n_l,n_e)
                elif 2 == chans[i]:
                    output[:,:,i] = copy_over_C(dRe3a_over_dCS,n_l,n_e)
    else:
        if inverse:
            output = np.zeros((len(chans),len(n_e),len(n_l)))
            for i in range(len(chans)):
                if 0 == chans[i]:
                    output[i,:,:] = copy_over_C(dRe1_over_dCS,n_l,n_e,inverse=True)
                elif 1 == chans[i]:
                    output[i,:,:] = copy_over_C(dRe2_over_dCS,n_l,n_e,inverse=True)
                elif 2 == chans[i]:
                    output[i,:,:] = copy_over_C(dRe3a_over_dCS,n_l,n_e,inverse=True)
        else:
            output = np.zeros((len(chans),len(n_l),len(n_e)))
            for i in range(len(chans)):
                if 0 == chans[i]:
                    output[i,:,:] = copy_over_C(dRe1_over_dCS,n_l,n_e)
                elif 1 == chans[i]:
                    output[i,:,:] = copy_over_C(dRe2_over_dCS,n_l,n_e)
                elif 2 == chans[i]:
                    output[i,:,:] = copy_over_C(dRe3a_over_dCS,n_l,n_e)
                    
    return output
#
# No 1.6 micron channel case
#
def copy_C2_vis(dRe1_over_dCS,dRe2_over_dCS,\
                n_l,n_e,chans,inverse=False,xchan=False):

    if xchan:
        if inverse:
            output = np.zeros((len(n_e),len(n_l),len(chans)))
            for i in range(len(chans)):
                if 0 == chans[i]:
                    output[:,:,i] = copy_over_C(dRe1_over_dCS,n_l,n_e,inverse=True)
                elif 1 == chans[i]:
                    output[:,:,i] = copy_over_C(dRe2_over_dCS,n_l,n_e,inverse=True)
                else:
                    raise Exception('chans out of range in copy_C2_vis')
        else:
            output = np.zeros((len(n_l),len(n_e),len(chans)))
            for i in range(len(chans)):
                if 0 == chans[i]:
                    output[:,:,i] = copy_over_C(dRe1_over_dCS,n_l,n_e)
                elif 1 == chans[i]:
                    output[:,:,i] = copy_over_C(dRe2_over_dCS,n_l,n_e)
                else:
                    raise Exception('chans out of range in copy_C2_vis')
    else:
        if inverse:
            output = np.zeros((len(chans),len(n_e),len(n_l)))
            for i in range(len(chans)):
                if 0 == chans[i]:
                    output[i,:,:] = copy_over_C(dRe1_over_dCS,n_l,n_e,inverse=True)
                elif 1 == chans[i]:
                    output[i,:,:] = copy_over_C(dRe2_over_dCS,n_l,n_e,inverse=True)
                else:
                    raise Exception('chans out of range in copy_C2_vis')
        else:
            output = np.zeros((len(chans),len(n_l),len(n_e)))
            for i in range(len(chans)):
                if 0 == chans[i]:
                    output[i,:,:] = copy_over_C(dRe1_over_dCS,n_l,n_e)
                elif 1 == chans[i]:
                    output[i,:,:] = copy_over_C(dRe2_over_dCS,n_l,n_e)
                else:
                    raise Exception('chans out of range in copy_C2_vis')

    return output

#
# Check for bad data at the array level
#
def check_for_bad_data_TL(array):

    test1 = ~np.isfinite(array)
    gd = np.isfinite(array)
    test2 = (array[gd] < -1e20)
    if np.ma.is_masked(array):
        test3 = ~array.mask 
        error = np.sum(test1) + np.sum(test2) + np.sum(test3)
    else:
        error = np.sum(test1) + np.sum(test2)
    if 0 == error:
        return False
    else:
        return True

#
# Check for bad data to mask out in CURUC routines
#
def check_bad_data_5(lines,elems,quality,dBT3_over_dT,dBT4_over_dT,\
                         dBT5_over_dT,dBT3_over_dCS,\
                         dBT4_over_dCS,dBT5_over_dCS,\
                         dBT3_over_dCICT,dBT4_over_dCICT,\
                         dBT5_over_dCICT,single=False):

    index = np.zeros(len(lines),dtype=np.int32)
    
    j=0
    for i in range(len(lines)): 
        if not single:
            test = (quality[lines[i]] > 0) or \
                check_for_bad_data_TL(dBT3_over_dT[lines[i],elems]) or \
                check_for_bad_data_TL(dBT4_over_dT[lines[i],elems]) or \
                check_for_bad_data_TL(dBT5_over_dT[lines[i],elems]) or \
                check_for_bad_data_TL(dBT3_over_dCS[lines[i],elems]) or \
                check_for_bad_data_TL(dBT4_over_dCS[lines[i],elems]) or \
                check_for_bad_data_TL(dBT5_over_dCS[lines[i],elems]) or \
                check_for_bad_data_TL(dBT3_over_dCICT[lines[i],elems]) or \
                check_for_bad_data_TL(dBT4_over_dCICT[lines[i],elems]) or \
                check_for_bad_data_TL(dBT5_over_dCICT[lines[i],elems])
        else:
            test = (quality[lines[i]] > 0) or \
                check_for_bad_data_TL(dBT3_over_dT[lines[i],elems]) or \
                check_for_bad_data_TL(dBT4_over_dT[lines[i],elems]) or \
                check_for_bad_data_TL(dBT5_over_dT[lines[i],elems])
        if test:
            index[j] = i
            j=j+1

    index = index[0:j]

    return index

def check_bad_data_no5(lines,elems,quality,dBT3_over_dT,dBT4_over_dT,\
                         dBT3_over_dCS,\
                         dBT4_over_dCS,\
                         dBT3_over_dCICT,dBT4_over_dCICT,\
                         single=False):

    index = np.zeros(len(lines),dtype=np.int32)
    
    j=0
    for i in range(len(lines)): 
        if not single:
            test = (quality[lines[i]] > 0) or \
                check_for_bad_data_TL(dBT3_over_dT[lines[i],elems]) or \
                check_for_bad_data_TL(dBT4_over_dT[lines[i],elems]) or \
                check_for_bad_data_TL(dBT3_over_dCS[lines[i],elems]) or \
                check_for_bad_data_TL(dBT4_over_dCS[lines[i],elems]) or \
                check_for_bad_data_TL(dBT3_over_dCICT[lines[i],elems]) or \
                check_for_bad_data_TL(dBT4_over_dCICT[lines[i],elems])
        else:
            test = (quality[lines[i]] > 0) or \
                check_for_bad_data_TL(dBT3_over_dT[lines[i],elems]) or \
                check_for_bad_data_TL(dBT4_over_dT[lines[i],elems])
        if test:
            index[j] = i
            j=j+1

    index = index[0:j]

    return index

def run_CURUC(data,inchans,line_skip=5,elem_skip=25):

    #
    # Check chans in right order - note chans goes from 0 to 5
    #
    chans = np.sort(inchans)
    #
    # Make two sets of chans - one visible, one IR
    #
    gd = (chans <= 2)
    vis_chans = chans[gd]
    ir_start = vis_chans[-1]+1
    gd = (chans >= 3)
    ir_chans = chans[gd]
    #
    # Allocate CURUC arrays
    #
    # Note three systematic effects considered:
    #
    #    1) Uncertainty in the ICT temperature
    #    2) Uncertainty in the space view counts
    #    3) Uncertainty in the ICT view counts
    #
    # Effect 1-3 are all spatially correlated due to the averaging
    #
    # Effect 1 will cause a cross channel correlation
    #
    # Just noise for the independent effect
    #
    nlines = data.ch1.shape[0]
    nelems = data.ch1.shape[1]
    neffects_s = 3
    neffects_i = 1
    nchans = len(chans)
    R_xelem_s, R_xline_s, R_xchan_i, R_xchan_s,\
        U_xelem_s, U_xline_s, U_xchan_s, U_xchan_i,\
        C_xelem_s, C_xline_s, C_xchan_s, C_xchan_i, coords = \
        met.allocate_curuc(nchans,nlines,nelems,neffects_s,neffects_i,\
                               line_skip,elem_skip)

    # Fill CURUC arrays
    #
    # First define R (correlation) matrices
    #
    # Cross element
    R_xelem_s.values[...] = 1  # Fully systematic across elements
    
    # Cross line - block diagonal
    # Uses temp array for spatial correlation which is the copied
    # to all channels etc.
    Temp_array = np.zeros((len(R_xline_s.coords['n_l']),\
                               len(R_xline_s.coords['n_l'])))

    N = data.spatial_correlation_scale
    Nsize = N//line_skip
    step_size = line_skip*1./N

    for step in range(-Nsize,Nsize+1):
        diagonal = np.diagonal(Temp_array,step)
        diagonal.setflags(write=True)
        diagonal[:] = 1.-np.abs(step)*step_size

    R_xline_s.values[:,:,:,:,:] = \
        Temp_array[np.newaxis,np.newaxis,np.newaxis,:,:]

    # Cross channel (independent)
    # effects, lines, elements, chan, chan
    R_xchan_i.values[0,:,:,:,:] = 1

    # Cross channel (systematic)
    # Set cross channel only for effect which is cross channel
    # Other effects have Identity matrix
    # Take into account visible channels which are not correlated with
    # the first effect (Tict)
    # Diagnonal for the vis chans, full matrix for the IR chans
    for i in range(ir_start-1):
        R_xchan_s.values[0,:,:,i,i] = 1
    R_xchan_s.values[0,:,:,ir_start:,ir_start:] = 1
    Temp_array = np.zeros((len(R_xline_s.coords['n_c']),\
                                   len(R_xline_s.coords['n_c'])))
    Temp_array[:,:] = np.identity(len(R_xline_s.coords['n_c']))
    R_xchan_s.values[1,:,:,:,:] = Temp_array[np.newaxis,np.newaxis,:,:]
    R_xchan_s.values[2,:,:,:,:] = Temp_array[np.newaxis,np.newaxis,:,:]

    # Map uncertainties (which are constant in the T, Counts space).
    # All these are averaged over the smoothing scale (+/-)
    # Elements
    # n_c, n_s, n_l, n_e
    #
    # First one is the ICT temperature one
    U_xelem_s.values[ir_start:,0,:,:] = data.ICT_Temperature_Uncertainty

    # Second is the Space view one
    for i in range(len(chans)):
        U_xelem_s.values[i,1,:,:] = data.cal_cnts_noise[chans[i]]
        
    # Third is the ICT view one
    for i in range(ir_start,len(chans)):
        U_xelem_s.values[i,2,:,:] = data.cal_cnts_noise[chans[i]]

    # Cross line case
    # Elements
    # n_c, n_s, n_e, n_l
    #
    # First one is the ICT temperature one
    U_xline_s.values[:,0,:,:] = data.ICT_Temperature_Uncertainty

    # Second is the Space view one
    for i in range(len(chans)):
        U_xline_s.values[i,1,:,:] = data.cal_cnts_noise[chans[i]]

    # Third is the ICT view one
    for i in range(len(chans)):
        U_xline_s.values[i,2,:,:] = data.cal_cnts_noise[chans[i]]

    # Uncertainties for channel to channel (systematic)
    # Elements
    # n_l, n_e, n_s, n_c
    #
    # Only ICT temperature is cross channel
    #
    # But in CURUC all are considered to match to stored uncertainties
    # in FCDR file
    #
    # Note have to remember first 2/3 channels are vis chans so no
    # temperature uncertainty
    #
    U_xchan_s.values[:,:,0,ir_start:] = data.ICT_Temperature_Uncertainty
    for i in range(len(chans)):
        U_xchan_s.values[:,:,1,i] = data.cal_cnts_noise[chans[i]]
    # Only IR channels as this is noise on ICT counts
    for i in range(len(chans)):
        if chans[i] >= 3:
            U_xchan_s.values[:,:,2,i] = data.cal_cnts_noise[chans[i]]

    # Uncertainties for channel to channel (independent) - zero
    # Note single pixel noise here
    for i in range(len(chans)):
        U_xchan_i.values[:,:,:,i] = data.cnts_noise[chans[i]]
        
    # Map sensitivities using coordinates
    # This has to use the coordinates to make correctly as the sensitivities
    # change pixel to pixel
    # 3.7 micron channel
    if data.ch5_there:
        C_xelem_s.values[ir_start:,0,:,:] = \
            copy_C3(data.dBT3_over_dT,data.dBT4_over_dT,data.dBT5_over_dT,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        ir_chans)
        C_xelem_s.values[ir_start:,1,:,:] = \
            copy_C3(data.dBT3_over_dCS,data.dBT4_over_dCS,data.dBT5_over_dCS,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        ir_chans)
        C_xelem_s.values[ir_start:,2,:,:] = \
            copy_C3(data.dBT3_over_dCICT,data.dBT4_over_dCICT,\
                        data.dBT5_over_dCICT,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        ir_chans)

        C_xline_s.values[ir_start:,0,:,:] = \
            copy_C3(data.dBT3_over_dT,data.dBT4_over_dT,data.dBT5_over_dT,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        ir_chans,inverse=True)
        C_xline_s.values[ir_start:,1,:,:] = \
            copy_C3(data.dBT3_over_dCS,data.dBT4_over_dCS,data.dBT5_over_dCS,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        ir_chans,inverse=True)
        C_xline_s.values[ir_start:,2,:,:] = \
            copy_C3(data.dBT3_over_dCICT,data.dBT4_over_dCICT,\
                        data.dBT5_over_dCICT,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        ir_chans,inverse=True)

        # Now X channel
        C_xchan_s.values[0,:,:,ir_start:] = \
            copy_C3(data.dBT3_over_dT,data.dBT4_over_dT,data.dBT5_over_dT,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        ir_chans,xchan=True)
        C_xchan_s.values[1,:,:,ir_start:] = \
            copy_C3(data.dBT3_over_dCS,data.dBT4_over_dCS,data.dBT5_over_dCS,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        ir_chans,xchan=True)
        C_xchan_s.values[2,:,:,ir_start:] = \
            copy_C3(data.dBT3_over_dCICT,data.dBT4_over_dCICT,\
                        data.dBT5_over_dCICT,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        ir_chans,xchan=True)

        #
        # Independent effects
        #
        C_xchan_i.values[0,:,:,ir_start:] = \
            copy_C3(data.dBT3_over_dT,data.dBT4_over_dT,data.dBT5_over_dT,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        ir_chans,xchan=True)

    else: # Don't have channel 5
        C_xelem_s.values[:,0,:,ir_start:] = \
            copy_C2(data.dBT3_over_dT,data.dBT4_over_dT,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        ir_chans)
        C_xelem_s.values[:,1,:,ir_start:] = \
            copy_C2(data.dBT3_over_dCS,data.dBT4_over_dCS,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        ir_chans)
        C_xelem_s.values[:,2,:,ir_start:] = \
            copy_C2(data.dBT3_over_dCICT,data.dBT4_over_dCICT,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        ir_chans)

        C_xline_s.values[:,0,:,ir_start:] = \
            copy_C2(data.dBT3_over_dT,data.dBT4_over_dT,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        ir_chans,inverse=True)
        C_xline_s.values[:,1,:,ir_start:] = \
            copy_C2(data.dBT3_over_dCS,data.dBT4_over_dCS,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        ir_chans,inverse=True)
        C_xline_s.values[:,2,:,ir_start:] = \
            copy_C2(data.dBT3_over_dCICT,data.dBT4_over_dCICT,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        ir_chans,inverse=True)

        # Now X channel
        C_xchan_s.values[0,:,:,ir_start:] = \
            copy_C2(data.dBT3_over_dT,data.dBT4_over_dT,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        ir_chans,xchan=True)
        C_xchan_s.values[1,:,:,ir_start:] = \
            copy_C2(data.dBT3_over_dCS,data.dBT4_over_dCS,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        ir_chans,xchan=True)
        C_xchan_s.values[2,:,:,ir_start:] = \
            copy_C2(data.dBT3_over_dCICT,data.dBT4_over_dCICT,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        ir_chans,xchan=True)

        C_xchan_i.values[0,:,:,ir_start:] = \
            copy_C2(data.dBT3_over_dT,data.dBT4_over_dT,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        ir_chans,xchan=True)

    #
    # Visible channel case now - only for CS averaging no X channel
    #
    if data.ch3a_there:
        C_xelem_s.values[0:ir_start,1,:,:] = \
            copy_C3_vis(data.dRe1_over_dCS,data.dRe2_over_dCS,data.dRe3a_over_dCS,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        vis_chans)

        C_xline_s.values[0:ir_start,1,:,:] = \
            copy_C3_vis(data.dRe1_over_dCS,data.dRe1_over_dCS,data.dRe3a_over_dCS,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        vis_chans,inverse=True)

    else: # Don't have channel 1.6 and no X channel
        C_xelem_s.values[0:ir_start,1,:,:] = \
            copy_C2_vis(data.dRe1_over_dCS,data.dRe2_over_dCS,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        vis_chans)

        C_xline_s.values[0:ir_start,1,:,:] = \
            copy_C2_vis(data.dRe1_over_dCS,data.dRe2_over_dCS,\
                        C_xelem_s.coords['n_l'],C_xelem_s.coords['n_e'],\
                        vis_chans,inverse=True)

    # make mask data (xarray bool)
    # Keep all channels (dtype='?' means boolean)
    mask1 = xarray.DataArray(np.zeros(nchans,dtype='?'),dims=['n_c'])

    # Setup mask for bad lines
    mask2 = xarray.DataArray(np.zeros(len(R_xline_s.coords['n_l']),dtype='?'),\
                                 dims=['n_l'])
    
    #
    # Check for bad scan lines and lines with bad TLs
    #
    if data.ch5_there:
        index = check_bad_data_5(R_xline_s.coords['n_l'].values,\
                                     R_xline_s.coords['n_e'].values,\
                                     data.scan_qual,\
                                     data.dBT3_over_dT,\
                                     data.dBT4_over_dT,\
                                     data.dBT5_over_dT,\
                                     data.dBT3_over_dCS,\
                                     data.dBT4_over_dCS,\
                                     data.dBT5_over_dCS,\
                                     data.dBT3_over_dCICT,\
                                     data.dBT4_over_dCICT,\
                                     data.dBT5_over_dCICT)
    else:
        index = check_bad_data_no5(R_xline_s.coords['n_l'].values,\
                                       R_xline_s.coords['n_e'].values,\
                                       data.scan_qual,\
                                       data.dBT3_over_dT,\
                                       data.dBT4_over_dT,\
                                       data.dBT3_over_dCS,\
                                       data.dBT4_over_dCS,\
                                       data.dBT3_over_dCICT,\
                                       data.dBT4_over_dCICT)
    mask2.values[index] = True
    
    # Now run CURUC
    ny = int(data.ch1.shape[0]/line_skip)
    if 0 == ny:
        ny = 1
    nx = int(data.ch1.shape[1]/elem_skip)
    if 0 == nx:
        nx = 1
#    np.seterr(divide="raise",over="raise",invalid="raise",under="ignore")
    xline_length, xelem_length, xchan_corr_i, xchan_corr_s, xl_all, xe_all = \
        met.apply_curuc(R_xelem_s, R_xline_s, R_xchan_i, R_xchan_s,\
                            U_xelem_s, U_xline_s, U_xchan_s, U_xchan_i,\
                            C_xelem_s, C_xline_s, C_xchan_s, C_xchan_i,\
                            coords,mask1,mask2,return_vectors=True,\
                            interpolate_lengths=True,\
                            cutoff_l=ny,cutoff_e=nx)    

    return xline_length, xelem_length, xchan_corr_i, xchan_corr_s, \
        xl_all, xe_all

#
# Get SRF information for a given AVHRR
#
def get_srf(noaa,allchans):

    srf_dir = \
        '/group_workspaces/cems2/fiduceo/Users/jmittaz/FCDR/Mike/FCDR_AVHRR/SRF/data/'

    #
    # Remove visible channels as they are not controlled by a FIDUCEO
    # process and therefore we can't know what was used
    #
    gd = (allchans >= 3)
    chans = allchans[gd]

    #
    # Get filenames
    #
    if noaa == 'TIROSN':
        srf_file_wave = srf_dir+'tiros-n_wave.dat'
        srf_file_srf = srf_dir+'tiros-n_srf.dat'
        lut_radiance = srf_dir+'tiros-n_rad.dat'
        lut_bt = srf_dir+'tiros-n_bt.dat'
        inchans = [3,4]
    elif noaa == 'NOAA06':
        srf_file_wave = srf_dir+'noaa06_wave.dat'
        srf_file_srf = srf_dir+'noaa06_srf.dat'
        lut_radiance = srf_dir+'noaa06_rad.dat'
        lut_bt = srf_dir+'noaa06_bt.dat'
        inchans = [3,4]
    elif noaa == 'NOAA07':
        srf_file_wave = srf_dir+'noaa07_wave.dat'
        srf_file_srf = srf_dir+'noaa07_srf.dat'
        lut_radiance = srf_dir+'noaa07_rad.dat'
        lut_bt = srf_dir+'noaa07_bt.dat'
        inchans = [3,4,5]
    elif noaa == 'NOAA08':
        srf_file_wave = srf_dir+'noaa08_wave.dat'
        srf_file_srf = srf_dir+'noaa08_srf.dat'
        lut_radiance = srf_dir+'noaa08_rad.dat'
        lut_bt = srf_dir+'noaa08_bt.dat'
        inchans = [3,4]
    elif noaa == 'NOAA09':
        srf_file_wave = srf_dir+'noaa09_wave.dat'
        srf_file_srf = srf_dir+'noaa09_srf.dat'
        lut_radiance = srf_dir+'noaa09_rad.dat'
        lut_bt = srf_dir+'noaa09_bt.dat'
        inchans = [3,4,5]
    elif noaa == 'NOAA10':
        srf_file_wave = srf_dir+'noaa10_wave.dat'
        srf_file_srf = srf_dir+'noaa10_srf.dat'
        lut_radiance = srf_dir+'noaa10_rad.dat'
        lut_bt = srf_dir+'noaa10_bt.dat'
        inchans = [3,4]
    elif noaa == 'NOAA11':
        srf_file_wave = srf_dir+'noaa11_wave.dat'
        srf_file_srf = srf_dir+'noaa11_srf.dat'
        lut_radiance = srf_dir+'noaa11_rad.dat'
        lut_bt = srf_dir+'noaa11_bt.dat'
        inchans = [3,4,5]
    elif noaa == 'NOAA12':
        srf_file_wave = srf_dir+'noaa12_wave.dat'
        srf_file_srf = srf_dir+'noaa12_srf.dat'
        lut_radiance = srf_dir+'noaa12_rad.dat'
        lut_bt = srf_dir+'noaa12_bt.dat'
        inchans = [3,4,5]
    elif noaa == 'NOAA14':
        srf_file_wave = srf_dir+'noaa14_wave.dat'
        srf_file_srf = srf_dir+'noaa14_srf.dat'
        lut_radiance = srf_dir+'noaa14_rad.dat'
        lut_bt = srf_dir+'noaa14_bt.dat'
        inchans = [3,4,5]
    elif noaa == 'NOAA15':
        srf_file_wave = srf_dir+'noaa15_wave.dat'
        srf_file_srf = srf_dir+'noaa15_srf.dat'
        lut_radiance = srf_dir+'noaa15_rad.dat'
        lut_bt = srf_dir+'noaa15_bt.dat'
        inchans = [3,4,5]
    elif noaa == 'NOAA16':
        srf_file_wave = srf_dir+'noaa16_wave.dat'
        srf_file_srf = srf_dir+'noaa16_srf.dat'
        lut_radiance = srf_dir+'noaa16_rad.dat'
        lut_bt = srf_dir+'noaa16_bt.dat'
        inchans = [3,4,5]
    elif noaa == 'NOAA17':
        srf_file_wave = srf_dir+'noaa17_wave.dat'
        srf_file_srf = srf_dir+'noaa17_srf.dat'
        lut_radiance = srf_dir+'noaa17_rad.dat'
        lut_bt = srf_dir+'noaa17_bt.dat'
        inchans = [3,4,5]
    elif noaa == 'NOAA18':
        srf_file_wave = srf_dir+'noaa18_wave.dat'
        srf_file_srf = srf_dir+'noaa18_srf.dat'
        lut_radiance = srf_dir+'noaa18_rad.dat'
        lut_bt = srf_dir+'noaa18_bt.dat'
        inchans = [3,4,5]
    elif noaa == 'NOAA19':
        srf_file_wave = srf_dir+'noaa19_wave.dat'
        srf_file_srf = srf_dir+'noaa19_srf.dat'
        lut_radiance = srf_dir+'noaa19_rad.dat'
        lut_bt = srf_dir+'noaa19_bt.dat'
        inchans = [3,4,5]
    elif noaa == 'METOPA':
        srf_file_wave = srf_dir+'metopa_wave.dat'
        srf_file_srf = srf_dir+'metopa_srf.dat'
        lut_radiance = srf_dir+'metopa_rad.dat'
        lut_bt = srf_dir+'metopa_bt.dat'
        inchans = [3,4,5]
    elif noaa == 'METOPB':
        srf_file_wave = srf_dir+'metopb_wave.dat'
        srf_file_srf = srf_dir+'metopb_srf.dat'
        lut_radiance = srf_dir+'metopb_rad.dat'
        lut_bt = srf_dir+'metopb_bt.dat'
        inchans = [3,4,5]
    else:
        raise Exception('Cannot find noaa name for SRF')

    #
    # Read in SRF values themselves
    #    
    wave = np.loadtxt(srf_file_wave)
    srf = np.loadtxt(srf_file_srf)

    #
    # Read in lookup tables
    #
    radiance = np.loadtxt(lut_radiance)
    bt = np.loadtxt(lut_bt)

    #
    # Only select channels that are being used
    #
    start,end = inchans.index(chans[0]),len(chans)

    out_wave = wave[start:end,:]
    out_srf = srf[start:end,:]
    out_radiance = radiance[start:end,:]
    out_radiance = out_radiance.transpose()
    out_bt = bt[start:end,:]
    out_bt = out_bt.transpose()

    # Set 0's in srf to NaN
    gd = (wave == 0)
    out_wave[gd] = float('nan')
    out_srf[gd] = float('nan')

    return out_wave,out_srf,out_radiance,out_bt

#
# Calculate CURUC etc. and output file. Note changes behaviour
# dependent on channel set
#
def main_outfile(data,ch3a_version,fileout='None',split=False,gbcs_l1c=False):

    # Run CURUC to get CURUC values (lenths, vectors and chan cross 
    # correlations)
    if data.noaa_string in ['NOAA06','NOAA08','NOAA10']:
        # Setup CURUC for 2 channel IR AVHRR
        chans = np.array([0,1,3,4],dtype=np.int8)
    elif data.noaa_string in ['NOAA07','NOAA09','NOAA11','NOAA12','NOAA14']:
        # Setup CURUC for 3 channel IR AVHRR
        chans = np.array([0,1,3,4,5],dtype=np.int8)
    elif data.noaa_string in ['NOAA15','NOAA16','NOAA17','NOAA18',\
                                  'NOAA19','METOPA','METOPB','METOPC']:
        # Setup CURUC for 2 or 3 channel IR AVHRR dependnet on ch3a
        if ch3a_version:
            chans = np.array([0,1,2,4,5],dtype=np.int8)
        else:
            chans = np.array([0,1,3,4,5],dtype=np.int8)
    else:
        raise Exception('noaa_string not found')
    xline_length, xelem_length, xchan_corr_i, xchan_corr_s, \
        xl_all, xe_all = run_CURUC(data,chans,line_skip=5,elem_skip=25)

    gd = (xl_all.values > 0.)
    corr_l = xl_all.values[gd]
    corr_l = np.append(corr_l,np.array([0.]))

    # Get SRF and lookup tables
    srf_x,srf_y,lut_rad,lut_bt = get_srf(data.noaa_string,chans)


# MT: 09-11-2017: define sensor specific channel_correlation_matrix (ccm)
    noch3a=False
    noch5=False
    if data.noaa_string in ['NOAA06','NOAA08','NOAA10']:
        # Run CURUC for 2 channel IR AVHRR
        inS_s = np.array([[xchan_corr_s.values[0,0],\
                             xchan_corr_s.values[0,1]],\
                            [xchan_corr_s.values[1,0],\
                                 xchan_corr_s.values[1,1]]])
        inS_i = np.array([[1.,0.,],[0.,1.]])
        start=3
        stop=5
        noch3a=True
        noch5=True
    elif data.noaa_string in ['NOAA07','NOAA09','NOAA11','NOAA12','NOAA14']:
        inS_s = np.array([[xchan_corr_s.values[0,0],\
                             xchan_corr_s.values[0,1],\
                             xchan_corr_s.values[0,2]],\
                            [xchan_corr_s.values[1,0],\
                                 xchan_corr_s.values[1,1],\
                                 xchan_corr_s.values[1,2]],\
                            [xchan_corr_s.values[2,0],\
                                 xchan_corr_s.values[2,1],\
                                 xchan_corr_s.values[2,2]]])
        inS_i = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
        start=3
        stop=6
        noch3a=True
    elif data.noaa_string in ['NOAA15','NOAA16','NOAA17','NOAA18','NOAA19','METOPA','METOPB','METOPC']:
        if ch3a_version:
            inS_i = np.array([[0.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
            inS_s = np.diag([0,\
                               xchan_corr_s.values[0,1],\
                               xchan_corr_s.values[0,2],\
                               xchan_corr_s.values[1,0],\
                               xchan_corr_s.values[1,1],\
                               xchan_corr_s.values[1,2],\
                               xchan_corr_s.values[2,0],\
                               xchan_corr_s.values[2,1],\
                               xchan_corr_s.values[2,2]])
        else:
            inS_i = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
            inS_s = np.array([[0,0,0],\
                                [0,xchan_corr_s.values[0,0],\
                                     xchan_corr_s.values[0,1]],\
                                [0,xchan_corr_s.values[1,0],\
                                     xchan_corr_s.values[1,1]]])
        start=3
        stop=6

    S_s = np.identity(6)
    S_i = np.identity(6)
    if noch3a:
        S_s[2,2]=0.
        S_i[2,2]=0.
    if noch5:
        S_s[5,5]=0.
        S_i[5,5]=0.
    S_s[start:stop,start:stop] = inS_s
    S_i[start:stop,start:stop] = inS_i

    # Either write L1C with chan covariance or FIDUCEO easy FCDR
    if gbcs_l1c:
        l1c.write_gbcs_l1c(fileout,data,S_s)
    else:
        writer = FCDRWriter()

        # get a template for sensor name in EASY format, supply product height
        # The scan-width is set automatically
        # Update now gives srf_size for SRF
        #                  corr_dx for correlation vector for elements
        #                  corr_dy for correlation vector scanline direction
        #                  lut_size for radiance conversion to BT lokkup table
        dataset = writer.createTemplateEasy("AVHRR", data.ny, \
                                                srf_size=srf_x.shape[1],\
                                                corr_dx=data.lat.shape[1],\
                                                corr_dy=len(corr_l),\
                                                lut_size=lut_rad.shape[0])

        # set some mandatory global attributes. Writing will fail if not all of them are filled
        dataset.attrs["institution"] = "University of Reading"
        dataset.attrs["title"] = "pre-B version of AVHRR Fundamental Climate Data Records"
        dataset.attrs["source"] = data.sources
        dataset.attrs["history"] = ""
        dataset.attrs["references"] = "CDF_FCDR_File Spec"
        dataset.attrs["comment"] = "This version is a pre-B one and aims at showing the kind of uncertainties we are aiming to deliver within FIDUCEO. The values are not final ones and should not be used for science purposes."
        dataset.attrs['sensor'] = "AVHRR"
        dataset.attrs['platform'] = data.noaa_string
        dataset.attrs['software_version'] = data.version
        # write real data to the variables. All variables initially contain "_FillValue".
        # Not writing to the whole array is completely OK
        dataset.variables["latitude"].data = data.lat
        dataset.variables["longitude"].data = data.lon
        dataset.variables["Time"].data = data.time
        dataset.variables["satellite_zenith_angle"].data = data.satza
        dataset.variables["solar_zenith_angle"].data = data.solza
        dataset.variables["relative_azimuth_angle"].data = data.relaz
        dataset.variables["Ch1"].data = data.ch1
        dataset.variables["Ch2"].data = data.ch2
        if data.ch3a_there:
            dataset.variables["Ch3a"].data = data.ch3a
        dataset.variables["Ch3b"].data = data.ch3b
        dataset.variables["Ch4"].data = data.ch4
        if data.ch5_there:
            dataset.variables["Ch5"].data = data.ch5
# MT: 30-10-2017: uncertainty variable name change in line with FCDR-CDR file format spec fv1.1.1
        dataset.variables["u_independent_Ch1"].data = data.u_random_ch1
        dataset.variables["u_independent_Ch2"].data = data.u_random_ch2
        if data.ch3a_there:
            dataset.variables["u_independent_Ch3a"].data = data.u_random_ch3a
        dataset.variables["u_independent_Ch3b"].data = data.u_random_ch3b
        dataset.variables["u_independent_Ch4"].data = data.u_random_ch4
        if data.ch5_there:
            dataset.variables["u_independent_Ch5"].data = data.u_random_ch5   
        dataset.variables["u_structured_Ch1"].data = data.u_non_random_ch1
        dataset.variables["u_structured_Ch2"].data = data.u_non_random_ch2
        if data.ch3a_there:
            dataset.variables["u_structured_Ch3a"].data = data.u_non_random_ch3a
        dataset.variables["u_structured_Ch3b"].data = data.u_non_random_ch3b
        dataset.variables["u_structured_Ch4"].data = data.u_non_random_ch4
        if data.ch5_there:
            dataset.variables["u_structured_Ch5"].data = data.u_non_random_ch5
        dataset.variables["u_common_Ch1"].data = data.u_common_ch1
        dataset.variables["u_common_Ch2"].data = data.u_common_ch2
        if data.ch3a_there:
            dataset.variables["u_common_Ch3a"].data = data.u_common_ch3a
        dataset.variables["u_common_Ch3b"].data = data.u_common_ch3b
        dataset.variables["u_common_Ch4"].data = data.u_common_ch4
        if data.ch5_there:
            dataset.variables["u_common_Ch5"].data = data.u_common_ch5
# MT: 18-10-2017: added quality flag fields
        dataset.variables["quality_scanline_bitmask"].data = data.scan_qual
        dataset.variables["quality_channel_bitmask"].data = data.chan_qual

# Update for reader version 1.1.5
        dataset.variables["channel_correlation_matrix_structured"].data = S_s
        dataset.variables["channel_correlation_matrix_independent"].data = S_i
        elem_corr = np.zeros((data.lat.shape[1],stop-start))
        elem_corr[:,:] = 1.
        dataset.variables["cross_element_correlation_coefficients"].\
            data[:,start:stop] = elem_corr
        line_corr = np.zeros((corr_l.shape[0],stop-start))
        for i in range(line_corr.shape[1]):
            line_corr[:,i] = corr_l
        dataset.variables["cross_line_correlation_coefficients"].\
            data[:,start:stop] = line_corr
        
# JM: 05/07/2018: Added in SRF data and rad->BT etc. lookup tables
        dataset.variables["SRF_weights"].data[start:stop,:] = srf_y
        dataset.variables["SRF_frequencies"].data[start:stop,:] = srf_x
        dataset.variables["lookup_table_BT"].data[:,start:stop] = lut_bt
        dataset.variables["lookup_table_radiance"].data[:,start:stop] = lut_rad

        # Traceability back to original Level 1B file and scanline
        dataset.variables["scanline_map_to_origl1bfile"].data[:] = \
            data.orig_scanline
        dataset.variables["scanline_origl1b"].data[:] = data.scanline

        # dump it to disk, netcdf4, medium compression
        # writing will fail when the target file already exists
        if 'None' == fileout:
            file_out = writer.create_file_name_FCDR_easy('AVHRR',data.noaa_string,\
                                                             data.date_time[0],\
                                                             data.date_time[-1],\
                                                             data.version)
        else:
            if split:
                if ch3a_version:
                    file_out = 'ch3a_'+fileout
                else:
                    file_out = 'ch3b_'+fileout
            else:
                file_out = fileout
        writer.write(dataset, file_out)

#
# Copy data into data class based on filter
#
class copy_to_data(object):

    def copy(self,data,gd):
        self.ch3a_there = data.ch3a_there
        self.ch5_there = data.ch5_there
        self.lat = data.lat[gd,:]
        self.lon = data.lon[gd,:]
        self.time = data.time[gd]
        self.satza = data.satza[gd,:]
        self.solza = data.solza[gd,:]
        self.relaz = data.relaz[gd,:]
        self.ch1 = data.ch1[gd,:]
        self.ch2 = data.ch2[gd,:]
        if self.ch3a_there:
            self.ch3a = data.ch3a[gd,:]
        self.ch3b = data.ch3b[gd,:]
        self.ch4 = data.ch4[gd,:]
        if self.ch5_there:
            self.ch5 = data.ch5[gd,:]
        self.u_random_ch1 = data.u_random_ch1[gd,:]
        self.u_random_ch2 = data.u_random_ch2[gd,:]
        if self.ch3a_there:
            self.u_random_ch3a = data.u_random_ch3a[gd,:]
        self.u_random_ch3b = data.u_random_ch3b[gd,:]
        self.u_random_ch4 = data.u_random_ch4[gd,:]
        if self.ch5_there:
            self.u_random_ch5 = data.u_random_ch5[gd,:]
        self.u_non_random_ch1 = data.u_non_random_ch1[gd,:]
        self.u_non_random_ch2 = data.u_non_random_ch2[gd,:]
        if self.ch3a_there:
            self.u_non_random_ch3a = data.u_non_random_ch3a[gd,:]
        self.u_non_random_ch3b = data.u_non_random_ch3b[gd,:]
        self.u_non_random_ch4 = data.u_non_random_ch4[gd,:]
        if self.ch5_there:
            self.u_non_random_ch5 = data.u_non_random_ch5[gd,:]
        self.dBT3_over_dT = data.dBT3_over_dT[gd,:]
        self.dBT4_over_dT = data.dBT4_over_dT[gd,:]
        self.dBT5_over_dT = data.dBT5_over_dT[gd,:]
        self.dBT3_over_dCS = data.dBT3_over_dCS[gd,:]
        self.dBT4_over_dCS = data.dBT4_over_dCS[gd,:]
        self.dBT5_over_dCS = data.dBT5_over_dCS[gd,:]
        self.dBT3_over_dCICT = data.dBT3_over_dCICT[gd,:]
        self.dBT4_over_dCICT = data.dBT4_over_dCICT[gd,:]
        self.dBT5_over_dCICT = data.dBT5_over_dCICT[gd,:]

        self.cal_cnts_noise = data.cal_cnts_noise[:]
        self.cnts_noise = data.cnts_noise[:]
        self.spatial_correlation_scale = data.spatial_correlation_scale
        self.ICT_Temperature_Uncertainty = data.ICT_Temperature_Uncertainty
        self.noaa_string = data.noaa_string

    def __init__(self,data,gd):

        self.copy(data,gd)

#
# Split out data with and without ch3a data
#
def get_split_data(data,ch3a=False):

    try:
        gd = np.zeros(data.ch4.shape[0],dtype='?')
        if ch3a:
            for i in range(len(gd)):
                gd[i] = np.any(np.isfinite(data.ch3a[i,:])) 
        else:
            for i in range(len(gd)):
                gd[i] = np.any(np.isfinite(data.ch3b[i,:])) 
        
    except:
        raise Exception('Cannot find c3a data when trying to split')

    newdata = copy_to_data(data,gd)
    
    return newdata

#
# Top level routine to output FCDR
#
def main(file_in,fileout='None'):

    data = read_netcdf(file_in)

    #
    # If we have c3a data then have to split file to 2 channel and 3 channel
    # cases within data
    #
    if data.ch3a_there:
        #
        # Have to split orbit into two to ensure CURUC works
        #
        data1 = get_split_data(data,ch3a=True)
        main_outfile(data1,ch3a_version=True,fileout=fileout,split=True)
        data2 = get_split_data(data,ch3a=False)
        main_outfile(data1,ch3a_version=False,fileout=fileout,split=True)
    else:
        main_outfile(data,ch3a_version=False,fileout=fileout)

if __name__ == "__main__":

    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()
    
    if len(args) != 1 and len(args) != 2:
        parser.error("incorrect number of arguments")
    
    if len(args) == 1:
        main(args[0])
    elif len(args) == 2:
        main(args[0],fileout=args[1])


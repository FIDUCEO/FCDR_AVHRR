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

from fiduceo.fcdr.writer.fcdr_writer import FCDRWriter
from fiduceo.fcdr.writer.templates import avhrr
import sys
import netCDF4
import numpy as np
import datetime
import xarray
from optparse import OptionParser

class read_netcdf(object):

    def add_nan_values(self,values):
        with np.errstate(invalid='ignore'):
            gd = np.isfinite(values) & (values < -1e20)
        if np.sum(gd) > 0:
            values[gd] = np.NaN

    def scale_values(self,values):
        with np.errstate(invalid='ignore'):
            gd = (values > -1e20) & np.isfinite(values)
        if np.sum(gd) > 0:
            values[gd] = values[gd]*100.

    def read_data(self,filename):

        ncid = netCDF4.Dataset(filename,'r')

        self.noaa_string = ncid.noaa_string
        self.version = ncid.version
        self.spatial_correlation_scale = ncid.spatial_correlation_scale
        year = ncid.variables['year'][:]
        month = ncid.variables['month'][:]
        day = ncid.variables['day'][:]
        hours = ncid.variables['hours'][:]
        gd = (year > 1900) & (day > 0) & (month > 0) & (hours >= 0)
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
        try:
            ch3a = ncid.variables['ch3a'][:,:]
            self.ch3a = ch3a[gd,:]
            self.ch3a_there = True
        except:
            self.ch3a_there = False
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
        scan_qual = ncid.variables['quality_scanline_bitmask'][:]
        chan_qual = ncid.variables['quality_channel_bitmask'][:,:]
        self.scan_qual = scan_qual[gd]
        self.chan_qual = chan_qual[gd,:]

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

        self.add_nan_values(self.lat)
        self.add_nan_values(self.lon)
        self.add_nan_values(self.time)
        self.add_nan_values(self.satza)
        self.add_nan_values(self.relaz)
        self.add_nan_values(self.ch1)
        self.add_nan_values(self.ch2)
        if self.ch3a_there:
            self.add_nan_values(self.ch3a)
        self.add_nan_values(self.ch3b)
        self.add_nan_values(self.ch4)
        if self.ch5_there:
            self.add_nan_values(self.ch5)
        self.add_nan_values(self.u_random_ch1)
        self.add_nan_values(self.u_random_ch2)
        if self.ch3a_there:
            self.add_nan_values(self.u_random_ch3a)
        self.add_nan_values(self.u_random_ch3b)
        self.add_nan_values(self.u_random_ch4)
        if self.ch5_there:
            self.add_nan_values(self.u_non_random_ch5)
        self.add_nan_values(self.u_non_random_ch1)
        self.add_nan_values(self.u_non_random_ch2)
        if self.ch3a_there:
            self.add_nan_values(self.u_non_random_ch3a)
        self.add_nan_values(self.u_non_random_ch3b)
        self.add_nan_values(self.u_non_random_ch4)
        if self.ch5_there:
            self.add_nan_values(self.u_non_random_ch5)
            
#        self.scale_values(self.ch1)
#        self.scale_values(self.ch2)
#        if self.ch3a_there:
#            self.scale_values(self.ch3a)
#        self.scale_values(self.u_random_ch1)
#        self.scale_values(self.u_random_ch2)
#        if self.ch3a_there:
#            self.scale_values(self.u_random_ch3a)
#        self.scale_values(self.u_non_random_ch1)
#        self.scale_values(self.u_non_random_ch2)
#        if self.ch3a_there:
#            self.scale_values(self.u_non_random_ch3a)

    def __init__(self,filename):

        self.read_data(filename)

def main(file_in,fileout='None'):
    data = read_netcdf(file_in)

    writer = FCDRWriter()

    # get a template for sensor name in EASY format, supply product height
    # The scan-width is set automatically
    dataset = writer.createTemplateEasy("AVHRR", data.ny)

    # set some mandatory global attributes. Writing will fail if not all of them are filled
    dataset.attrs["institution"] = "University of Reading"
    dataset.attrs["title"] = "pre-B version of AVHRR Fundamental Climate Data Records"
    dataset.attrs["source"] = "FIDUCEO"
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
    dataset.variables["Ch1_Ref"].data = data.ch1
    dataset.variables["Ch2_Ref"].data = data.ch2
    if data.ch3a_there:
        dataset.variables["Ch3a_Ref"].data = data.ch3a
    dataset.variables["Ch3b_Bt"].data = data.ch3b
    dataset.variables["Ch4_Bt"].data = data.ch4
    if data.ch5_there:
        dataset.variables["Ch5_Bt"].data = data.ch5
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
# MT: 18-10-2017: added quality flag fields
    dataset.variables["quality_scanline_bitmask"].data = data.scan_qual
    dataset.variables["quality_channel_bitmask"].data = data.chan_qual

# MT: 09-11-2017: define sensor specific channel_correlation_matrix (ccm)
    if data.noaa_string in ['NOAA06','NOAA08','NOAA10']:
        S = np.diag([1,1,0,1,1,0])
    elif data.noaa_string in ['NOAA07','NOAA09','NOAA11','NOAA12','NOAA14']:
        S = np.diag([1,1,0,1,1,1])
    elif data.noaa_string in ['NOAA15','NOAA16','NOAA17','NOAA18','NOAA19','METOPA','METOPB','METOPC']:
        S = np.diag([1,1,1,1,1,1])
    ccm = xarray.DataArray(S,dims=("nchan", "nchan"))
    ccm.name = "channel_correlation_matrix"
    ccm.attrs["units"] = "1"
    ccm.attrs["valid_max"] = 1
    ccm.attrs["valid_min"] = 0
    ccm.encoding["dtype"] = np.int16
#    dataset=xarray.merge([dataset,ccm])

# MT: 10-11-2017: spatial_correlation_scale (scs)
    scs = xarray.DataArray(data.spatial_correlation_scale)
    scs.name = "spatial_correlation_scale"
    scs.attrs["units"] = "1"
    scs.attrs["valid_max"] = 1000
    scs.attrs["valid_min"] = 0
    scs.encoding["dtype"] = np.int16
# MT: 23-11-2017: store and reset global variables broken by xarray.merge [thanks, Gerrit Holl]
    attr = dataset.attrs
    dataset=xarray.merge([dataset,ccm,scs])
    dataset.attrs = attr

    #avhrr.AVHRR._create_channel_refl_variable(12835, "Channel 6 Reflectance") 
    # dump it to disk, netcdf4, medium compression
    # writing will fail when the target file already exists
    if 'None' == fileout:
        file_out = writer.create_file_name_FCDR_easy('AVHRR',data.noaa_string,\
                                                         data.date_time[0],\
                                                         data.date_time[-1],\
                                                         data.version)
    else:
        file_out = fileout
    writer.write(dataset, file_out)

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


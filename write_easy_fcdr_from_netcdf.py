# * Copyright (C) 2017 J.Mittaz University of Reading
# * This code was developed for the EC project “Fidelity and Uncertainty in   
# * Climate Data Records from Earth Observations (FIDUCEO)”. 
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

from fiduceo.fcdr.writer.fcdr_writer import FCDRWriter
from fiduceo.fcdr.writer.templates import avhrr
import sys
import netCDF4
import numpy as np
import datetime

class read_netcdf(object):

    def read_data(self,filename):

        ncid = netCDF4.Dataset(filename,'r')

        self.lat = ncid.variables['latitude'][:,:]
        self.lat = (self.lat*1000+0.5).astype(np.int32)/1000.
        self.lon = ncid.variables['longitude'][:,:]
        self.lon = (self.lon*1000+0.5).astype(np.int32)/1000.
        year = ncid.variables['year'][:]
        month = ncid.variables['month'][:]
        day = ncid.variables['day'][:]
        hours = ncid.variables['hours'][:]
        self.time = np.zeros(len(year))
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
        date_time = []
        for i in range(len(hour)):
            date_time = datetime.datetime(yr[i],mn[i],dy[i],hour[i],minute[i],second[i],microsec[i])
        self.time[gd] = netCDF4.date2num(date_time,'seconds since 1970-01-01')
        gd = (year < 1900) | (day < 0) | (month < 0) | (hours < 0)
        self.time[gd] = float('nan')
        self.satza = ncid.variables['satza'][:,:]
        self.solza = ncid.variables['solza'][:,:]
        self.relaz = ncid.variables['relaz'][:,:]
        self.ch1 = ncid.variables['ch1'][:,:]
        self.ch2 = ncid.variables['ch2'][:,:]
        try:
            self.ch3a = ncid.variables['ch3a'][:,:]
            self.ch3a_there = True
        except:
            self.ch3a_there = False
        self.ch3b = ncid.variables['ch3b'][:,:]
        self.ch4 = ncid.variables['ch4'][:,:]
        try:
            self.ch5 = ncid.variables['ch5'][:,:]
            self.ch5_there = True
        except:
            self.ch5_there = False
        self.u_random_ch1 = ncid.variables['ch1_random'][:,:]
        self.u_random_ch2 = ncid.variables['ch2_random'][:,:]
        if self.ch3a_there:
            self.u_random_ch3a = ncid.variables['ch3a_random'][:,:]
        self.u_random_ch3b = ncid.variables['ch3b_random'][:,:]
        self.u_random_ch4 = ncid.variables['ch4_random'][:,:]
        if self.ch5_there:
            self.u_random_ch5 = ncid.variables['ch5_random'][:,:]
        self.u_non_random_ch1 = ncid.variables['ch1_non_random'][:,:]
        self.u_non_random_ch2 = ncid.variables['ch2_non_random'][:,:]
        if self.ch3a_there:
            self.u_non_random_ch3a = ncid.variables['ch3a_non_random'][:,:]
        self.u_non_random_ch3b = ncid.variables['ch3b_non_random'][:,:]
        self.u_non_random_ch4 = ncid.variables['ch4_non_random'][:,:]
        if self.ch5_there:
            self.u_non_random_ch5 = ncid.variables['ch5_non_random'][:,:]

        self.nx = self.lat.shape[1]
        self.ny = self.lat.shape[0]

        ncid.close()

        self.lat = np.ma.filled(self.lat,float('nan'))
        self.lon = np.ma.filled(self.lon,float('nan'))
        self.time = np.ma.filled(self.time,float('nan'))
        self.satza = np.ma.filled(self.satza,float('nan'))
        self.solza = np.ma.filled(self.solza,float('nan'))
        self.relaz = np.ma.filled(self.relaz,float('nan'))
        self.ch1 = np.ma.filled(self.ch1,float('nan'))
        self.ch2 = np.ma.filled(self.ch2,float('nan'))
        if self.ch3a_there:
            self.ch3a = np.ma.filled(self.ch3a,float('nan'))
        self.ch3b = np.ma.filled(self.ch3b,float('nan'))
        self.ch4 = np.ma.filled(self.ch4,float('nan'))
        if self.ch5_there:
            self.ch5 = np.ma.filled(self.ch5,float('nan'))
        self.u_random_ch1 = np.ma.filled(self.u_random_ch1,float('nan'))
        self.u_random_ch2 = np.ma.filled(self.u_random_ch2,float('nan'))
        if self.ch3a_there:
            self.u_random_ch3a = np.ma.filled(self.u_random_ch3a,float('nan'))
        self.u_random_ch3b = np.ma.filled(self.u_random_ch3b,float('nan'))
        self.u_random_ch4 = np.ma.filled(self.u_random_ch4,float('nan'))
        if self.ch5_there:
            self.u_random_ch5 = np.ma.filled(self.u_random_ch5,float('nan'))
        self.u_non_random_ch1 = np.ma.filled(self.u_non_random_ch1,float('nan'))
        self.u_non_random_ch2 = np.ma.filled(self.u_non_random_ch2,float('nan'))
        if self.ch3a_there:
            self.u_non_random_ch3a = np.ma.filled(self.u_non_random_ch3a,float('nan'))
        self.u_non_random_ch3b = np.ma.filled(self.u_non_random_ch3b,float('nan'))
        self.u_non_random_ch4 = np.ma.filled(self.u_non_random_ch4,float('nan'))
        if self.ch5_there:
            self.u_non_random_ch5 = np.ma.filled(self.u_non_random_ch5,float('nan'))


    def __init__(self,filename):

        self.read_data(filename)

def main(file_in,file_out):
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
    print "File created"
    # write real data to the variables. All variables initially contain "_FillValue".
    # Not writing to the whole array is completely OK
    gd = np.isfinite(data.lat)
    dataset.variables["latitude"].data[gd] = data.lat[gd]
    gd = np.isfinite(data.lon)
    dataset.variables["longitude"].data[gd] = data.lon[gd]
    gd = np.isfinite(data.time)
    dataset.variables["Time"].data[gd] = data.time[gd]
    
    gd = np.isfinite(data.satza)
    dataset.variables["satellite_zenith_angle"].data[gd] = data.satza[gd]
    gd = np.isfinite(data.solza)
    dataset.variables["solar_zenith_angle"].data[gd] = data.solza[gd]
    gd = np.isfinite(data.relaz)
    dataset.variables["relative_azimuth_angle"].data[gd] = data.relaz[gd]

    gd = np.isfinite(data.ch1)
    dataset.variables["Ch1_Ref"].data[gd] = data.ch1[gd]
    gd = np.isfinite(data.ch2)
    dataset.variables["Ch2_Ref"].data[gd] = data.ch2[gd]
    if data.ch3a_there:
        gd = np.isfinite(data.ch3a)
        dataset.variables["Ch3a_Ref"].data[gd] = data.ch3a[gd]
    gd = np.isfinite(data.ch3b)
    dataset.variables["Ch3b_Bt"].data[gd] = data.ch3b[gd]
    gd = np.isfinite(data.ch4)
    dataset.variables["Ch4_Bt"].data[gd] = data.ch4[gd]
    if data.ch5_there:
        gd = np.isfinite(data.ch5)
        dataset.variables["Ch5_Bt"].data[gd] = data.ch5[gd]

    gd = np.isfinite(data.u_random_ch1)
    dataset.variables["u_random_Ch1"].data[gd] = data.u_random_ch1[gd]
    gd = np.isfinite(data.u_random_ch2)
    dataset.variables["u_random_Ch2"].data[gd] = data.u_random_ch2[gd]
    if data.ch3a_there:
        gd = np.isfinite(data.u_random_ch3a)
        dataset.variables["u_random_Ch3a"].data[gd] = data.u_random_ch3a[gd]
    gd = np.isfinite(data.u_random_ch3b)
    dataset.variables["u_random_Ch3b"].data[gd] = data.u_random_ch3b[gd]
    gd = np.isfinite(data.u_random_ch4)
    dataset.variables["u_random_Ch4"].data[gd] = data.u_random_ch4[gd]
    if data.ch5_there:
        gd = np.isfinite(data.u_random_ch5)
        dataset.variables["u_random_Ch5"].data[gd] = data.u_random_ch5[gd]
   
    gd = np.isfinite(data.u_non_random_ch1)
    dataset.variables["u_non_random_Ch1"].data[gd] = data.u_non_random_ch1[gd]
    gd = np.isfinite(data.u_non_random_ch2)
    dataset.variables["u_non_random_Ch2"].data[gd] = data.u_non_random_ch2[gd]
    if data.ch3a_there:
        gd = np.isfinite(data.u_non_random_ch3a)
        dataset.variables["u_non_random_Ch3a"].data[gd] = data.u_non_random_ch3a[gd]
    gd = np.isfinite(data.u_non_random_ch3b)
    dataset.variables["u_non_random_Ch3b"].data[gd] = data.u_non_random_ch3b[gd]
    gd = np.isfinite(data.u_non_random_ch4)
    dataset.variables["u_non_random_Ch4"].data[gd] = data.u_non_random_ch4[gd]
    if data.ch5_there:
        gd = np.isfinite(data.u_non_random_ch5)
        dataset.variables["u_non_random_Ch5"].data[gd] = data.u_non_random_ch5[gd]
   
    #avhrr.AVHRR._create_channel_refl_variable(12835, "Channel 6 Reflectance") 
    # dump it to disk, netcdf4, medium compression
    # writing will fail when the target file already exists
    writer.write(dataset, file_out)

if __name__ == "__main__":

    file_in=sys.argv[1]
    file_out=sys.argv[2]

    main(file_in,file_out)

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
        self.lon = ncid.variables['longitude'][:,:]
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

    dataset.variables["u_random_Ch1"].data = data.u_random_ch1
    dataset.variables["u_random_Ch2"].data = data.u_random_ch2
    if data.ch3a_there:
        dataset.variables["u_random_Ch3a"].data = data.u_random_ch3a
    dataset.variables["u_random_Ch3b"].data = data.u_random_ch3b
    dataset.variables["u_random_Ch4"].data = data.u_random_ch4
    if data.ch5_there:
        dataset.variables["u_random_Ch5"].data = data.u_random_ch5
   
    dataset.variables["u_non_random_Ch1"].data = data.u_non_random_ch1
    dataset.variables["u_non_random_Ch2"].data = data.u_non_random_ch2
    if data.ch3a_there:
        dataset.variables["u_non_random_Ch3a"].data = data.u_non_random_ch3a
    dataset.variables["u_non_random_Ch3b"].data = data.u_non_random_ch3b
    dataset.variables["u_non_random_Ch4"].data = data.u_non_random_ch4
    if data.ch5_there:
        dataset.variables["u_non_random_Ch5"].data = data.u_non_random_ch5
   
    #avhrr.AVHRR._create_channel_refl_variable(12835, "Channel 6 Reflectance") 
    # dump it to disk, netcdf4, medium compression
    # writing will fail when the target file already exists
    writer.write(dataset, file_out)

if __name__ == "__main__":

    file_in=sys.argv[1]
    file_out=sys.argv[2]

    main(file_in,file_out)

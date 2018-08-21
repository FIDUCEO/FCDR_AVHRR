# * Copyright (C) 2018 J.Mittaz University of Reading
# * This code was developed for the EC project "Fidelity and Uncertainty in
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

#
# Code to get general statistics from a single orbit of easy FCDR input
# data
#
# Written by J.Mittaz UoR
# --------------------------------

import numpy as np
import netCDF4 as nc
import datetime
from  optparse import OptionParser
import matplotlib.pyplot as plt

class read_stats_file(object):

    def __init__(self,filename):

        self.date=[]
        data = np.loadtxt(filename)
        chdata=np.zeros((1,6))
        for i in range(0,data.shape[0],18):
            year=data[i,0].astype(np.int32)
            month=data[i,1].astype(np.int32)
            day=data[i,2].astype(np.int32)
            hour=data[i,3].astype(np.int32)
            minute=data[i,4].astype(np.int32)
            second=data[i,5].astype(np.int32)
            self.date.append(datetime.datetime(year,month,day,hour,minute,second))
            if 0 == i:
                self.ch1_r=np.zeros((1,6))
                self.ch1_r[0,:] = data[i,6:]
                self.ch2_r=np.zeros((1,6))
                self.ch2_r[0,:] = data[i+1,6:]
                self.ch3a_r=np.zeros((1,6))
                self.ch3a_r[0,:] = data[i+2,6:]
                self.ch3b_r=np.zeros((1,6))
                self.ch3b_r[0,:] = data[i+3,6:]
                self.ch4_r=np.zeros((1,6))
                self.ch4_r[0,:] = data[i+4,6:]
                self.ch5_r=np.zeros((1,6))
                self.ch5_r[0,:] = data[i+5,6:]
                self.ch1_s=np.zeros((1,6))
                self.ch1_s[0,:] = data[i+6,6:]
                self.ch2_s=np.zeros((1,6))
                self.ch2_s[0,:] = data[i+7,6:]
                self.ch3a_s=np.zeros((1,6))
                self.ch3a_s[0,:] = data[i+8,6:]
                self.ch3b_s=np.zeros((1,6))
                self.ch3b_s[0,:] = data[i+9,6:]
                self.ch4_s=np.zeros((1,6))
                self.ch4_s[0,:] = data[i+10,6:]
                self.ch5_s=np.zeros((1,6))
                self.ch5_s[0,:] = data[i+11,6:]
                self.ch1_c=np.zeros((1,6))
                self.ch1_c[0,:] = data[i+12,6:]
                self.ch2_c=np.zeros((1,6))
                self.ch2_c[0,:] = data[i+13,6:]
                self.ch3a_c=np.zeros((1,6))
                self.ch3a_c[0,:] = data[i+14,6:]
                self.ch3b_c=np.zeros((1,6))
                self.ch3b_c[0,:] = data[i+15,6:]
                self.ch4_c=np.zeros((1,6))
                self.ch4_c[0,:] = data[i+16,6:]
                self.ch5_c=np.zeros((1,6))
                self.ch5_c[0,:] = data[i+17,6:]
            else:
                chdata[0,:]=data[i,6:]
                self.ch1_r = np.append(self.ch1_r,chdata,axis=0)
                chdata[0,:]=data[i+1,6:]
                self.ch2_r = np.append(self.ch2_r,chdata,axis=0)
                chdata[0,:]=data[i+2,6:]
                self.ch3a_r = np.append(self.ch3a_r,chdata,axis=0)
                chdata[0,:]=data[i+3,6:]
                self.ch3b_r = np.append(self.ch3b_r,chdata,axis=0)
                chdata[0,:]=data[i+4,6:]
                self.ch4_r = np.append(self.ch4_r,chdata,axis=0)
                chdata[0,:]=data[i+5,6:]
                self.ch5_r = np.append(self.ch5_r,chdata,axis=0)
                chdata[0,:]=data[i+6,6:]
                self.ch1_s = np.append(self.ch1_s,chdata,axis=0)
                chdata[0,:]=data[i+7,6:]
                self.ch2_s = np.append(self.ch2_s,chdata,axis=0)
                chdata[0,:]=data[i+8,6:]
                self.ch3a_s = np.append(self.ch3a_s,chdata,axis=0)
                chdata[0,:]=data[i+9,6:]
                self.ch3b_s = np.append(self.ch3b_s,chdata,axis=0)
                chdata[0,:]=data[i+10,6:]
                self.ch4_s = np.append(self.ch4_s,chdata,axis=0)
                chdata[0,:]=data[i+11,6:]
                self.ch5_s = np.append(self.ch5_s,chdata,axis=0)
                chdata[0,:]=data[i+12,6:]
                self.ch1_c = np.append(self.ch1_c,chdata,axis=0)
                chdata[0,:]=data[i+13,6:]
                self.ch2_c = np.append(self.ch2_c,chdata,axis=0)
                chdata[0,:]=data[i+14,6:]
                self.ch3a_c = np.append(self.ch3a_c,chdata,axis=0)
                chdata[0,:]=data[i+15,6:]
                self.ch3b_c = np.append(self.ch3b_c,chdata,axis=0)
                chdata[0,:]=data[i+16,6:]
                self.ch4_c = np.append(self.ch4_c,chdata,axis=0)
                chdata[0,:]=data[i+17,6:]
                self.ch5_c = np.append(self.ch5_c,chdata,axis=0)

def plot_comp(filename,filt=2):

    d=read_stats_file(filename)

    if 1 == filt:
        plt.figure(1)
        plt.plot(d.date,d.ch3b_r[:,2],',')
        plt.plot(d.date,d.ch3b_r[:,0],',')
        plt.plot(d.date,d.ch3b_r[:,1],',')
        plt.title('$3.7\mu$m Random')
        
        plt.figure(2)
        plt.plot(d.date,d.ch3b_s[:,2],',')
        plt.plot(d.date,d.ch3b_s[:,0],',')
        plt.plot(d.date,d.ch3b_s[:,1],',')
        plt.title('$3.7\mu$m Systematic')
        
        plt.figure(3)
        plt.plot(d.date,d.ch3b_c[:,2],',')
        plt.plot(d.date,d.ch3b_c[:,0],',')
        plt.plot(d.date,d.ch3b_c[:,1],',')
        plt.title('$3.7\mu$m Common')
    elif 2 == filt:
        plt.figure(1)
        plt.plot(d.date,d.ch4_r[:,2],',')
        plt.plot(d.date,d.ch4_r[:,0],',')
        plt.plot(d.date,d.ch4_r[:,1],',')
        plt.title('$11\mu$m Random')
        
        plt.figure(2)
        plt.plot(d.date,d.ch4_s[:,2],',')
        plt.plot(d.date,d.ch4_s[:,0],',')
        plt.plot(d.date,d.ch4_s[:,1],',')
        plt.title('$11\mu$m Systematic')
        
        plt.figure(3)
        plt.plot(d.date,d.ch4_c[:,2],',')
        plt.plot(d.date,d.ch4_c[:,0],',')
        plt.plot(d.date,d.ch4_c[:,1],',')
        plt.title('$11\mu$m Common')
    elif 3 == filt:
        plt.figure(1)
        plt.plot(d.date,d.ch5_r[:,2],',')
        plt.plot(d.date,d.ch5_r[:,0],',')
        plt.plot(d.date,d.ch5_r[:,1],',')
        plt.title('$12\mu$m Random')
        
        plt.figure(2)
        plt.plot(d.date,d.ch5_s[:,2],',')
        plt.plot(d.date,d.ch5_s[:,0],',')
        plt.plot(d.date,d.ch5_s[:,1],',')
        plt.title('$12\mu$m Systematic')
        
        plt.figure(3)
        plt.plot(d.date,d.ch5_c[:,2],',')
        plt.plot(d.date,d.ch5_c[:,0],',')
        plt.plot(d.date,d.ch5_c[:,1],',')
        plt.title('$12\mu$m Common')

    plt.show()

class read_temp_file(object):

    def read_data(self,filename):

        ncid = nc.Dataset(filename,'r')
	
        self.lat = ncid.variables['latitude'][:,:]
        self.lon = ncid.variables['longitude'][:,:]
        self.time = ncid.variables['time'][:]
        self.date=[]
        gd = (self.time > -1e20)
        for i in range(len(self.time)):
            if gd[i]:
                self.date.append(nc.num2date(self.time[i],\
                                                 'seconds since 1975-01-01'))
            else:
                self.date.append(None)
        self.satza = ncid.variables['satza'][:,:]
        self.solza = ncid.variables['solza'][:,:]
        self.relaz = ncid.variables['relaz'][:,:]
        self.ch1 = ncid.variables['ch1'][:,:]
        self.ch2 = ncid.variables['ch2'][:,:]
        self.ch3a = ncid.variables['ch3a'][:,:]
        self.ch3b = ncid.variables['ch3b'][:,:]
        self.ch4 = ncid.variables['ch4'][:,:]
        self.ch5 = ncid.variables['ch5'][:,:]
        self.ch1_random = ncid.variables['ch1_random'][:,:]
        self.ch2_random = ncid.variables['ch2_random'][:,:]
        self.ch3a_random = ncid.variables['ch3a_random'][:,:]
        self.ch3b_random = ncid.variables['ch3b_random'][:,:]
        self.ch4_random = ncid.variables['ch4_random'][:,:]
        self.ch5_random = ncid.variables['ch5_random'][:,:]
        self.ch1_non_random = ncid.variables['ch1_non_random'][:,:]
        self.ch2_non_random = ncid.variables['ch2_non_random'][:,:]
        self.ch3a_non_random = ncid.variables['ch3a_non_random'][:,:]
        self.ch3b_non_random = ncid.variables['ch3b_non_random'][:,:]
        self.ch4_non_random = ncid.variables['ch4_non_random'][:,:]
        self.ch5_non_random = ncid.variables['ch5_non_random'][:,:]
        self.ch1_common = ncid.variables['ch1_common'][:,:]
        self.ch2_common = ncid.variables['ch2_common'][:,:]
        self.ch3a_common = ncid.variables['ch3a_common'][:,:]
        self.ch3b_common = ncid.variables['ch3b_common'][:,:]
        self.ch4_common = ncid.variables['ch4_common'][:,:]
        self.ch5_common = ncid.variables['ch5_common'][:,:]

        ncid.close()

    def __init__(self,filename):

        self.read_data(filename)

def get_values(data):

    gd = np.isfinite(data) & (data > -1e20)
    if 0 < np.sum(gd):
        minval = data[gd].min()
        maxval = data[gd].max()
        mean = np.mean(data[gd])
        q25 = np.percentile(data[gd],25)
        q75 = np.percentile(data[gd],75)
        median = np.median(data[gd])
    else:
        minval = float('nan')
        maxval = float('nan')
        mean = float('nan')
        q25 = float('nan')
        q75 = float('nan')
        median = float('nan')

    return minval,maxval,mean,q25,q75,median

def write_stats(fp,date,indata):

    minval,maxval,mean,q25,q75,median = get_values(indata)
    fp.write('{0:d} {1:d} {2:d} {3:d} {4:d} {5:d} {6:e} {7:e} {8:e} {9:e} {10:e} {11:e}\n'.\
                 format(date.year,date.month,date.day,date.hour,\
                            date.minute,date.second,minval,maxval,mean,\
                            q25,q75,median))

def get_stats(uuid,ofile):

    infile = uuid+'.nc'
    data = read_temp_file(infile)

    with open(ofile,'w') as fp:

        # Get simple stats
        write_stats(fp,data.date[0],data.ch1_random)
        write_stats(fp,data.date[0],data.ch2_random)
        write_stats(fp,data.date[0],data.ch3a_random)
        write_stats(fp,data.date[0],data.ch3b_random)
        write_stats(fp,data.date[0],data.ch4_random)
        write_stats(fp,data.date[0],data.ch5_random)

        write_stats(fp,data.date[0],data.ch1_non_random)
        write_stats(fp,data.date[0],data.ch2_non_random)
        write_stats(fp,data.date[0],data.ch3a_non_random)
        write_stats(fp,data.date[0],data.ch3b_non_random)
        write_stats(fp,data.date[0],data.ch4_non_random)
        write_stats(fp,data.date[0],data.ch5_non_random)

        write_stats(fp,data.date[0],data.ch1_common)
        write_stats(fp,data.date[0],data.ch2_common)
        write_stats(fp,data.date[0],data.ch3a_common)
        write_stats(fp,data.date[0],data.ch3b_common)
        write_stats(fp,data.date[0],data.ch4_common)
        write_stats(fp,data.date[0],data.ch5_common)

if "__main__" == __name__:
    
    parser = OptionParser("usage: %prog uuid outfile")
    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.error("incorrect number of arguments")

    get_stats(args[0],args[1])

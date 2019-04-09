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
import dateutil
import os
import re
import matplotlib.dates as mdates

class read_stats_file(object):

    def read_netcdf(self,filename):

        ncid = nc.Dataset(filename,'r')

        ch1_rand = ncid.variables['ch1_random'][:,:]
        ch2_rand = ncid.variables['ch2_random'][:,:]
        ch3a_rand = ncid.variables['ch3a_random'][:,:]
        ch3b_rand = ncid.variables['ch3b_random'][:,:]
        ch4_rand = ncid.variables['ch4_random'][:,:]
        ch5_rand = ncid.variables['ch5_random'][:,:]
        ch1_non_rand = ncid.variables['ch1_non_random'][:,:]
        ch2_non_rand = ncid.variables['ch2_non_random'][:,:]
        ch3a_non_rand = ncid.variables['ch3a_non_random'][:,:]
        ch3b_non_rand = ncid.variables['ch3b_non_random'][:,:]
        ch4_non_rand = ncid.variables['ch4_non_random'][:,:]
        ch5_non_rand = ncid.variables['ch5_non_random'][:,:]
        ch1_common = ncid.variables['ch1_common'][:,:]
        ch2_common = ncid.variables['ch2_common'][:,:]
        ch3a_common = ncid.variables['ch3a_common'][:,:]
        ch3b_common = ncid.variables['ch3b_common'][:,:]
        ch4_common = ncid.variables['ch4_common'][:,:]
        ch5_common = ncid.variables['ch5_common'][:,:]
        data_map = ncid.variables['map'][:,:,:]
        ndata_map = ncid.variables['nmap'][:,:,:]
        datestr = ncid.variables['date'][:]

        date = []
        for i in range(len(datestr)):
            date.append(dateutil.parser.parse(datestr[i]))

        ncid.close()

        data = np.zeros((len(datestr),14))

        return date,data,data_map,ndata_map

    def __init__(self,filename,ascii=True):

        self.date=[]
        if ascii:
            data = np.loadtxt(filename)
        else:
            data,data_map,ndata_map = self.read_netcdf(filename)
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
                if not ascii:
                    self.total=np.zeros((1,3,18))
                    self.total[0,:,:] = data_map
                    self.ntotal=np.zeros((1,3,18),dtype=np.int32)
                    self.ntotal[0,:,:] = ndata_map
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
                if not ascii:
                    self.total = np.append(self.total,data_map,axis=0)
                    self.ntotal = np.append(self.ntotal,ndata_map,axis=0)
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

def get_max_date(date,array):

    gd = np.isfinite(array)
    newarray = array[gd]
    newdate=[]
    for i in range(len(gd)):
        if gd[i]:
            newdate.append(date[i])
    gd = (newarray == newarray.max())
    for i in range(len(gd)):
        if gd[i]:
            maxdate = newdate[i]
    return maxdate.ctime()

def plot_mean(filename,filt=2,ascii=True,utype=0):

    d=read_stats_file(filename,ascii=ascii)
    
    if filename == 'avhrr11.stats':
        name = 'NOAA11'
    elif filename == 'avhrr12.stats':
        name = 'NOAA12'
    elif filename == 'avhrr14.stats':
        name = 'NOAA14'
    elif filename == 'avhrr15.stats':
        name = 'NOAA15'
    elif filename == 'avhrr16.stats':
        name = 'NOAA16'
    elif filename == 'avhrr17.stats':
        name = 'NOAA17'
    elif filename == 'avhrr18.stats':
        name = 'NOAA18'
    elif filename == 'avhrr19.stats':
        name = 'NOAA19'
    elif filename == 'avhrrmta.stats':
        name = 'MetOp-A'


    fig,ax=plt.subplots()

    if 1 == filt:
        plt.plot(d.date,d.ch3b_r[:,utype],'k,')
        plt.plot(d.date,d.ch3b_s[:,utype],'b,')
        plt.plot(d.date,d.ch3b_c[:,utype],'g,')
        name = name + ' 3.7$\mu$m'
    elif 2 == filt:
        plt.plot(d.date,d.ch4_r[:,utype],'k,')
        plt.plot(d.date,d.ch4_s[:,utype],'b,')
        plt.plot(d.date,d.ch4_c[:,utype],'g,')
        name = name + ' 11$\mu$m'
    elif 3 == filt:
        plt.plot(d.date,d.ch5_r[:,utype],'k,')
        plt.plot(d.date,d.ch5_s[:,utype],'b,')
        plt.plot(d.date,d.ch5_c[:,utype],'g,')
        name = name + ' 12$\mu$m'

    plt.title(name)
    if 0 == utype:
        plt.ylabel('Mean Uncert. / K')
    elif 1 == utype:
        plt.ylabel('Min Uncert. / K')
    elif 2 == utype:
        plt.ylabel('Max Uncert. / K')
    plt.xlabel('Date')

    plt.text(0.98, 0.95, 'Independent', horizontalalignment='right',\
                 verticalalignment='center', transform=ax.transAxes,\
                 color='black')
    plt.text(0.98, 0.9, 'Structured', horizontalalignment='right',\
                 verticalalignment='center', transform=ax.transAxes,\
                 color='blue')
    plt.text(0.98, 0.85, 'Common', horizontalalignment='right',\
                 verticalalignment='center', transform=ax.transAxes,\
                 color='green')
    plt.text(0.98, 0.85, 'Oper. Ne$\Delta$T', horizontalalignment='right',\
                 verticalalignment='center', transform=ax.transAxes,\
                 color='ref')
    plt.plot([d.date[0],d.date[-1]],[0.12,0.12],'r--')

    plt.show()

def plot_comp(filename,filt=2,ascii=True):

    d=read_stats_file(filename,ascii=ascii)

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
        print '$3.7\mu$m Systematic (max value date: {0}'.\
                      format(get_max_date(d.date,d.ch3b_s[:,1]))
        
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
        print '$11\mu$m Systematic (max value date: {0}'.\
                      format(get_max_date(d.date,d.ch4_s[:,1]))
        
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
        print '$12\mu$m Systematic (max value date: {0}'.\
                      format(get_max_date(d.date,d.ch5_s[:,1]))
        
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

def get_values(indata):

    if np.ma.is_masked(indata):
        data=indata.filled(-1e30)
    else:
        data=indata
    gd = np.isfinite(data) & (data > -1e20)
    if 0 < np.sum(gd):
        minval = data[gd].min()
        maxval = data[gd].max()
        mean = np.mean(data[gd])
        q25 = np.percentile(data[gd],25)
        q75 = np.percentile(data[gd],75)
        median = np.median(data[gd])
    else:
        minval = -1e30
        maxval = -1e30
        mean = -1e30
        q25 = -1e30
        q75 = -1e30
        median = -1e30

    return minval,maxval,mean,q25,q75,median

def write_stats(fp,date,indata):

    minval,maxval,mean,q25,q75,median = get_values(indata)
    fp.write('{0:d} {1:d} {2:d} {3:d} {4:d} {5:d} {6:e} {7:e} {8:e} {9:e} {10:e} {11:e}\n'.\
                 format(date.year,date.month,date.day,date.hour,\
                            date.minute,date.second,minval,maxval,mean,\
                            q25,q75,median))

def get_map_data(data):

    data_map = np.zeros((6,18),dtype=np.float32)
    ndata_map  = np.zeros((6,18),dtype=np.int32)

    for i in range(18):
        minlat = i*10.-90.
        maxlat = (i+1)*10.-90.
        gd = (data.lat >= minlat) & (data.lat < maxlat) & (data.ch1 >= 0.) & (data.solza >=0) & (data.solza < 90.)
        if np.sum(gd) > 0:
            data_map[0,i] = np.sum(data.ch1[gd])
            ndata_map[0,i] = np.sum(gd)
        gd = (data.lat >= minlat) & (data.lat < maxlat) & (data.ch2 >= 0.) & (data.solza >=0) & (data.solza < 90.)
        if np.sum(gd) > 0:
            data_map[1,i] = np.sum(data.ch2[gd])
            ndata_map[1,i] = np.sum(gd)
        gd = (data.lat >= minlat) & (data.lat < maxlat) & (data.ch3a >= 0.) & (data.solza >=0) & (data.solza < 90.)
        if np.sum(gd) > 0:
            data_map[2,i] = np.sum(data.ch3a[gd])
            ndata_map[2,i] = np.sum(gd)
        gd = (data.lat >= minlat) & (data.lat < maxlat) & (data.ch3b >= 0.) & (data.solza > 90.)
        if np.sum(gd) > 0:
            data_map[3,i] = np.sum(data.ch3b[gd])
            ndata_map[3,i] = np.sum(gd)
        gd = (data.lat >= minlat) & (data.lat < maxlat) & (data.ch4 >= 0.) & (data.solza > 90.)
        if np.sum(gd) > 0:
            data_map[4,i] = np.sum(data.ch4[gd])
            ndata_map[4,i] = np.sum(gd)
        gd = (data.lat >= minlat) & (data.lat < maxlat) & (data.ch5 >= 0.) & (data.solza > 90.)
        if np.sum(gd) > 0:
            data_map[5,i] = np.sum(data.ch5[gd])
            ndata_map[5,i] = np.sum(gd)

    return data_map,ndata_map

def write_stats_netcdf(ofile,data):

    outfile = ofile+'.nc'

    for i in range(len(data.date)):
        if data.date[i].year > 0 and data.date[i].year < 2100:
            date = data.date[i]
            date_str = date.isoformat()
            break

    minval,maxval,mean,q35,q75,median = get_values(data.ch1_random)
    ch1_rand = np.array([minval,maxval,mean,q35,q75,median])
    minval,maxval,mean,q35,q75,median = get_values(data.ch2_random)
    ch2_rand = np.array([minval,maxval,mean,q35,q75,median])
    minval,maxval,mean,q35,q75,median = get_values(data.ch3a_random)
    ch3a_rand = np.array([minval,maxval,mean,q35,q75,median])
    minval,maxval,mean,q35,q75,median = get_values(data.ch3b_random)
    ch3b_rand = np.array([minval,maxval,mean,q35,q75,median])
    minval,maxval,mean,q35,q75,median = get_values(data.ch4_random)
    ch4_rand = np.array([minval,maxval,mean,q35,q75,median])
    minval,maxval,mean,q35,q75,median = get_values(data.ch5_random)
    ch5_rand = np.array([minval,maxval,mean,q35,q75,median])
    
    minval,maxval,mean,q35,q75,median = get_values(data.ch1_non_random)
    ch1_non_rand = np.array([minval,maxval,mean,q35,q75,median])
    minval,maxval,mean,q35,q75,median = get_values(data.ch2_non_random)
    ch2_non_rand = np.array([minval,maxval,mean,q35,q75,median])
    minval,maxval,mean,q35,q75,median = get_values(data.ch3a_non_random)
    ch3a_non_rand = np.array([minval,maxval,mean,q35,q75,median])
    minval,maxval,mean,q35,q75,median = get_values(data.ch3b_non_random)
    ch3b_non_rand = np.array([minval,maxval,mean,q35,q75,median])
    minval,maxval,mean,q35,q75,median = get_values(data.ch4_non_random)
    ch4_non_rand = np.array([minval,maxval,mean,q35,q75,median])
    minval,maxval,mean,q35,q75,median = get_values(data.ch5_non_random)
    ch5_non_rand = np.array([minval,maxval,mean,q35,q75,median])
    
    minval,maxval,mean,q35,q75,median = get_values(data.ch1_common)
    ch1_common = np.array([minval,maxval,mean,q35,q75,median])
    minval,maxval,mean,q35,q75,median = get_values(data.ch2_common)
    ch2_common = np.array([minval,maxval,mean,q35,q75,median])
    minval,maxval,mean,q35,q75,median = get_values(data.ch3a_common)
    ch3a_common = np.array([minval,maxval,mean,q35,q75,median])
    minval,maxval,mean,q35,q75,median = get_values(data.ch3b_common)
    ch3b_common = np.array([minval,maxval,mean,q35,q75,median])
    minval,maxval,mean,q35,q75,median = get_values(data.ch4_common)
    ch4_common = np.array([minval,maxval,mean,q35,q75,median])
    minval,maxval,mean,q35,q75,median = get_values(data.ch5_common)
    ch5_common = np.array([minval,maxval,mean,q35,q75,median])
    
    data_map,ndata_map = get_map_data(data)

    ncid = nc.Dataset(outfile,'w')

    nelem = ncid.createDimension('nelem',size=6)
    nlat = ncid.createDimension('nlat',size=18)
    nchan = ncid.createDimension('nchan',size=6)

    ch1rand = ncid.createVariable('ch1_random','f4',dimensions=('nelem'),\
                                      zlib=True,fill_value=-1e30)
    ch2rand = ncid.createVariable('ch2_random','f4',dimensions=('nelem'),\
                                      zlib=True,fill_value=-1e30)
    ch3arand = ncid.createVariable('ch3a_random','f4',dimensions=('nelem'),\
                                      zlib=True,fill_value=-1e30)
    ch3brand = ncid.createVariable('ch3b_random','f4',dimensions=('nelem'),\
                                      zlib=True,fill_value=-1e30)
    ch4rand = ncid.createVariable('ch4_random','f4',dimensions=('nelem'),\
                                      zlib=True,fill_value=-1e30)
    ch5rand = ncid.createVariable('ch5_random','f4',dimensions=('nelem'),\
                                      zlib=True,fill_value=-1e30)

    ch1nonrand = ncid.createVariable('ch1_non_random','f4',\
                                         dimensions=('nelem'),\
                                         zlib=True,fill_value=-1e30)
    ch2nonrand = ncid.createVariable('ch2_non_random','f4',\
                                         dimensions=('nelem'),\
                                         zlib=True,fill_value=-1e30)
    ch3anonrand = ncid.createVariable('ch3a_non_random','f4',\
                                          dimensions=('nelem'),\
                                          zlib=True,fill_value=-1e30)
    ch3bnonrand = ncid.createVariable('ch3b_non_random','f4',\
                                          dimensions=('nelem'),\
                                          zlib=True,fill_value=-1e30)
    ch4nonrand = ncid.createVariable('ch4_non_random','f4',\
                                         dimensions=('nelem'),\
                                         zlib=True,fill_value=-1e30)
    ch5nonrand = ncid.createVariable('ch5_non_random','f4',\
                                         dimensions=('nelem'),\
                                         zlib=True,fill_value=-1e30)

    ch1common = ncid.createVariable('ch1_common','f4',\
                                        dimensions=('nelem'),\
                                        zlib=True,fill_value=-1e30)
    ch2common = ncid.createVariable('ch2_common','f4',\
                                        dimensions=('nelem'),\
                                        zlib=True,fill_value=-1e30)
    ch3acommon = ncid.createVariable('ch3a_common','f4',\
                                         dimensions=('nelem'),\
                                         zlib=True,fill_value=-1e30)
    ch3bcommon = ncid.createVariable('ch3b_common','f4',\
                                        dimensions=('nelem'),\
                                        zlib=True,fill_value=-1e30)
    ch4common = ncid.createVariable('ch4_common','f4',\
                                        dimensions=('nelem'),\
                                        zlib=True,fill_value=-1e30)
    ch5common = ncid.createVariable('ch5_common','f4',\
                                        dimensions=('nelem'),\
                                        zlib=True,fill_value=-1e30)

    out_datamap = ncid.createVariable('map','f4',\
                                          dimensions=('nchan','nlat'),\
                                          zlib=True,fill_value=-1e30)

    out_ndatamap = ncid.createVariable('nmap','i4',\
                                           dimensions=('nchan','nlat'),\
                                           zlib=True,fill_value=-999)

    ch1rand[:] = ch1_rand
    ch2rand[:] = ch2_rand
    ch3arand[:] = ch3a_rand
    ch3brand[:] = ch3b_rand
    ch4rand[:] = ch4_rand
    ch5rand[:] = ch5_rand

    ch1nonrand[:] = ch1_non_rand
    ch2nonrand[:] = ch2_non_rand
    ch3anonrand[:] = ch3a_non_rand
    ch3bnonrand[:] = ch3b_non_rand
    ch4nonrand[:] = ch4_non_rand
    ch5nonrand[:] = ch5_non_rand

    ch1common[:] = ch1_common
    ch2common[:] = ch2_common
    ch3acommon[:] = ch3a_common
    ch3bcommon[:] = ch3b_common
    ch4common[:] = ch4_common
    ch5common[:] = ch5_common

    out_datamap[:,:] = data_map
    out_ndatamap[:,:] = ndata_map

    ncid.date = date_str

    ncid.close()

def get_stats(uuid,ofile,netcdf_out):

    infile = uuid+'.nc'
    data = read_temp_file(infile)

    if netcdf_out:
        write_stats_netcdf(ofile,data)
    else:
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

class merge_stats(object):

    def add(self,filename):

        ncid = nc.Dataset(filename,'r')

        try:
            date = dateutil.parser.parse(ncid.date)
        
            ch1_rand = ncid.variables['ch1_random'][:]
            ch2_rand = ncid.variables['ch2_random'][:]
            ch3a_rand = ncid.variables['ch3a_random'][:]
            ch3b_rand = ncid.variables['ch3b_random'][:]
            ch4_rand = ncid.variables['ch4_random'][:]
            ch5_rand = ncid.variables['ch5_random'][:]
            
            ch1_sys = ncid.variables['ch1_non_random'][:]
            ch2_sys = ncid.variables['ch2_non_random'][:]
            ch3a_sys = ncid.variables['ch3a_non_random'][:]
            ch3b_sys = ncid.variables['ch3b_non_random'][:]
            ch4_sys = ncid.variables['ch4_non_random'][:]
            ch5_sys = ncid.variables['ch5_non_random'][:]
            
            ch1_com = ncid.variables['ch1_common'][:]
            ch2_com = ncid.variables['ch2_common'][:]
            ch3a_com = ncid.variables['ch3a_common'][:]
            ch3b_com = ncid.variables['ch3b_common'][:]
            ch4_com = ncid.variables['ch4_common'][:]
            ch5_com = ncid.variables['ch5_common'][:]
            
            map_vals = ncid.variables['map'][:,:]
            nmap_vals = ncid.variables['nmap'][:,:]
            
            ncid.close()
            
            if not self.init:
                
                self.ndata = 0
                self.ch1_rand = np.zeros((1,len(ch1_rand)))
                self.ch2_rand = np.zeros((1,len(ch2_rand)))
                self.ch3a_rand = np.zeros((1,len(ch3a_rand)))
                self.ch3b_rand = np.zeros((1,len(ch3b_rand)))
                self.ch4_rand = np.zeros((1,len(ch4_rand)))
                self.ch5_rand = np.zeros((1,len(ch5_rand)))
                
                self.ch1_sys = np.zeros((1,len(ch1_rand)))
                self.ch2_sys = np.zeros((1,len(ch2_rand)))
                self.ch3a_sys = np.zeros((1,len(ch3a_rand)))
                self.ch3b_sys = np.zeros((1,len(ch3b_rand)))
                self.ch4_sys = np.zeros((1,len(ch4_rand)))
                self.ch5_sys = np.zeros((1,len(ch5_rand)))
                
                self.ch1_com = np.zeros((1,len(ch1_rand)))
                self.ch2_com = np.zeros((1,len(ch2_rand)))
                self.ch3a_com = np.zeros((1,len(ch3a_rand)))
                self.ch3b_com = np.zeros((1,len(ch3b_rand)))
                self.ch4_com = np.zeros((1,len(ch4_rand)))
                self.ch5_com = np.zeros((1,len(ch5_rand)))
                
                self.map_vals = np.zeros((1,map_vals.shape[0],map_vals.shape[1]),\
                                             dtype=np.float32)
                self.nmap_vals = np.zeros((1,map_vals.shape[0],map_vals.shape[1]),\
                                              dtype=np.int32)
                
                self.datetime = []
                
                self.ch1_rand[0,:] = ch1_rand
                self.ch2_rand[0,:] = ch2_rand
                self.ch3a_rand[0,:] = ch3a_rand
                self.ch3b_rand[0,:] = ch3b_rand
                self.ch4_rand[0,:] = ch4_rand
                self.ch5_rand[0,:] = ch5_rand
                
                self.ch1_sys[0,:] = ch1_sys
                self.ch2_sys[0,:] = ch2_sys
                self.ch3a_sys[0,:] = ch3a_sys
                self.ch3b_sys[0,:] = ch3b_sys
                self.ch4_sys[0,:] = ch4_sys
                self.ch5_sys[0,:] = ch5_sys
                
                self.ch1_com[0,:] = ch1_com
                self.ch2_com[0,:] = ch2_com
                self.ch3a_com[0,:] = ch3a_com
                self.ch3b_com[0,:] = ch3b_com
                self.ch4_com[0,:] = ch4_com
                self.ch5_com[0,:] = ch5_com
                
                self.map_vals[0,:,:] = map_vals
                self.nmap_vals[0,:,:] = nmap_vals
                
                self.datetime.append(date)
                
            else:
                
                self.ch1_rand = np.append(self.ch1_rand,\
                                              ch1_rand[np.newaxis,...],axis=0)
                self.ch2_rand = np.append(self.ch2_rand,\
                                              ch2_rand[np.newaxis,...],axis=0)
                self.ch3a_rand = np.append(self.ch3a_rand,\
                                               ch3a_rand[np.newaxis,...],axis=0)
                self.ch3b_rand = np.append(self.ch3b_rand,\
                                               ch3b_rand[np.newaxis,...],axis=0)
                self.ch4_rand = np.append(self.ch4_rand,\
                                              ch4_rand[np.newaxis,...],axis=0)
                self.ch5_rand = np.append(self.ch5_rand,\
                                              ch5_rand[np.newaxis,...],axis=0)
                
                self.ch1_sys = np.append(self.ch1_sys,\
                                             ch1_sys[np.newaxis,...],axis=0)
                self.ch2_sys = np.append(self.ch2_sys,\
                                             ch2_sys[np.newaxis,...],axis=0)
                self.ch3a_sys = np.append(self.ch3a_sys,\
                                              ch3a_sys[np.newaxis,...],axis=0)
                self.ch3b_sys = np.append(self.ch3b_sys,\
                                              ch3b_sys[np.newaxis,...],axis=0)
                self.ch4_sys = np.append(self.ch4_sys,\
                                             ch4_sys[np.newaxis,...],axis=0)
                self.ch5_sys = np.append(self.ch5_sys,\
                                             ch5_sys[np.newaxis,...],axis=0)
                
                self.ch1_com = np.append(self.ch1_com,\
                                             ch1_com[np.newaxis,...],axis=0)
                self.ch2_com = np.append(self.ch2_com,\
                                             ch2_com[np.newaxis,...],axis=0)
                self.ch3a_com = np.append(self.ch3a_com,\
                                              ch3a_com[np.newaxis,...],axis=0)
                self.ch3b_com = np.append(self.ch3b_com,\
                                              ch3b_com[np.newaxis,...],axis=0)
                self.ch4_com = np.append(self.ch4_com,\
                                             ch4_com[np.newaxis,...],axis=0)
                self.ch5_com = np.append(self.ch5_com,\
                                             ch5_com[np.newaxis,...],axis=0)
                
                self.map_vals = np.append(self.map_vals,\
                                              map_vals[np.newaxis,...],axis=0)
                self.nmap_vals = np.append(self.nmap_vals,\
                                               nmap_vals[np.newaxis,...],axis=0)
            
                self.datetime.append(date)
            
            self.init = True
            self.ndata = self.ndata + 1                                    
        
        except:

            ncid.close()

    def write(self,ofile):

        ncid = nc.Dataset(ofile,'w')

        ndata = ncid.createDimension('ndata',size=self.ndata)
        nelem = ncid.createDimension('nelem',size=self.ch1_rand.shape[1])
        nlat = ncid.createDimension('nlat',size=18)
        nchan = ncid.createDimension('nchan',size=6)

        time = ncid.createVariable('time','f8',\
                                       dimensions=('ndata'),\
                                       zlib=True,fill_value=-1e30)
        time.units='seconds since 1970-01-01'

        ch1rand = ncid.createVariable('ch1_random','f4',\
                                          dimensions=('ndata','nelem'),\
                                          zlib=True,fill_value=-1e30)
        ch2rand = ncid.createVariable('ch2_random','f4',\
                                          dimensions=('ndata','nelem'),\
                                          zlib=True,fill_value=-1e30)
        ch3arand = ncid.createVariable('ch3a_random','f4',\
                                          dimensions=('ndata','nelem'),\
                                          zlib=True,fill_value=-1e30)
        ch3brand = ncid.createVariable('ch3b_random','f4',\
                                          dimensions=('ndata','nelem'),\
                                          zlib=True,fill_value=-1e30)
        ch4rand = ncid.createVariable('ch4_random','f4',\
                                          dimensions=('ndata','nelem'),\
                                          zlib=True,fill_value=-1e30)
        ch5rand = ncid.createVariable('ch5_random','f4',\
                                          dimensions=('ndata','nelem'),\
                                          zlib=True,fill_value=-1e30)

        ch1sys = ncid.createVariable('ch1_structured','f4',\
                                         dimensions=('ndata','nelem'),\
                                         zlib=True,fill_value=-1e30)
        ch2sys = ncid.createVariable('ch2_structured','f4',\
                                         dimensions=('ndata','nelem'),\
                                         zlib=True,fill_value=-1e30)
        ch3asys = ncid.createVariable('ch3a_structured','f4',\
                                          dimensions=('ndata','nelem'),\
                                          zlib=True,fill_value=-1e30)
        ch3bsys = ncid.createVariable('ch3b_structured','f4',\
                                          dimensions=('ndata','nelem'),\
                                          zlib=True,fill_value=-1e30)
        ch4sys = ncid.createVariable('ch4_structured','f4',\
                                         dimensions=('ndata','nelem'),\
                                         zlib=True,fill_value=-1e30)
        ch5sys = ncid.createVariable('ch5_structured','f4',\
                                         dimensions=('ndata','nelem'),\
                                         zlib=True,fill_value=-1e30)

        ch1com = ncid.createVariable('ch1_common','f4',\
                                         dimensions=('ndata','nelem'),\
                                         zlib=True,fill_value=-1e30)
        ch2com = ncid.createVariable('ch2_common','f4',\
                                         dimensions=('ndata','nelem'),\
                                         zlib=True,fill_value=-1e30)
        ch3acom = ncid.createVariable('ch3a_common','f4',\
                                          dimensions=('ndata','nelem'),\
                                          zlib=True,fill_value=-1e30)
        ch3bcom = ncid.createVariable('ch3b_common','f4',\
                                          dimensions=('ndata','nelem'),\
                                          zlib=True,fill_value=-1e30)
        ch4com = ncid.createVariable('ch4_common','f4',\
                                         dimensions=('ndata','nelem'),\
                                         zlib=True,fill_value=-1e30)
        ch5com = ncid.createVariable('ch5_common','f4',\
                                         dimensions=('ndata','nelem'),\
                                         zlib=True,fill_value=-1e30)

        mapdata = ncid.createVariable('map_data','f4',\
                                          dimensions=('ndata','nchan','nlat'),\
                                          zlib=True,fill_value=-1e30)

        nmapdata = ncid.createVariable('nmap_data','i4',\
                                           dimensions=('ndata','nchan','nlat'),\
                                           zlib=True,fill_value=-999)

        timeval = nc.date2num(self.datetime,'seconds since 1970-01-01')

        index = np.argsort(timeval)

        time[:] = timeval[index]
        ch1rand[:,:] = self.ch1_rand[index,:]
        ch2rand[:,:] = self.ch2_rand[index,:]
        ch3arand[:,:] = self.ch3a_rand[index,:]
        ch3brand[:,:] = self.ch3b_rand[index,:]
        ch4rand[:,:] = self.ch4_rand[index,:]
        ch5rand[:,:] = self.ch5_rand[index,:]

        ch1sys[:,:] = self.ch1_sys[index,:]
        ch2sys[:,:] = self.ch2_sys[index,:]
        ch3asys[:,:] = self.ch3a_sys[index,:]
        ch3bsys[:,:] = self.ch3b_sys[index,:]
        ch4sys[:,:] = self.ch4_sys[index,:]
        ch5sys[:,:] = self.ch5_sys[index,:]

        ch1com[:,:] = self.ch1_com[index,:]
        ch2com[:,:] = self.ch2_com[index,:]
        ch3acom[:,:] = self.ch3a_com[index,:]
        ch3bcom[:,:] = self.ch3b_com[index,:]
        ch4com[:,:] = self.ch4_com[index,:]
        ch5com[:,:] = self.ch5_com[index,:]

        mapdata[:,:,:] = self.map_vals[index,:,:]
        nmapdata[:,:,:] = self.nmap_vals[index,:,:]

        ncid.close()

    def __init__(self,directory='.'):

        prog=re.compile('.*/stats.*dat.*nc')
        self.init=False
        for dirpath,dirs,files in os.walk(directory):
            for filename in files:
                fname = os.path.join(dirpath,filename)
                if prog.match(fname):
                    self.add(fname)

class anal_merge(object):

    def add_to_arrays(self,date,values,nvalues):

        if not self.init_single:
            self.init_single = True
            self.stored_time_single.append(date)
            self.stored_value_single[0,:] = values
            self.stored_nvalue_single[0,:] = values
        else:
            self.stored_time_single.append(date)
            self.stored_value_single = np.append(self.stored_value_single,\
                                                     values[np.newaxis,:],\
                                                     axis=0)
            self.stored_nvalue_single = np.append(self.stored_nvalue_single,\
                                                      nvalues[np.newaxis,:],\
                                                      axis=0)

    def read_single(self,filename,channel,day=True,week=False):

        self.stored_time_single=[]
        self.stored_value_single=np.zeros((1,18))
        self.stored_nvalue_single=np.zeros((1,18),dtype=np.int32)
        self.init_single=False

        ncid = nc.Dataset(filename)
            
        date = nc.num2date(ncid.variables['time'][:],\
                               ncid.variables['time'].units)
        
        values = ncid.variables['map_data'][:,channel-1,:]
        nvalues = ncid.variables['nmap_data'][:,channel-1,:]                   

        for i in range(values.shape[0]):
            self.add_to_arrays(date[i],values[i,:],nvalues[i,:])

        ncid.close()

    def check_date(self,old_day,old_week,time,day=True,week=False):

        if day:
            if time.timetuple().tm_yday == old_day:
                return True
            else:
                return False
        else:
            raise Exception,"Week averaging not yet supported"

    def plot_merge(self,filelist,channel,day=True,week=False,latpos=-1,\
                       yrange=[-1,-1]):    

        for filename in filelist:
            self.read_single(filename,channel,day=day,week=week)
            init=False
            average = np.zeros(1)
            naverage = np.zeros(1,dtype=np.int32)
            timeval = []
            old_day = -1
            old_week = -1
            nday=0
            for i in range(self.stored_value_single.shape[0]):
                if self.check_date(old_day,old_week,\
                                       self.stored_time_single[i],\
                                       day=day,week=week):
                    if -1 == latpos:
                        naverage[nday] = naverage[nday]+\
                            np.sum(self.stored_nvalue_single[i,:])
                        average[nday] = average[nday]+\
                            np.sum(self.stored_value_single[i,:])
                    else:
                        print latpos,i,self.stored_nvalue_single.shape
                        naverage[nday] = naverage[nday]+\
                            self.stored_nvalue_single[i,latpos]
                        average[nday] = average[nday]+\
                            self.stored_value_single[i,latpos]
                else:
                    old_day = self.stored_time_single[i].timetuple().tm_yday
                    old_week = self.stored_time_single[i].timetuple().tm_wday
   
                    if -1 == latpos:
                        if init:
                            naverage = np.append(naverage,\
                                                     np.sum(self.stored_nvalue_single[i,:]))
                            average = np.append(average,\
                                                    np.sum(self.stored_value_single[i,:]))
                        else:
                            naverage[0] = np.sum(self.stored_nvalue_single[i,:])
                            average[0] = np.sum(self.stored_value_single[i,:])
                    else:
                        if init:
                            naverage = np.append(naverage,\
                                                     self.stored_nvalue_single[i,latpos])
                            average = np.append(average,\
                                                    self.stored_value_single[i,latpos])
                        else:
                            naverage[0] = self.stored_nvalue_single[i,latpos]
                            average[0] = self.stored_value_single[i,latpos]
                    if not init:
                        init=True
                    else:
                        nday = nday + 1                        
                    timeval.append(\
                        datetime.datetime(\
                            self.stored_time_single[i].year,\
                            self.stored_time_single[i].month,\
                            self.stored_time_single[i].day))

            gd = (naverage > 0)
            if np.sum(gd) > 0:
                average[gd] = average[gd]/naverage[gd]

            gd = (average > 180)
            if np.sum(gd) > 0:
                newtime=[]
                for i in range(len(gd)):
                    if gd[i]:
                        newtime.append(timeval[i])
            plt.plot(newtime,average[gd],'.')            

        if yrange[0] != yrange[1]:
            plt.ylim(yrange)
        if 4 == channel:
            plt.title('3.7$\mu$m')
        elif 5 == channel:
            plt.title('11$\mu$m')
        elif 6 == channel:
            plt.title('12$\mu$m')
        plt.show()

    def __init__(self,filelist,channel,day=True,week=False,latpos=-1,\
                     yrange=[-1,-1]):

        self.plot_merge(filelist,channel,day=day,week=week,latpos=latpos,\
                            yrange=yrange)

def plot_all(channel,yrange=[-1,-1]):

    d = anal_merge(['avhrr11_stats.nc',\
                        'avhrr12_stats.nc',\
                        'avhrr14_stats.nc',\
                        'avhrr15_stats.nc',\
                        'avhrr16_stats.nc',\
                        'avhrr17_stats.nc',\
                        'avhrr18_stats.nc',\
                        'avhrr19_stats.nc',\
                        'avhrrmta_stats.nc'],channel,yrange=yrange)

class read_uncertainty_file(object):

    def read(self,filename):

        ncid = nc.Dataset(filename,'r')

        self.time = nc.num2date(ncid.variables['time'][:],\
                                    ncid.variables['time'].units)
        self.ch1_rand = ncid.variables['ch1_random'][:,:]        
        self.ch2_rand = ncid.variables['ch2_random'][:,:]
        self.ch3a_rand = ncid.variables['ch3a_random'][:,:]
        self.ch3b_rand = ncid.variables['ch3b_random'][:,:]
        self.ch4_rand = ncid.variables['ch4_random'][:,:]
        self.ch5_rand = ncid.variables['ch5_random'][:,:]
        self.ch1_str = ncid.variables['ch1_structured'][:,:]
        self.ch2_str = ncid.variables['ch2_structured'][:,:]
        self.ch3a_str = ncid.variables['ch3a_structured'][:,:]
        self.ch3b_str = ncid.variables['ch3b_structured'][:,:]
        self.ch4_str = ncid.variables['ch4_structured'][:,:]
        self.ch5_str = ncid.variables['ch5_structured'][:,:]
        self.ch1_common = ncid.variables['ch1_common'][:,:]
        self.ch2_common = ncid.variables['ch2_common'][:,:]
        self.ch3a_common = ncid.variables['ch3a_common'][:,:]
        self.ch3b_common = ncid.variables['ch3b_common'][:,:]
        self.ch4_common = ncid.variables['ch4_common'][:,:]
        self.ch5_common = ncid.variables['ch5_common'][:,:]

        ncid.close()

    def __init__(self,filename):

        self.read(filename)

def plot_uncertainties(filename,avhrr_type,chan=1,plotall=True,mean=True,\
                           maxval=False):

    if maxval:
        pos = 1
    else:
        if mean:
            pos = 2
        else:
            pos = 5

    max_val=1e35
    min_val=-1e35
    if 1 == chan:
        if maxval:
            ytitle = ' Max 0.6$\mu$m'
        else:
            if mean:
                ytitle = ' Mean 0.6$\mu$m'
            else:
                ytitle = ' Median 0.6$\mu$m'
        if not plotall:
            min_val=0.
            max_val=1.5
    elif 2 == chan:
        if maxval:
            ytitle = ' Max 0.8$\mu$m'
        else:
            if mean:
                ytitle = ' Mean 0.8$\mu$m'
            else:
                ytitle = ' Median 0.8$\mu$m'
        if not plotall:
            min_val=0.
            max_val=1.5
    elif 3 == chan:
        if maxval:
            ytitle = ' Max 1.6$\mu$m'
        else:
            if mean:
                ytitle = ' Mean 1.6$\mu$m'
            else:
                ytitle = ' Median 1.6$\mu$m'
        if not plotall:
            min_val=0.
            max_val=1.5
    elif 4 == chan:
        if maxval:
            ytitle = ' Max 3.7$\mu$m'
        else:
            if mean:
                ytitle = ' Mean 3.7$\mu$m'
            else:
                ytitle = ' Median 3.7$\mu$m'
        if not plotall:
            min_val=0.
            max_val=100.
    elif 5 == chan:
        if maxval:
            ytitle = ' Max 11$\mu$m'
        else:
            if mean:
                ytitle = ' Mean 11$\mu$m'
            else:
                ytitle = ' Median 11$\mu$m'
        if not plotall:
            min_val=0.
            max_val=100.
    elif 6 == chan:
        if maxval:
            ytitle = ' Max 12$\mu$m'
        else:
            if mean:
                ytitle = ' Mean 12$\mu$m'
            else:
                ytitle = ' Median 12$\mu$m'
        if not plotall:
            min_val=0.
            max_val=100.

    if 11 == avhrr_type:
        title = 'NOAA11'
    elif 12 == avhrr_type:
        title = 'NOAA12'
    elif 14 == avhrr_type:
        title = 'NOAA14'
    elif 15 == avhrr_type:
        title = 'NOAA15'
    elif 16 == avhrr_type:
        title = 'NOAA16'
    elif 17 == avhrr_type:
        title = 'NOAA17'
    elif 18 == avhrr_type:
        title = 'NOAA18'
    elif 19 == avhrr_type:
        title = 'NOAA19'
    elif -1 == avhrr_type:
        title = 'MetOp-A'
    else:
        raise Exception,"avhrr_type not found: plot_uncertainties"

    d = read_uncertainty_file(filename)

    if 1 == chan:
        u_rand = d.ch1_rand[:,pos]
        u_str = d.ch1_str[:,pos]
        u_common = d.ch1_common[:,pos]
    elif 2 == chan:
        u_rand = d.ch2_rand[:,pos]
        u_str = d.ch2_str[:,pos]
        u_common = d.ch2_common[:,pos]
    elif 3 == chan:
        u_rand = d.ch3a_rand[:,pos]
        u_str = d.ch3a_str[:,pos]
        u_common = d.ch3a_common[:,pos]
    elif 4 == chan:
        u_rand = d.ch3b_rand[:,pos]
        u_str = d.ch3b_str[:,pos]
        u_common = d.ch3b_common[:,pos]
    elif 5 == chan:
        u_rand = d.ch4_rand[:,pos]
        u_str = d.ch4_str[:,pos]
        u_common = d.ch4_common[:,pos]
    elif 6 == chan:
        u_rand = d.ch5_rand[:,pos]
        u_str = d.ch5_str[:,pos]
        u_common = d.ch5_common[:,pos]
    gd = (u_rand < max_val) & (u_str < max_val) & (u_common < max_val) & \
        (u_rand >= min_val) & (u_str >= min_val) & (u_common >= min_val)
    if np.sum(gd) > 0:
        time=[]
        urand=[]
        ustr=[]
        ucommon=[]
        minyear=3000
        maxyear=0
        for i in range(len(gd)):
            if gd[i]:
                minyear = min([d.time[i].year,minyear])
                maxyear = max([d.time[i].year,maxyear])
                time.append(d.time[i])
                urand.append(u_rand[i])
                ustr.append(u_str[i])
                ucommon.append(u_common[i])
        urand=np.array(urand)
        ustr=np.array(ustr)
        ucommon=np.array(ucommon)
        
        plt.plot(time,urand,',',color='red')
        plt.plot(time,ustr,',',color='blue')
        plt.plot(time,ucommon,',',color='green')
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        if maxyear-minyear > 8:
            plt.gca().xaxis.set_major_locator(mdates.YearLocator(4))
        else:
            plt.gca().xaxis.set_major_locator(mdates.YearLocator(3))
        plt.ylabel(ytitle)
        plt.title(title)

def uncert_plots(chan=1,plotall=True,mean=True,maxval=False):

    plt.subplot(331)
    plot_uncertainties('avhrr11_stats.nc',11,chan=chan,plotall=plotall,\
                           mean=mean,maxval=maxval)
    plt.subplot(332)
    plot_uncertainties('avhrr12_stats.nc',12,chan=chan,plotall=plotall,\
                           mean=mean,maxval=maxval)
    plt.subplot(333)
    plot_uncertainties('avhrr14_stats.nc',14,chan=chan,plotall=plotall,\
                           mean=mean,maxval=maxval)
    plt.subplot(334)
    plot_uncertainties('avhrr15_stats.nc',15,chan=chan,plotall=plotall,\
                           mean=mean,maxval=maxval)
    plt.subplot(335)
    plot_uncertainties('avhrr16_stats.nc',16,chan=chan,plotall=plotall,\
                           mean=mean,maxval=maxval)
    plt.subplot(336)
    plot_uncertainties('avhrr17_stats.nc',17,chan=chan,plotall=plotall,\
                           mean=mean,maxval=maxval)
    plt.subplot(337)
    plot_uncertainties('avhrr18_stats.nc',18,chan=chan,plotall=plotall,\
                           mean=mean,maxval=maxval)
    plt.subplot(338)
    plot_uncertainties('avhrr19_stats.nc',19,chan=chan,plotall=plotall,\
                           mean=mean,maxval=maxval)
    plt.subplot(339)
    plot_uncertainties('avhrrmta_stats.nc',-1,chan=chan,plotall=plotall,\
                           mean=mean,maxval=maxval)
    plt.tight_layout()
    plt.show()


if "__main__" == __name__:
    
    parser = OptionParser("usage: %prog uuid outfile netcdf_out")
    (options, args) = parser.parse_args()

    if len(args) != 3:
        parser.error("incorrect number of arguments")

    if 'Y' == args[2] or 'y' ==args[2]:
        netcdf_out=True
    else:
        netcdf_out=False

    get_stats(args[0],args[1],netcdf_out)

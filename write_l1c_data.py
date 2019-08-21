# * Copyright (C) 2018 J.Mittaz University of Reading
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
#
# Write GBCS Level 1C format data
#
#
# Write l1c float values
#
# Routines to write all metadata and scale etc. if needed
#
def Write_GBCS_Float(ncid,name,dim_nx,dim_ny,indata,fill_value,\
                         long_name,standard_name,units,\
                         valid_min,valid_max,reference_datum=None):

    newdata = data
    gd = (newdata < valid_min) | (newdata > valid_max)
    if np.sum(gd) > 0:
        newdata[gd] = fill_value

    var = ncid.createVariable(name,'f4',(dim_ny,dim_nx),zlib=False,\
                                  complevel=5,shuffle=True,\
                                  fill_value=fill_value)
    var.long_name = long_name
    var.standard_name = standard_name
    var.units = units
    var.valid_min = valid_min
    var.valid_max = valid_max
    if reference_datum:
        var.reference_datum = reference_datum

    var[:,:] = newdata

def Write_GBCS_Float_1d(ncid,name,dim_ny,indata,fill_value,\
                            long_name,standard_name,units,valid_min,valid_max):

    newdata = data
    gd = (newdata < valid_min) | (newdata > valid_max)
    if np.sum(gd) > 0:
        newdata[gd] = fill_value

    var = ncid.createVariable(name,'f4',(dim_ny),zlib=False,\
                                  complevel=5,shuffle=True,\
                                  fill_value=fill_value)
    var.long_name = long_name
    var.standard_name = standard_name
    var.units = units
    var.valid_min = valid_min
    var.valid_max = valid_max

    var[:] = newdata

def Write_GBCS_time(ncid,dim_time,data,units):

    newdata = nc.date2num(data,units)
    var = ncid.createVariable(name,'i4',(dim_time),zlib=False,\
                                  complevel=5,shuffle=True)

    var.long_name = 'reference time of sst file'
    var.standard_name = 'time'
    var.units = units
    var.calendar = 'gregorian'

    var[:] = newdata

def Write_GBCS_dtime(ncid,dim_nx,dim_ny,nx,ny,dim_time,time_step,\
                         long_name,units,coordinates,fill_value):

    data = np.zeros((1,ny,nx),dtype=np.float32)
    for i in range(ny):
        data[0,i,:] = i*time_step
    var = ncid.createVariable(name,'f4',(dim_time,dim_ny,dim_nx),zlib=False,\
                                  complevel=5,shuffle=True,\
                                  fill_value=fill_value)
    var.long_name = long_name
    var.units = units
    var.coordinates = coordinates

    var[:,:,:] = data

def Write_GBCS_Float_Time(ncid,name,dim_ny,dim_time,ny,\
                              indata,fill_value,long_name,units,\
                              valid_min,valid_max,scale,offset):

    data = np.zeros((1,ny),dtype=np.int16)
    gd = (indata >= valid_min_f) & (indata <= valid_max_f)
    if np.sum(gd) > 0:
        data[0,gd] = ((indata[gd]-offset)/scale).astype(dtype=np.int16)
    gd = (indata < valid_min_f) | (indata > valid_max_f)
    if np.sum(gd) > 0:
        data[0,gd] = fill_value

    var = ncid.createVariable(name,'i2',(time_time,dim_ny),zlib=False,\
                                  complevel=5,shuffle=True,\
                                  fill_value=fill_value)
    var.long_name = long_name
    var.standard_name = standard_name
    var.units = units
    var.valid_min = valid_min
    var.valid_max = valid_max
    var.add_offset = offset
    var.scale_factor = scale_factor

    var[:,:] = data

def Write_GBCS_Float_to_Int(ncid,name,dim_nx,dim_ny,nx,ny,\
                                indata,fill_value,long_name,\
                                standard_name,units,valid_min,valid_max,\
                                scale,offset,coordinates=None,\
                                comment=None):

    data = np.zeros((ny,nx),dtype=np.int16)
    gd = (indata >= valid_min_f) & (indata <= valid_max_f)
    if np.sum(gd) > 0:
        data[gd] = ((indata[gd]-offset)/scale).astype(dtype=np.int16)
    gd = (indata < valid_min_f) | (indata > valid_max_f)
    if np.sum(gd) > 0:
        data[gd] = fill_value

    var = ncid.createVariable(name,'i2',(dim_ny,dim_nx),zlib=False,\
                                  complevel=5,shuffle=True,\
                                  fill_value=fill_value)
    var.long_name = long_name
    var.standard_name = standard_name
    var.units = units
    var.valid_min = valid_min
    var.valid_max = valid_max
    var.add_offset = offset
    var.scale_factor = scale_factor
    if coordinates:
        var.coordinates = coordinates
    if comment:
        var.comment = comment

    var[:,:] = data

def Write_GBCS_Float_to_Int_3D(ncid,name,dim_nx,dim_ny,dim_time,nx,ny,\
                                   indata,fill_value,long_name,\
                                   units,valid_min,valid_max,\
                                   scale,offset,standard_name=None,\
                                   coordinates=None,\
                                   comment=None):

    data = np.zeros((ny,nx),dtype=np.int16)
    gd = (indata >= valid_min_f) & (indata <= valid_max_f)
    if np.sum(gd) > 0:
        data[gd] = ((indata[gd]-offset)/scale).astype(dtype=np.int16)
    gd = (indata < valid_min_f) | (indata > valid_max_f)
    if np.sum(gd) > 0:
        data[gd] = fill_value

    var = ncid.createVariable(name,'i2',(dim_time,dim_ny,dim_nx),zlib=False,\
                                  complevel=5,shuffle=True,\
                                  fill_value=fill_value)
    var.long_name = long_name
    if standard_name:
        var.standard_name = standard_name
    var.units = units
    var.valid_min = valid_min
    var.valid_max = valid_max
    var.add_offset = offset
    var.scale_factor = scale_factor
    if coordinates:
        var.coordinates = coordinates
    if comment:
        var.comment = comment

    var[:,:,:] = data[np.newaxis,:,:]

def Write_GBCS_Byte(ncid,name,dim_nx,dim_ny,dim_time,nx,ny,\
                        indata,fill_value,long_name,flag_meanings,nflags,\
                        flag_values,valid_min,valid_max,comment=None,\
                        scale=None,offset=None):

    data = np.zeros((1,ny,nx),dtype=np.int8)
    for i in range(ny):
        data[0,i,:] = indata[i]

    var.long_name = long_name
    var.standard_name = standard_name
    var.units = units
    var.valid_min = valid_min
    var.valid_max = valid_max
    if offset: 
        var.add_offset = offset
    if scale:
        var.scale_factor = scale
    if comment:
        var.comment = comment    

def ibset_flags(value,bitpos):

    bit = 1 << bitpos
    value &= bit

    return value

def make_flags(data):

    nflags = 7
    flags_meanings = 'bad_navigation bad_calibration bad_timing missing_line solar_contamination_3B solar_contamination_4 solar_contamination_5'
    valid_min = 0
    flags = np.zeros(data.lat.shape[0],dtype=np.int8)
    flags_values = np.zeros(data.lat.shape[0],dtype=np.int8)
    for i in range(nflags):
        flag_values[I] = IBSET(0,I-1)       

    valid_max = flag_values[nflags-1]

    flags = 0
    for i in range(data.lat.shape[0]):
       if data.badNav[i] == 1:
           flags[i] = ibset_flag(flags[i],0)
       if data.badCal[i] == 1:
           flags[i] = ibset_flag(flags[i],1)
       if data.badTime[i] == 1:
           flags[i] = ibset_flag(flags[i],2)
       if data.missingLines[i] == 1:
           flags[i] = ibset_flag(flags[i],3)
       if data.solar3[i] == 1:
           flags[i] = ibset_flag(flags[i],4)
       if data.solar4[i] == 1:
           flags[i] = ibset_flag(flags[i],5)
       if data.solar5[i] == 1:
           flags[i] = ibset_flag(flags[i],6)

    return nflags,flags,flags_meaning,flag_values,valid_min,valid_max

#
# Write GBCS L1C output including the channel to channel covariances
#
def write_gbcs_l1c(data,S_s,himawari=False,fiduceo=False):

    name = '{0:04d}{1:02d}{2:02d}{3:02d}{4:02d}{5:02d}-ESACCI-L1C-{6}-fv01.0.nc'.\
        format(data.date_time[0].year,\
                   data.date_time[0].month,\
                   data.date_time[0].day,\
                   data.date_time[0].hour,\
                   data.date_time[0].minute,\
                   data.date_time[0].second,\
                   data.noaa_string)

    ncid = nc.Dataset(name,'w',format='NETCDF4')

    # Create dimensions
    dim_ni = ncid.createDimension('ni',size=data.lat.shape[1])
    dim_nj = ncid.createDimension('nj',size=data.lat.shape[0])
    dim_time = ncid.createDimension('time',size=1)
    dim_band = ncid.createDimension('band',size=3)

    # Create and wrie variables
    Write_GBCS_Float(ncid,'lat','ni','nj',data.lat,\
                         -32768.,'Latitude coordinates','latitude',\
                         'degrees_north',-90.,90.,\
                         'geographical coordinates, WGS84 projection')
    lon = data.lon
    gd = (lon > 180.)
    if np.sum(gd) > 0:
        lon[gd] = lon-360.
    Write_GBCS_Float(ncid,'lon','ni','nj',lon,\
                         -32768.,'Longitude coordinates','longitude',\
                         'degrees_east',-180.,180.,\
                         'geographical coordinates, WGS84 projection')

    Write_GBCS_time(ncid,'time',data.date_time[0],\
                        'seconds since 1981-01-01 00:00:00')
    Write_GBCS_dtime(ncid,'ni','nj','time',data.date_time[0],\
                         0.5,'scanline time difference from start time',\
                         'seconds','lon lat',-32768.)

    Write_GBCS_Float_to_Int_3D(ncid,'ch1','ni','nj','time',\
                                   data.ch1,-32768,'Channel 1 Reflectance',\
                                   'reflectance',0,15000,\
                                   0.0001,0.,coordinates='lon lat')
    
    Write_GBCS_Float_to_Int_3D(ncid,'ch2','ni','nj','time',\
                                   data.ch2,-32768,'Channel 2 Reflectance',\
                                   'reflectance',0,15000,\
                                   0.0001,0.,coordinates='lon lat')
    
    if data.ch3a_there:
        Write_GBCS_Float_to_Int_3D(ncid,'ch3a','ni','nj','time',\
                                       data.ch3a,-32768,\
                                       'Channel 3A Reflectance',\
                                       'reflectance',0,15000,\
                                       0.0001,0.,coordinates='lon lat')
    
    Write_GBCS_Float_to_Int_3D(ncid,'ch3b','ni','nj','time',\
                                   data.ch3b,-32768,\
                                   'Channel 3B Brightness Temperature',\
                                   'kelvin',-20000,10000,\
                                   0.01,273.15,coordinates='lon lat')

    Write_GBCS_Float_to_Int_3D(ncid,'ch4','ni','nj','time',\
                                   data.ch4,-32768,\
                                   'Channel 4 Brightness Temperature',\
                                   'kelvin',-20000,10000,\
                                   0.01,273.15,coordinates='lon lat')
    
    if data.ch5_there:
        Write_GBCS_Float_to_Int_3D(ncid,'ch5','ni','nj','time',\
                                       data.ch5,-32768,\
                                       'Channel 5 Brightness Temperature',\
                                       'kelvin',-20000,10000,\
                                       0.01,273.15,coordinates='lon lat')
    
    
    #
    # Now noise estimates - note here we combine independent and structured
    # from FIDUCEO information
    #
    noise = np.zeros(data.ch1.shape,dtype=np.float32)
    noise[:,:] = -1e30
    gd = (data.ch1 > 0)&(data.u_random_Ch1 > 0)&\
        (data.u_non_random_Ch1 > 0)&\
        (data.u_common_Ch1 > 0)
    if np.sum(gd) > 0:
        noise[gd] = np.sqrt((data.u_common_ch1[gd]*data.ch1[gd])**2 +\
                                (data.u_random_ch1[gd]**2) + \
                                (data.u_non_random_ch1[gd]**2))
    Write_GBCS_Float_to_Int_3D(ncid,'ch1_noise','ni','nj','time',noise,\
                                   -32768,'Channel 1 noise estimate',\
                                   'reflectance',0,10000,1.e-5,0.,\
                                   coordinates='lon lat')
    
    noise[:,:] = -1e30
    gd = (data.ch2 > 0)&(data.u_random_Ch2 > 0)&(data.u_non_random_Ch2 > 0)
    if np.sum(gd) > 0:
        noise[gd] = np.sqrt((data.u_common_ch2[gd]*data.ch2[gd])**2 +\
                                (data.u_random_ch2[gd]**2)+\
                                (data.u_non_random_ch2[gd]**2))
    Write_GBCS_Float_to_Int_3D(ncid,'ch2_noise','ni','nj','time',noise,\
                                   -32768,'Channel 2 noise estimate',\
                                   'reflectance',0,10000,1.e-5,0.,\
                                   coordinates='lon lat')
    if data.ch3a_there:
        noise[:,:] = -1e30
        gd = (data.ch3a > 0)&(data.u_random_Ch3a > 0)&\
            (data.u_non_random_Ch3a > 0)
        if np.sum(gd) > 0:
            noise[gd] = np.sqrt((data.u_common_ch3a[gd]*\
                                     data.ch3a[gd])**2 +\
                                    (data.u_random_ch3a[gd]**2)+\
                                    (data.u_non_random_ch3a[gd]**2))
            Write_GBCS_Float_to_Int_3D(ncid,'ch3a_noise','ni','nj','time',\
                                           noise,-32768,\
                                           'Channel 3A noise estimate',\
                                           'reflectance',0,10000,1.e-5,0.,\
                                           coordinates='lon lat')
    
    noise = np.zeros(data.ch3b.shape,dtype=np.float32)
    noise[:,:] = -1e30
    gd = (data.ch3b > 0)&(data.u_random_Ch3b > 0)&(data.u_non_random_Ch3b > 0)
    if np.sum(gd) > 0:
        noise[gd] = np.sqrt((data.u_non_random_ch3b[gd])**2 +\
                                (data.u_random_ch3b[gd]**2)+\
                                (data.u_common_ch3b[gd]**2))
    Write_GBCS_Float_to_Int_3D(ncid,'ch3b_nedt','ni','nj','time',noise,\
                                   -32768,'Channel 3B noise estimate',\
                                   'kelvin',0,10000,0.001,0.,\
                                   coordinates='lon lat')
    
    noise = np.zeros(data.ch4.shape,dtype=np.float32)
    noise[:,:] = -1e30
    gd = (data.ch4 > 0)&(data.u_random_Ch4 > 0)&(data.u_non_random_Ch4 > 0)
    if np.sum(gd) > 0:
        noise[gd] = np.sqrt((data.u_non_random_ch4[gd])**2 +\
                                (data.u_random_ch4[gd]**2) +\
                                (data.u_common_ch4[gd]**2))
    Write_GBCS_Float_to_Int_3D(ncid,'ch4_nedt','ni','nj','time',noise,\
                                   -32768,'Channel 4 noise estimate',\
                                   'kelvin',0,10000,0.001,0.,\
                                   coordinates='lon lat')
    
    if data.ch5_there:
        noise = np.zeros(data.ch4.shape,dtype=np.float32)
        noise[:,:] = -1e30
        gd = (data.ch5 > 0)&(data.u_random_Ch5 > 0)&(data.u_non_random_Ch5 > 0)
        if np.sum(gd) > 0:
            noise[gd] = np.sqrt((data.u_non_random_ch5[gd]**2) +\
                                (data.u_random_ch5[gd]**2) + \
                                (data.u_common_ch5[gd]**2))
        Write_GBCS_Float_to_Int_3D(ncid,'ch5_nedt','ni','nj','time',noise,\
                                       -32768,'Channel 5 noise estimate',\
                                       'kelvin',0,10000,0.001,0.,\
                                       coordinates='lon lat')
    

    if fiduceo:
        #
        # Now provide separate uncertainties for more complex uncertainties
        # for IR channels
        #
        noise = np.zeros(data.ch1.shape,dtype=np.float32)
        noise[:,:] = -1e30
        gd = (data.ch1 > 0)&(data.u_random_Ch1 > 0)
        if np.sum(gd) > 0:
            noise[gd] = data.u_random_ch1[gd]
        Write_GBCS_Float_to_Int_3D(ncid,'ch1_u_random','ni','nj','time',noise,\
                                       -32768,'Channel 1 random uncertainty estimate',\
                                       'reflectance',0,10000,1.e-5,0.,\
                                       coordinates='lon lat')
    
        noise[:,:] = -1e30
        gd = (data.ch2 > 0)&(data.u_random_Ch2 > 0)
        if np.sum(gd) > 0:
            noise[gd] = data.u_random_ch2[gd]
        Write_GBCS_Float_to_Int_3D(ncid,'ch2_u_random','ni','nj','time',noise,\
                                       -32768,'Channel 2 random uncertainty estimate',\
                                       'reflectance',0,10000,1.e-5,0.,\
                                       coordinates='lon lat')
        if data.ch3a_there:
            noise[:,:] = -1e30
            gd = (data.ch3a > 0)&(data.u_random_Ch3a > 0)
            if np.sum(gd) > 0:
                noise[gd] = data.u_random_ch3a[gd]
                Write_GBCS_Float_to_Int_3D(ncid,'ch3a_u_random','ni','nj','time',\
                                               noise,-32768,\
                                               'Channel 3A random uncertainty estimate',\
                                               'reflectance',0,10000,1.e-5,0.,\
                                               coordinates='lon lat')
    
        noise = np.zeros(data.ch3b.shape,dtype=np.float32)
        noise[:,:] = -1e30
        gd = (data.ch3b > 0)&(data.u_random_Ch3b > 0)
        if np.sum(gd) > 0:
            noise[gd] = data.u_random_ch3b[gd]
        Write_GBCS_Float_to_Int_3D(ncid,'ch3b_u_random','ni','nj','time',noise,\
                                       -32768,'Channel 3B random uncertainty estimate',\
                                       'kelvin',0,10000,0.001,0.,\
                                       coordinates='lon lat')
    
        noise = np.zeros(data.ch4.shape,dtype=np.float32)
        noise[:,:] = -1e30
        gd = (data.ch4 > 0)&(data.u_random_Ch4 > 0)
        if np.sum(gd) > 0:
            noise[gd] = data.u_random_ch4[gd]
        Write_GBCS_Float_to_Int_3D(ncid,'ch4_u_random','ni','nj','time',noise,\
                                       -32768,'Channel 4 random uncertainty estimate',\
                                       'kelvin',0,10000,0.001,0.,\
                                       coordinates='lon lat')
    
        if data.ch5_there:
            noise = np.zeros(data.ch4.shape,dtype=np.float32)
            noise[:,:] = -1e30
            gd = (data.ch5 > 0)&(data.u_random_Ch5 > 0)
            if np.sum(gd) > 0:
                noise[gd] = data.u_random_ch5[gd]
            Write_GBCS_Float_to_Int_3D(ncid,'ch5_u_random','ni','nj','time',noise,\
                                           -32768,'Channel 5 random uncertainty estimate',\
                                           'kelvin',0,10000,0.001,0.,\
                                           coordinates='lon lat')
    
        noise = np.zeros(data.ch1.shape,dtype=np.float32)
        noise[:,:] = -1e30
        gd = (data.ch1 > 0)&(data.u_non_random_Ch1 > 0)
        if np.sum(gd) > 0:
            noise[gd] = data.u_non_random_ch1[gd]
        Write_GBCS_Float_to_Int_3D(ncid,'ch1_u_non_random','ni','nj','time',noise,\
                                           -32768,'Channel 1 non random uncertainty estimate',\
                                           'reflectance',0,10000,1.e-5,0.,\
                                           coordinates='lon lat')
    
        noise[:,:] = -1e30
        gd = (data.ch2 > 0)&(data.u_non_random_Ch2 > 0)
        if np.sum(gd) > 0:
            noise[gd] = data.u_non_random_ch2[gd]
        Write_GBCS_Float_to_Int_3D(ncid,'ch2_u_non_random','ni','nj','time',noise,\
                                       -32768,'Channel 2 non random uncertainty estimate',\
                                       'reflectance',0,10000,1.e-5,0.,\
                                       coordinates='lon lat')
        if data.ch3a_there:
            noise[:,:] = -1e30
            gd = (data.ch3a > 0)&(data.u_non_random_Ch3a > 0)
            if np.sum(gd) > 0:
                noise[gd] = data.u_non_random_ch3a[gd]
            Write_GBCS_Float_to_Int_3D(ncid,'ch3a_u_non_random','ni','nj','time',\
                                           noise,-32768,\
                                           'Channel 3A non random uncertainty estimate',\
                                           'reflectance',0,10000,1.e-5,0.,\
                                           coordinates='lon lat')
    
        noise = np.zeros(data.ch3b.shape,dtype=np.float32)
        noise[:,:] = -1e30
        gd = (data.ch3b > 0)&(data.u_non_random_Ch3b > 0)
        if np.sum(gd) > 0:
            noise[gd] = data.u_non_random_ch3b[gd]
        Write_GBCS_Float_to_Int_3D(ncid,'ch3b_u_non_random','ni','nj','time',noise,\
                                       -32768,'Channel 3B non random uncertainty estimate',\
                                       'kelvin',0,10000,0.001,0.,\
                                       coordinates='lon lat')
    
        noise = np.zeros(data.ch4.shape,dtype=np.float32)
        noise[:,:] = -1e30
        gd = (data.ch4 > 0)&(data.u_non_random_Ch4 > 0)
        if np.sum(gd) > 0:
            noise[gd] = data.u_non_random_ch4[gd]
        Write_GBCS_Float_to_Int_3D(ncid,'ch4_u_non_random','ni','nj','time',noise,\
                                       -32768,'Channel 4 non random uncertainty estimate',\
                                       'kelvin',0,10000,0.001,0.,\
                                       coordinates='lon lat')
        
        if data.ch5_there:
            noise = np.zeros(data.ch4.shape,dtype=np.float32)
            noise[:,:] = -1e30
            gd = (data.ch5 > 0)&(data.u_non_random_Ch5 > 0)
            if np.sum(gd) > 0:
                noise[gd] = data.u_non_random_ch5[gd]
            Write_GBCS_Float_to_Int_3D(ncid,'ch5_u_non_random','ni','nj','time',noise,\
                                           -32768,'Channel 5 non random uncertainty estimate',\
                                           'kelvin',0,10000,0.001,0.,\
                                           coordinates='lon lat')
    

        noise = np.zeros(data.ch1.shape,dtype=np.float32)
        noise[:,:] = -1e30
        gd = (data.ch1 > 0)&(data.u_common_Ch1 > 0)
        if np.sum(gd) > 0:
            noise[gd] = data.u_common_ch1[gd]*data.ch1[gd]
        Write_GBCS_Float_to_Int_3D(ncid,'ch1_u_common','ni','nj','time',noise,\
                                       -32768,'Channel 1 common uncertainty estimate',\
                                       'reflectance',0,10000,1.e-5,0.,\
                                       coordinates='lon lat')
    
        noise[:,:] = -1e30
        gd = (data.ch2 > 0)&(data.u_common_Ch2 > 0)
        if np.sum(gd) > 0:
            noise[gd] = data.u_common_ch2[gd]*data.ch2[gd]
        Write_GBCS_Float_to_Int_3D(ncid,'ch2_u_common','ni','nj','time',noise,\
                                       -32768,'Channel 2 common uncertainty estimate',\
                                       'reflectance',0,10000,1.e-5,0.,\
                                       coordinates='lon lat')
        if data.ch3a_there:
            noise[:,:] = -1e30
            gd = (data.ch3a > 0)&(data.u_common_Ch3a > 0)
            if np.sum(gd) > 0:
                noise[gd] = data.u_common_ch3a[gd]*data.cha[gd]
            Write_GBCS_Float_to_Int_3D(ncid,'ch3a_u_common','ni','nj','time',\
                                           noise,-32768,\
                                           'Channel 3A common uncertainty estimate',\
                                           'reflectance',0,10000,1.e-5,0.,\
                                           coordinates='lon lat')
    
        noise = np.zeros(data.ch3b.shape,dtype=np.float32)
        noise[:,:] = -1e30
        gd = (data.ch3b > 0)&(data.u_common_Ch3b > 0)
        if np.sum(gd) > 0:
            noise[gd] = data.u_common_ch3b[gd]
        Write_GBCS_Float_to_Int_3D(ncid,'ch3b_u_common','ni','nj','time',noise,\
                                       -32768,'Channel 3B common uncertainty estimate',\
                                       'kelvin',0,10000,0.001,0.,\
                                       coordinates='lon lat')
    
        noise = np.zeros(data.ch4.shape,dtype=np.float32)
        noise[:,:] = -1e30
        gd = (data.ch4 > 0)&(data.u_common_Ch4 > 0)
        if np.sum(gd) > 0:
            noise[gd] = data.u_common_ch4[gd]
        Write_GBCS_Float_to_Int_3D(ncid,'ch4_u_common','ni','nj','time',noise,\
                                       -32768,'Channel 4 non common uncertainty estimate',\
                                       'kelvin',0,10000,0.001,0.,\
                                       coordinates='lon lat')
    
        if data.ch5_there:
            noise = np.zeros(data.ch4.shape,dtype=np.float32)
            noise[:,:] = -1e30
            gd = (data.ch5 > 0)&(data.u_common_Ch5 > 0)
            if np.sum(gd) > 0:
                noise[gd] = data.u_common_ch5[gd]
            Write_GBCS_Float_to_Int_3D(ncid,'ch5_common_random','ni','nj','time',noise,\
                                           -32768,'Channel 5 non common uncertainty estimate',\
                                           'kelvin',0,10000,0.001,0.,\
                                           coordinates='lon lat')
    

    #
    # Angles
    #
    Write_GBCS_Float_to_Int_3D(ncid,'satellite_zenith_angle','ni','nj','time',\
                                   data.satza,-32768,'satellite zenith angle',\
                                   'angular degree',0,10000,0.01,0.,\
                                   standard_name='zenith angle',\
                                   coordinates='lon lat',\
                                   comment='The satellite zenith angle at time of the observations')

    Write_GBCS_Float_to_Int_3D(ncid,'solar_zenith_angle','ni','nj','time',\
                                   data.solza,-32768,'solar zenith angle',\
                                   'angular degree',0,18000,0.01,0.,\
                                   standard_name='zenith angle',\
                                   coordinates='lon lat',\
                                   comment='The solar zenith angle at time of the observations')

    Write_GBCS_Float_to_Int_3D(ncid,'relative_azimuth_angle','ni','nj','time',\
                                   data.solaz,-32768,'relative azimuth angle',\
                                   'angular degree',0,18000,0.01,0.,\
                                   standard_name='zenith angle',\
                                   coordinates='lon lat',\
                                   comment='The relative azimuth angle at time of the observations')

    #
    # ICT temperature
    #
    Write_GBCS_Float_Time(ncid,'ict_temp','nj','time',\
                              data.smoothPRT,-32768,\
                              'Temperature of internal calibration target',\
                              'kelvin',-20000,10000,0.01,273.15)

#
# Quality flags
#
    nflags,flags,flags_meaning,flag_values,\
        valid_min,valid_max = make_flags(data)
#    CALL Make_Flags(AVHRR,nflags,flags,flags_meaning,flag_values,&
#         valid_min,valid_max)
#    CALL Write_GBCS_Byte(ncid,'qual_flags',dimid_nx,dimid_ny,dimid_time,nx,ny,&
#       flags,-128_GbcsInt1,'Quality Flags',flags_meaning,nflags,flag_values,&
#       valid_min,valid_max,'Bad data is flagged when any one of the other &
#         &categories (bad navigation, calibration, timing or missing line) &
#         &are flagged')
    Write_GBCS_Byte(ncid,'qual_flags','ni','nj','time',\
                        flags,-128,'Quality Flags',\
                        flags_meaning,nflags,flag_values,\
                        valid_min,valid_max,\
                        'Bad data is flagged when any one of the other categories (bad navigation, calibration, timing or missing line) are flagged')

#
#    !
#    ! Cloud masks
#    !
#    flags = 0
#    flag_values(1) = 0_GbcsInt2
#    CALL Write_GBCS_Byte(ncid,'cloud_mask',dimid_nx,dimid_ny,dimid_time,nx,ny,&
#         flags,-128_GbcsInt1,'No Cloud Mask','None',1,flag_values,&
#         0_GbcsInt1,0_GbcsInt1,'')
#    CALL Write_GBCS_Byte(ncid,'cloud_probaility',dimid_nx,dimid_ny,&
#         dimid_time,nx,ny,flags,-128_GbcsInt1,'No Cloud Mask',&
#         'None',1,flag_values,0_GbcsInt1,0_GbcsInt1,'',scale=0.,offset=0.)
#    DEALLOCATE(flags,flag_values)
#
#    !
#    ! Line numbers
#    !
#    CALL Write_GBCS_Float_Time_Int(ncid,'l1b_line_number',dimid_ny,dimid_time,&
#         ny,AVHRR%scanLineNumber,-32768_GbcsInt2,&
#         'FIDUCEO Level 1 line number',0_GbcsInt2)

    # Global attributes
    if himawari:
        ncid.title = 'GBCS Pre-Processing: {0} L1C product'.\
            format('Himawari')
    else:
        ncid.title = 'AVHRR Pre-Processing: {0} L1C product (FIDUCEO based)'.\
            format(data.noaa_string)
        ncid.id = '{0}-ESACCI-L1C-vFIDUCEO'.format(data.noaa_string)
        ncid.summary = '{0} L1C product from the FIDUCEO project'.\
            format(data.noaa_string)
        ncid.reference = 'http://www.fiduceo.eu'
        ncid.institution = 'FIDUCEO'
        ncid.history = 'Created using FIDUCEO code'
        ncid.license = 'This dataset is released for use under CC-BY licence (https://creativecommons.org/licenses/by/4.0/) and was developed in the EC FIDUCEO project \"Fidelity and Uncertainty in Climate Data Records from Earth Observations\". Grant Agreement: 638822.'
    ncid.product_version = data.version
    uuidval = uuid.uuid4()
    ncid.uuid = uuidval
    ncid.tracking_id = uuidval
    ncid.netcdf_version_id = nc.getlibversion()
    ncid.date_created = datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    ncid.file_quality_level = 3
    ncid.spatial_resolution = '4.0 km at nadir'
    ncid.start_time = data.date_time[0].strftime('%Y/%m/%d %H:%M:%S')
    ncid.time_coverage_start = data.date_time[0].strftime('%Y/%m/%d %H:%M:%S')
    ncid.stop_time = data.date_time[-1].strftime('%Y/%m/%d %H:%M:%S')
    ncid.time_coverage_end = data.date_time[-1].strftime('%Y/%m/%d %H:%M:%S')
    ncid.time_coverage_duration = \
        (data.date_time[-1]-data.date_time[0]).strftime('%Y/%m/%d %H:%M:%S')
    ncid.time_coverage_resolution = \
        (data.date_time[1]-data.date_time[0]).strftime('%Y/%m/%d %H:%M:%-S')
#    ncid.source = data.source_string
    if himawari:
        ncid.platform = 'Himawari'
        ncid.sensor = 'AHI'
    else:
        ncid.platform = data.noaa_string
        ncid.sensor = 'AVHRR_GAC'
    ncid.Metadata_Conventions = 'Unidata Dataset Discovery v1.0'
    ncid.metadata_link = 'http://www.esa-cci.org'
    ncid.standard_name_vocabularly = 'NetCDF Climate and Forecast (CF) Metadata Convention'
    ncid.geospatial_lat_units = 'degrees_north'
    ncid.geospatial_lat_resolution = 0.04
    ncid.geospatial_lon_units = 'degrees_east'
    ncid.geospatial_lon_resolution = 0.04
    if not himawari:
        ncid.acknowledgement = 'NOAA GAC Data. Processing funded by H2020 (EC)'
    if himawari:
        ncid.creator_name = 'J.Mittaz University of Reading'
        ncid.creator_email = 'j.mittaz@reading.ac.uk'
        ncid.creator_url = 'http://www.rdg.ac.uk'
        ncid.creator_processing_institution = 'Bureau of Meteorology, Melbourne, Australia'
        ncid.project = 'Collaboration between UoR and the Bureau'
        ncid.northernmost_latitude = 90.
        ncid.southernmost_latitude = -90.
        ncid.easternmost_longitude = 180.
        ncid.westernmost_longitude = -180.
        ncid.geospatial_lat_min = -90.
        ncid.geospatial_lat_max = 90.
        ncid.geospatial_lon_min = -180.
        ncid.geospatial_lon_max = 180.
        ncid.processing_level = 'L1C'
    else:
        ncid.creator_name = 'FIDUCEO'
        ncid.creator_email = 'fiduceo-coordinator@lists.reading.ac.uk'
        ncid.creator_url = 'http://www.fiduceo.eu'
        ncid.creator_processing_institution = 'These data were produced on the JASMIN infrastructure at STFC as part of the H2020 FIDUCEO project'
        ncid.project = 'EC H2020 FIDUCEO project'
        ncid.northernmost_latitude = 90.
        ncid.southernmost_latitude = -90.
        ncid.easternmost_longitude = 180.
        ncid.westernmost_longitude = -180.
        ncid.geospatial_lat_min = -90.
        ncid.geospatial_lat_max = 90.
        ncid.geospatial_lon_min = -180.
        ncid.geospatial_lon_max = 180.
        ncid.processing_level = 'L1C'
    ncid.cdm_data_type = 'swath'
#    ncid.source_file = data.source_string
    ncid.gac_file = data.source_string
    ncid.lines_truncated_start = 1
    ncid.lines_truncated_end = data.lat.shape[0]
    if not himawari:
        ncid.orbital_temperature = data.orbital_temperature

    ncid.close()



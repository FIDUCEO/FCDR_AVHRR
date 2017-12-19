#!/usr/bin/env python
# works with xarray > 0.9.6
from __future__ import print_function
import xarray
# import netCDF4
import numpy as np 

file_in = 'FIDUCEO_FCDR_L1C_AVHRR_NOAA19_20100810005235_20100810023441_EASY_v0.3pre_fv1.1.1.nc'
# ncfile = netCDF4.Dataset(file_in,'r')
ncfile = xarray.open_dataset(file_in)
ncfile = ncfile.assign_coords(x=np.arange(ncfile.dims['x'])) 
ncfile = ncfile.assign_coords(y=np.arange(ncfile.dims['y']))
da = ncfile['Ch5_Bt']
idx = da.where(da==da.min(),drop=True)

print('x=',idx.x)
print('y=',idx.y)
ncfile.load()
print(ncfile.sel(x=idx.x,y=idx.y))







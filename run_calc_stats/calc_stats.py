#!/usr/bin/env python2.7

# call as: python2.7 calc_stats.py instrument, year, month, day, file_nc

# =======================================
# CALCULATE ROBUST STATISTICS
# =======================================
# Original code: J.Mittaz University of Reading 2017 
#
# Version 3.0
# 12 December, 2017
# michael.taylor@reading.ac.uk
# =======================================

import os
import os.path
import optparse
from  optparse import OptionParser
import sys
import numpy as np
import netCDF4
from netCDF4 import Dataset

class switch(object):
    value = None
    def __new__(class_, value):
        class_.value = value
        return True

def case(*args):
    return any((arg == switch.value for arg in args))

def main():

    instrument = sys.argv[1]
    year = int(sys.argv[2])
    month = int(sys.argv[3])
    day = int(sys.argv[4])
    file_nc = sys.argv[5]
    datestamp = file_nc[30:38]
    hhmms = file_nc[38:42]
    hhmme = file_nc[53:57]
    hh = int(file_nc[38:40])
    mm = int(file_nc[40:42])
    file_stem = instrument+"_"+str(datestamp)+"_"+str(hhmms)+"_"+str(hhmme)

    path = "/group_workspaces/cems2/fiduceo/Users/mtaylor/FCDR/make_fcdr_code/"
    path_in = path+instrument+"/"+str(year)+"/"+str("%02d" %month)+"/"+str("%02d" %day)+"/"
#    path_out = path+"stats/"+instrument+"/"+str(year)+"/"+str("%02d" %month)+"/"+str("%02d" %day)+"/"
    path_out = path_in

    file_in = path_in+file_nc         

    file_out_independent1 = path_out+file_stem+"_u_independent_Ch1.data"
    file_out_independent2 = path_out+file_stem+"_u_independent_Ch2.data"
    file_out_independent3a = path_out+file_stem+"_u_independent_Ch3a.data"
    file_out_independent3b = path_out+file_stem+"_u_independent_Ch3b.data"
    file_out_independent4 = path_out+file_stem+"_u_independent_Ch4.data"
    file_out_independent5 = path_out+file_stem+"_u_independent_Ch5.data"

    file_out_structured1 = path_out+file_stem+"_u_structured_Ch1.data"
    file_out_structured2 = path_out+file_stem+"_u_structured_Ch2.data"
    file_out_structured3a = path_out+file_stem+"_u_structured_Ch3a.data"
    file_out_structured3b = path_out+file_stem+"_u_structured_Ch3b.data"
    file_out_structured4 = path_out+file_stem+"_u_structured_Ch4.data"
    file_out_structured5 = path_out+file_stem+"_u_structured_Ch5.data"

# MT: 12-12-2017: add in channel data
    file_out_measurement1 = path_out+file_stem+"_measurement_Ch1.data"
    file_out_measurement2 = path_out+file_stem+"_measurement_Ch2.data"
    file_out_measurement3a = path_out+file_stem+"_measurement_Ch3a.data"
    file_out_measurement3b = path_out+file_stem+"_measurement_Ch3b.data"
    file_out_measurement4 = path_out+file_stem+"_measurement_Ch4.data"
    file_out_measurement5 = path_out+file_stem+"_measurement_Ch5.data"

# Remove old files if present
    if os.path.exists(file_out_independent1):
        os.remove(file_out_independent1)
    if os.path.exists(file_out_independent2):
        os.remove(file_out_independent2)
    if os.path.exists(file_out_independent3a):
        os.remove(file_out_independent3a)
    if os.path.exists(file_out_independent3b):
        os.remove(file_out_independent3b)
    if os.path.exists(file_out_independent4):
        os.remove(file_out_independent4)
    if os.path.exists(file_out_independent5):
        os.remove(file_out_independent5)

    if os.path.exists(file_out_structured1):
        os.remove(file_out_structured1)
    if os.path.exists(file_out_structured2):
        os.remove(file_out_structured2)
    if os.path.exists(file_out_structured3a):
        os.remove(file_out_structured3a)
    if os.path.exists(file_out_structured3b):
        os.remove(file_out_structured3b)
    if os.path.exists(file_out_structured4):
        os.remove(file_out_structured4)
    if os.path.exists(file_out_structured5):
        os.remove(file_out_structured5)

# MT: 12-12-2017: add in channel data
    if os.path.exists(file_out_measurement1):
        os.remove(file_out_measurement1)
    if os.path.exists(file_out_measurement2):
        os.remove(file_out_measurement2)
    if os.path.exists(file_out_measurement3a):
        os.remove(file_out_measurement3a)
    if os.path.exists(file_out_measurement3b):
        os.remove(file_out_measurement3b)
    if os.path.exists(file_out_measurement4):
        os.remove(file_out_measurement4)
    if os.path.exists(file_out_measurement5):
        os.remove(file_out_measurement5)

# Check to see if input file exists
    if not os.path.exists(file_in):
        return

# load variables from .nc output file    
    ncfile = netCDF4.Dataset(file_in,'r')

    u_independent_Ch1 = ncfile.variables['u_independent_Ch1'][:,:]
    gd = ~np.isfinite(u_independent_Ch1) | (u_independent_Ch1 < -1e20)
    if np.sum(gd) > 0:
        u_independent_Ch1[gd] = float('nan')
    u_independent_Ch2 = ncfile.variables['u_independent_Ch2'][:,:]
    gd = ~np.isfinite(u_independent_Ch2) | (u_independent_Ch2 < -1e20)
    if np.sum(gd) > 0:
        u_independent_Ch2[gd] = float('nan')
    u_independent_Ch3a = ncfile.variables['u_independent_Ch3a'][:,:]
    gd = ~np.isfinite(u_independent_Ch3a) | (u_independent_Ch3a < -1e20)
    if np.sum(gd) > 0:
        u_independent_Ch3a[gd] = float('nan')
    u_independent_Ch3b = ncfile.variables['u_independent_Ch3b'][:,:]
    gd = ~np.isfinite(u_independent_Ch3b) | (u_independent_Ch3b < -1e20)
    if np.sum(gd) > 0:
        u_independent_Ch3b[gd] = float('nan')
    u_independent_Ch4 = ncfile.variables['u_independent_Ch4'][:,:]
    gd = ~np.isfinite(u_independent_Ch4) | (u_independent_Ch4 < -1e20)
    if np.sum(gd) > 0:
        u_independent_Ch4[gd] = float('nan')
    u_independent_Ch5 = ncfile.variables['u_independent_Ch5'][:,:]
    gd = ~np.isfinite(u_independent_Ch5) | (u_independent_Ch5 < -1e20)
    if np.sum(gd) > 0:
        u_independent_Ch5[gd] = float('nan')

    u_structured_Ch1 = ncfile.variables['u_structured_Ch1'][:,:]
    gd = ~np.isfinite(u_structured_Ch1) | (u_structured_Ch1 < -1e20)
    if np.sum(gd) > 0:
        u_structured_Ch1[gd] = float('nan')
    u_structured_Ch2 = ncfile.variables['u_structured_Ch2'][:,:]
    gd = ~np.isfinite(u_structured_Ch2) | (u_structured_Ch2 < -1e20)
    if np.sum(gd) > 0:
        u_structured_Ch2[gd] = float('nan')
    u_structured_Ch3a = ncfile.variables['u_structured_Ch3a'][:,:]
    gd = ~np.isfinite(u_structured_Ch3a) | (u_structured_Ch3a < -1e20)
    if np.sum(gd) > 0:
        u_structured_Ch3a[gd] = float('nan')
    u_structured_Ch3b = ncfile.variables['u_structured_Ch3b'][:,:]
    gd = ~np.isfinite(u_structured_Ch3b) | (u_structured_Ch3b < -1e20)
    if np.sum(gd) > 0:
        u_structured_Ch3b[gd] = float('nan')
    u_structured_Ch4 = ncfile.variables['u_structured_Ch4'][:,:]
    gd = ~np.isfinite(u_structured_Ch4) | (u_structured_Ch4 < -1e20)
    if np.sum(gd) > 0:
        u_structured_Ch4[gd] = float('nan')
    u_structured_Ch5 = ncfile.variables['u_structured_Ch5'][:,:]
    gd = ~np.isfinite(u_structured_Ch5) | (u_structured_Ch5 < -1e20)
    if np.sum(gd) > 0:
        u_structured_Ch5[gd] = float('nan')

# MT: 12-12-2017: add in channel data
    measurement_Ch1 = ncfile.variables['Ch1_Ref'][:,:]
    gd = ~np.isfinite(measurement_Ch1) | (measurement_Ch1 < -1e20)
    if np.sum(gd) > 0:
        measurement_Ch1[gd] = float('nan')
    measurement_Ch2 = ncfile.variables['Ch2_Ref'][:,:]
    gd = ~np.isfinite(measurement_Ch2) | (measurement_Ch2 < -1e20)
    if np.sum(gd) > 0:
        measurement_Ch2[gd] = float('nan')
    measurement_Ch3a = ncfile.variables['Ch3a_Ref'][:,:]
    gd = ~np.isfinite(measurement_Ch3a) | (measurement_Ch3a < -1e20)
    if np.sum(gd) > 0:
        measurement_Ch3a[gd] = float('nan')
    measurement_Ch3b = ncfile.variables['Ch3b_Bt'][:,:]
    gd = ~np.isfinite(measurement_Ch3b) | (measurement_Ch3b < -1e20)
    if np.sum(gd) > 0:
        measurement_Ch3b[gd] = float('nan')
    measurement_Ch4 = ncfile.variables['Ch4_Bt'][:,:]
    gd = ~np.isfinite(measurement_Ch4) | (measurement_Ch4 < -1e20)
    if np.sum(gd) > 0:
        measurement_Ch4[gd] = float('nan')
    measurement_Ch5 = ncfile.variables['Ch5_Bt'][:,:]
    gd = ~np.isfinite(measurement_Ch5) | (measurement_Ch5 < -1e20)
    if np.sum(gd) > 0:
        measurement_Ch5[gd] = float('nan')
    ncfile.close()

    stat_array = np.zeros((18,13))
    stat_array[:,:] = -1e30

    for i in range(1,19):
        if 1 == i:
            i_var = u_independent_Ch1[:,:]
        elif 2 == i:
            i_var = u_independent_Ch2[:,:]
        elif 3 == i:
            i_var = u_independent_Ch3a[:,:]
        elif 4 == i:
            i_var = u_independent_Ch3b[:,:]
        elif 5 == i:
            i_var = u_independent_Ch4[:,:]
        elif 6 == i:
            i_var = u_independent_Ch5[:,:]

        elif 7 == i:
            i_var = u_structured_Ch1[:,:]
        elif 8 == i:
            i_var = u_structured_Ch2[:,:]
        elif 9 == i:
            i_var = u_structured_Ch3a[:,:]
        elif 10 == i:
            i_var = u_structured_Ch3b[:,:]
        elif 11 == i:
            i_var = u_structured_Ch4[:,:]
        elif 12 == i:
            i_var = u_structured_Ch5[:,:]

# MT: 12-12-2017: add in channel data
        elif 13 == i:
            i_var = measurement_Ch1[:,:]
        elif 14 == i:
            i_var = measurement_Ch2[:,:]
        elif 15 == i:
            i_var = measurement_Ch3a[:,:]
        elif 16 == i:
            i_var = measurement_Ch3b[:,:]
        elif 17 == i:
            i_var = measurement_Ch4[:,:]
        elif 18 == i:
            i_var = measurement_Ch5[:,:]

        stat_array[i-1,0]  = i                                   # variable
# MT: 03-12-2017: Calculation of fraction of good data
        gd = np.isfinite(i_var)
        bd = ~np.isfinite(i_var)
        n_gd = len(gd)
        n_bd = len(bd)
        n_i_var = n_gd + n_bd
        if np.sum(gd) > 1000:
#            stat_array[i-1,1]  = len(i_var[gd].flatten())       
            stat_array[i-1,1]  = n_gd                            # N
            stat_array[i-1,2]  = np.min(i_var[gd])               # min
            stat_array[i-1,3]  = np.max(i_var[gd])               # max
            stat_array[i-1,4]  = np.sum(i_var[gd])               # sum
            stat_array[i-1,5]  = np.sum(i_var[gd]**2)            # sum squared
            stat_array[i-1,6]  = np.mean(i_var[gd])              # mean
            stat_array[i-1,7]  = np.std(i_var[gd])               # st. dev.
            stat_array[i-1,8]  = np.var(i_var[gd])               # variance
            stat_array[i-1,9]  = np.percentile(i_var[gd], 25)    # Q1
            stat_array[i-1,10]  = np.percentile(i_var[gd], 50)   # Q2
# MT: 21-11-2017: Add number of NaN as a fraction of vector length
#            stat_array[i-1,10] = np.median(i_var[gd])           # median
            stat_array[i-1,11] = np.percentile(i_var[gd], 75)    # Q3
# MT: 21-11-2017: Add fraction of good data
            stat_array[i-1,12] = n_gd*1./n_i_var                 # N (good data fraction)

# Use numpy write to txt
    i=0
    with open(file_out_independent1,'w') as fp:
        fp.write('{0:04d} {1:02d} {2:02d} {3:02d} {4:02d} {5:2d} {6:e} {7:e} {8:e} {9:e} {10:e} {11:e} {12:e} {13:e} {14:e} {15:e} {16:e} {17:e}\n'.\
                     format(year,month,day,hh,mm,int(stat_array[i,0]),\
                                stat_array[i,1],\
                                stat_array[i,2],\
                                stat_array[i,3],\
                                stat_array[i,4],\
                                stat_array[i,5],\
                                stat_array[i,6],\
                                stat_array[i,7],\
                                stat_array[i,8],\
                                stat_array[i,9],\
                                stat_array[i,10],\
                                stat_array[i,11],\
                                stat_array[i,12]))

    i=1
    with open(file_out_independent2,'w') as fp:
        fp.write('{0:04d} {1:02d} {2:02d} {3:02d} {4:02d} {5:2d} {6:e} {7:e} {8:e} {9:e} {10:e} {11:e} {12:e} {13:e} {14:e} {15:e} {16:e} {17:e}\n'.\
                     format(year,month,day,hh,mm,int(stat_array[i,0]),\
                                stat_array[i,1],\
                                stat_array[i,2],\
                                stat_array[i,3],\
                                stat_array[i,4],\
                                stat_array[i,5],\
                                stat_array[i,6],\
                                stat_array[i,7],\
                                stat_array[i,8],\
                                stat_array[i,9],\
                                stat_array[i,10],\
                                stat_array[i,11],\
                                stat_array[i,12]))

    i=2
    with open(file_out_independent3a,'w') as fp:
        fp.write('{0:04d} {1:02d} {2:02d} {3:02d} {4:02d} {5:2d} {6:e} {7:e} {8:e} {9:e} {10:e} {11:e} {12:e} {13:e} {14:e} {15:e} {16:e} {17:e}\n'.\
                     format(year,month,day,hh,mm,int(stat_array[i,0]),\
                                stat_array[i,1],\
                                stat_array[i,2],\
                                stat_array[i,3],\
                                stat_array[i,4],\
                                stat_array[i,5],\
                                stat_array[i,6],\
                                stat_array[i,7],\
                                stat_array[i,8],\
                                stat_array[i,9],\
                                stat_array[i,10],\
                                stat_array[i,11],\
                                stat_array[i,12]))

    i=3
    with open(file_out_independent3b,'w') as fp:
        fp.write('{0:04d} {1:02d} {2:02d} {3:02d} {4:02d} {5:2d} {6:e} {7:e} {8:e} {9:e} {10:e} {11:e} {12:e} {13:e} {14:e} {15:e} {16:e} {17:e}\n'.\
                     format(year,month,day,hh,mm,int(stat_array[i,0]),\
                                stat_array[i,1],\
                                stat_array[i,2],\
                                stat_array[i,3],\
                                stat_array[i,4],\
                                stat_array[i,5],\
                                stat_array[i,6],\
                                stat_array[i,7],\
                                stat_array[i,8],\
                                stat_array[i,9],\
                                stat_array[i,10],\
                                stat_array[i,11],\
                                stat_array[i,12]))

    i=4
    with open(file_out_independent4,'w') as fp:
        fp.write('{0:04d} {1:02d} {2:02d} {3:02d} {4:02d} {5:2d} {6:e} {7:e} {8:e} {9:e} {10:e} {11:e} {12:e} {13:e} {14:e} {15:e} {16:e} {17:e}\n'.\
                     format(year,month,day,hh,mm,int(stat_array[i,0]),\
                                stat_array[i,1],\
                                stat_array[i,2],\
                                stat_array[i,3],\
                                stat_array[i,4],\
                                stat_array[i,5],\
                                stat_array[i,6],\
                                stat_array[i,7],\
                                stat_array[i,8],\
                                stat_array[i,9],\
                                stat_array[i,10],\
                                stat_array[i,11],\
                                stat_array[i,12]))

    i=5
    with open(file_out_independent5,'w') as fp:
        fp.write('{0:04d} {1:02d} {2:02d} {3:02d} {4:02d} {5:2d} {6:e} {7:e} {8:e} {9:e} {10:e} {11:e} {12:e} {13:e} {14:e} {15:e} {16:e} {17:e}\n'.\
                     format(year,month,day,hh,mm,int(stat_array[i,0]),\
                                stat_array[i,1],\
                                stat_array[i,2],\
                                stat_array[i,3],\
                                stat_array[i,4],\
                                stat_array[i,5],\
                                stat_array[i,6],\
                                stat_array[i,7],\
                                stat_array[i,8],\
                                stat_array[i,9],\
                                stat_array[i,10],\
                                stat_array[i,11],\
                                stat_array[i,12]))

    i=6
    with open(file_out_structured1,'w') as fp:
        fp.write('{0:04d} {1:02d} {2:02d} {3:02d} {4:02d} {5:2d} {6:e} {7:e} {8:e} {9:e} {10:e} {11:e} {12:e} {13:e} {14:e} {15:e} {16:e} {17:e}\n'.\
                     format(year,month,day,hh,mm,int(stat_array[i,0]),\
                                stat_array[i,1],\
                                stat_array[i,2],\
                                stat_array[i,3],\
                                stat_array[i,4],\
                                stat_array[i,5],\
                                stat_array[i,6],\
                                stat_array[i,7],\
                                stat_array[i,8],\
                                stat_array[i,9],\
                                stat_array[i,10],\
                                stat_array[i,11],\
                                stat_array[i,12]))

    i=7
    with open(file_out_structured2,'w') as fp:
        fp.write('{0:04d} {1:02d} {2:02d} {3:02d} {4:02d} {5:2d} {6:e} {7:e} {8:e} {9:e} {10:e} {11:e} {12:e} {13:e} {14:e} {15:e} {16:e} {17:e}\n'.\
                     format(year,month,day,hh,mm,int(stat_array[i,0]),\
                                stat_array[i,1],\
                                stat_array[i,2],\
                                stat_array[i,3],\
                                stat_array[i,4],\
                                stat_array[i,5],\
                                stat_array[i,6],\
                                stat_array[i,7],\
                                stat_array[i,8],\
                                stat_array[i,9],\
                                stat_array[i,10],\
                                stat_array[i,11],\
                                stat_array[i,12]))

    i=8
    with open(file_out_structured3a,'w') as fp:
        fp.write('{0:04d} {1:02d} {2:02d} {3:02d} {4:02d} {5:2d} {6:e} {7:e} {8:e} {9:e} {10:e} {11:e} {12:e} {13:e} {14:e} {15:e} {16:e} {17:e}\n'.\
                     format(year,month,day,hh,mm,int(stat_array[i,0]),\
                                stat_array[i,1],\
                                stat_array[i,2],\
                                stat_array[i,3],\
                                stat_array[i,4],\
                                stat_array[i,5],\
                                stat_array[i,6],\
                                stat_array[i,7],\
                                stat_array[i,8],\
                                stat_array[i,9],\
                                stat_array[i,10],\
                                stat_array[i,11],\
                                stat_array[i,12]))

    i=9
    with open(file_out_structured3b,'w') as fp:
        fp.write('{0:04d} {1:02d} {2:02d} {3:02d} {4:02d} {5:2d} {6:e} {7:e} {8:e} {9:e} {10:e} {11:e} {12:e} {13:e} {14:e} {15:e} {16:e} {17:e}\n'.\
                     format(year,month,day,hh,mm,int(stat_array[i,0]),\
                                stat_array[i,1],\
                                stat_array[i,2],\
                                stat_array[i,3],\
                                stat_array[i,4],\
                                stat_array[i,5],\
                                stat_array[i,6],\
                                stat_array[i,7],\
                                stat_array[i,8],\
                                stat_array[i,9],\
                                stat_array[i,10],\
                                stat_array[i,11],\
                                stat_array[i,12]))

    i=10
    with open(file_out_structured4,'w') as fp:
        fp.write('{0:04d} {1:02d} {2:02d} {3:02d} {4:02d} {5:2d} {6:e} {7:e} {8:e} {9:e} {10:e} {11:e} {12:e} {13:e} {14:e} {15:e} {16:e} {17:e}\n'.\
                     format(year,month,day,hh,mm,int(stat_array[i,0]),\
                                stat_array[i,1],\
                                stat_array[i,2],\
                                stat_array[i,3],\
                                stat_array[i,4],\
                                stat_array[i,5],\
                                stat_array[i,6],\
                                stat_array[i,7],\
                                stat_array[i,8],\
                                stat_array[i,9],\
                                stat_array[i,10],\
                                stat_array[i,11],\
                                stat_array[i,12]))

    i=11
    with open(file_out_structured5,'w') as fp:
        fp.write('{0:04d} {1:02d} {2:02d} {3:02d} {4:02d} {5:2d} {6:e} {7:e} {8:e} {9:e} {10:e} {11:e} {12:e} {13:e} {14:e} {15:e} {16:e} {17:e}\n'.\
                     format(year,month,day,hh,mm,int(stat_array[i,0]),\
                                stat_array[i,1],\
                                stat_array[i,2],\
                                stat_array[i,3],\
                                stat_array[i,4],\
                                stat_array[i,5],\
                                stat_array[i,6],\
                                stat_array[i,7],\
                                stat_array[i,8],\
                                stat_array[i,9],\
                                stat_array[i,10],\
                                stat_array[i,11],\
                                stat_array[i,12]))

    i=12
    with open(file_out_measurement1,'w') as fp:
        fp.write('{0:04d} {1:02d} {2:02d} {3:02d} {4:02d} {5:2d} {6:e} {7:e} {8:e} {9:e} {10:e} {11:e} {12:e} {13:e} {14:e} {15:e} {16:e} {17:e}\n'.\
                     format(year,month,day,hh,mm,int(stat_array[i,0]),\
                                stat_array[i,1],\
                                stat_array[i,2],\
                                stat_array[i,3],\
                                stat_array[i,4],\
                                stat_array[i,5],\
                                stat_array[i,6],\
                                stat_array[i,7],\
                                stat_array[i,8],\
                                stat_array[i,9],\
                                stat_array[i,10],\
                                stat_array[i,11],\
                                stat_array[i,12]))

    i=13
    with open(file_out_measurement2,'w') as fp:
        fp.write('{0:04d} {1:02d} {2:02d} {3:02d} {4:02d} {5:2d} {6:e} {7:e} {8:e} {9:e} {10:e} {11:e} {12:e} {13:e} {14:e} {15:e} {16:e} {17:e}\n'.\
                     format(year,month,day,hh,mm,int(stat_array[i,0]),\
                                stat_array[i,1],\
                                stat_array[i,2],\
                                stat_array[i,3],\
                                stat_array[i,4],\
                                stat_array[i,5],\
                                stat_array[i,6],\
                                stat_array[i,7],\
                                stat_array[i,8],\
                                stat_array[i,9],\
                                stat_array[i,10],\
                                stat_array[i,11],\
                                stat_array[i,12]))

    i=14
    with open(file_out_measurement3a,'w') as fp:
        fp.write('{0:04d} {1:02d} {2:02d} {3:02d} {4:02d} {5:2d} {6:e} {7:e} {8:e} {9:e} {10:e} {11:e} {12:e} {13:e} {14:e} {15:e} {16:e} {17:e}\n'.\
                     format(year,month,day,hh,mm,int(stat_array[i,0]),\
                                stat_array[i,1],\
                                stat_array[i,2],\
                                stat_array[i,3],\
                                stat_array[i,4],\
                                stat_array[i,5],\
                                stat_array[i,6],\
                                stat_array[i,7],\
                                stat_array[i,8],\
                                stat_array[i,9],\
                                stat_array[i,10],\
                                stat_array[i,11],\
                                stat_array[i,12]))

    i=15
    with open(file_out_measurement3b,'w') as fp:
        fp.write('{0:04d} {1:02d} {2:02d} {3:02d} {4:02d} {5:2d} {6:e} {7:e} {8:e} {9:e} {10:e} {11:e} {12:e} {13:e} {14:e} {15:e} {16:e} {17:e}\n'.\
                     format(year,month,day,hh,mm,int(stat_array[i,0]),\
                                stat_array[i,1],\
                                stat_array[i,2],\
                                stat_array[i,3],\
                                stat_array[i,4],\
                                stat_array[i,5],\
                                stat_array[i,6],\
                                stat_array[i,7],\
                                stat_array[i,8],\
                                stat_array[i,9],\
                                stat_array[i,10],\
                                stat_array[i,11],\
                                stat_array[i,12]))

    i=16
    with open(file_out_measurement4,'w') as fp:
        fp.write('{0:04d} {1:02d} {2:02d} {3:02d} {4:02d} {5:2d} {6:e} {7:e} {8:e} {9:e} {10:e} {11:e} {12:e} {13:e} {14:e} {15:e} {16:e} {17:e}\n'.\
                     format(year,month,day,hh,mm,int(stat_array[i,0]),\
                                stat_array[i,1],\
                                stat_array[i,2],\
                                stat_array[i,3],\
                                stat_array[i,4],\
                                stat_array[i,5],\
                                stat_array[i,6],\
                                stat_array[i,7],\
                                stat_array[i,8],\
                                stat_array[i,9],\
                                stat_array[i,10],\
                                stat_array[i,11],\
                                stat_array[i,12]))

    i=17
    with open(file_out_measurement5,'w') as fp:
        fp.write('{0:04d} {1:02d} {2:02d} {3:02d} {4:02d} {5:2d} {6:e} {7:e} {8:e} {9:e} {10:e} {11:e} {12:e} {13:e} {14:e} {15:e} {16:e} {17:e}\n'.\
                     format(year,month,day,hh,mm,int(stat_array[i,0]),\
                                stat_array[i,1],\
                                stat_array[i,2],\
                                stat_array[i,3],\
                                stat_array[i,4],\
                                stat_array[i,5],\
                                stat_array[i,6],\
                                stat_array[i,7],\
                                stat_array[i,8],\
                                stat_array[i,9],\
                                stat_array[i,10],\
                                stat_array[i,11],\
                                stat_array[i,12]))


if __name__ == "__main__":

    main()


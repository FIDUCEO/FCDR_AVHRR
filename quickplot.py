#!/usr/bin/env python2.7

# call as: python2.7 quickplot.py filename.nc varname

# =======================================
# Quickly plot a variable or channel data 
# =======================================
# Version 1.0
# 03 November, 2017
# michael.taylor AT reading DOT ac DOT uk
# =======================================

import os
import sys
from  optparse import OptionParser
import numpy as np
import netCDF4
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from cartopy import config
import cartopy.crs as ccrs
    
def plot_var(ncfile,filename,varname):

#    time = ncfile.variables['Time'][:]
#    lat = ncfile.variables['latitude'][:,:]
#    lon = ncfile.variables['longitude'][:,:]    
    val = ncfile.variables[varname][:,:]

    fig, ax = plt.subplots()
    plt.plot(val,'-',linewidth=0.2,marker='.',markersize=2,color='k')
    plt.title(varname)
    plt.savefig(str(varname+'.png'))

def plot_chx(ncfile,filename):

#    time = ncfile.variables['Time'][:]
#    lat = ncfile.variables['latitude'][:,:]
#    lon = ncfile.variables['longitude'][:,:]    

    try:
        ch1 = ncfile.variables['Ch1_Ref'][:,:]
        ch1_there = True
    except:
        ch1_there = False
    try:
        ch2 = ncfile.variables['Ch2_Ref'][:,:]
        ch2_there = True
    except:
        ch2_there = False
    try:
        ch3a = ncfile.variables['Ch3a_Ref'][:,:]
        ch3a_there = True
    except:
        ch3a_there = False
    try:
        ch3b = ncfile.variables['Ch3b_Bt'][:,:]
        ch3b_there = True
    except:
        ch3b_there = False
    try:
        ch4 = ncfile.variables['Ch4_Bt'][:,:]
        ch4_there = True
    except:
        ch4_there = False
    try:
        ch5 = ncfile.variables['Ch5_Bt'][:,:]
        ch5_there = True
    except:
        ch5_there = False

    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(6, 6), sharex=True, sharey=False)
    if ch1_there:
        axes[0, 0].plot(ch1,'-',linewidth=0.2,marker='.',markersize=2,color='k')
    if ch2_there:
        axes[0, 1].plot(ch2,'-',linewidth=0.2,marker='.',markersize=2,color='k')
    if ch3a_there:
        axes[1, 0].plot(ch3a,'-',linewidth=0.2,marker='.',markersize=2,color='k')
    if ch3b_there:
        axes[1, 1].plot(ch3b,'-',linewidth=0.2,marker='.',markersize=2,color='k')
    if ch4_there:
        axes[2, 0].plot(ch4,'-',linewidth=0.2,marker='.',markersize=2,color='k')
    if ch5_there:
        axes[2, 1].plot(ch5,'-',linewidth=0.2,marker='.',markersize=2,color='k')
    axes[0, 0].set_title('Ch1')
    axes[0, 1].set_title('Ch2')
    axes[1, 0].set_title('Ch3a')
    axes[1, 1].set_title('Ch3b')
    axes[2, 0].set_title('Ch4')
    axes[2, 1].set_title('Ch5')

    fig.subplots_adjust(hspace=0.4)
    #plt.show()
    plt.savefig('chx.png')

def plot_chx_independent(ncfile,filename):

#    time = ncfile.variables['Time'][:]
#    lat = ncfile.variables['latitude'][:,:]
#    lon = ncfile.variables['longitude'][:,:]    

    try:
        u_independent_Ch1 = ncfile.variables['u_independent_Ch1'][:,:]
        u_independent_Ch1_there = True
    except:
        u_independent_Ch1_there = False
    try:
        u_independent_Ch2 = ncfile.variables['u_independent_Ch2'][:,:]
        u_independent_Ch2_there = True
    except:
        u_independent_Ch2_there = False
    try:
        u_independent_Ch3a = ncfile.variables['u_independent_Ch3a'][:,:]
        u_independent_Ch3a_there = True
    except:
        u_independent_Ch3a_there = False
    try:
        u_independent_Ch3b = ncfile.variables['u_independent_Ch3b'][:,:]
        u_independent_Ch3b_there = True
    except:
        u_independent_Ch3b_there = False
    try:
        u_independent_Ch4 = ncfile.variables['u_independent_Ch4'][:,:]
        u_independent_Ch4_there = True
    except:
        u_independent_Ch4_there = False
    try:
        u_independent_Ch5 = ncfile.variables['u_independent_Ch5'][:,:]
        u_independent_Ch5_there = True
    except:
        u_independent_Ch5_there = False

    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(6, 6), sharex=True, sharey=False)
    if u_independent_Ch1_there:
        axes[0, 0].plot(u_independent_Ch1,'-',linewidth=0.2,marker='.',markersize=2,color='k')
    if u_independent_Ch2_there:
        axes[0, 1].plot(u_independent_Ch2,'-',linewidth=0.2,marker='.',markersize=2,color='k')
    if u_independent_Ch3a_there:
        axes[1, 0].plot(u_independent_Ch3a,'-',linewidth=0.2,marker='.',markersize=2,color='k')
    if u_independent_Ch3b_there:
        axes[1, 1].plot(u_independent_Ch3b,'-',linewidth=0.2,marker='.',markersize=2,color='k')
    if u_independent_Ch4_there:
        axes[2, 0].plot(u_independent_Ch4,'-',linewidth=0.2,marker='.',markersize=2,color='k')
    if u_independent_Ch5_there:
        axes[2, 1].plot(u_independent_Ch5,'-',linewidth=0.2,marker='.',markersize=2,color='k')
    axes[0, 0].set_title('Ch1_independent')
    axes[0, 1].set_title('Ch2_independent')
    axes[1, 0].set_title('Ch3a_independent')
    axes[1, 1].set_title('Ch3b_independent')
    axes[2, 0].set_title('Ch4_independent')
    axes[2, 1].set_title('Ch5_independent')

    fig.subplots_adjust(hspace=0.4)
    #plt.show()
    plt.savefig('chx_independent.png')

def plot_chx_structured(ncfile,filename):

#    time = ncfile.variables['Time'][:]
#    lat = ncfile.variables['latitude'][:,:]
#    lon = ncfile.variables['longitude'][:,:]    

    try:
        u_structured_Ch1 = ncfile.variables['u_structured_Ch1'][:,:]
        u_structured_Ch1_there = True
    except:
        u_structured_Ch1_there = False
    try:
        u_structured_Ch2 = ncfile.variables['u_structured_Ch2'][:,:]
        u_structured_Ch2_there = True
    except:
        u_structured_Ch2_there = False
    try:
        u_structured_Ch3a = ncfile.variables['u_structured_Ch3a'][:,:]
        u_structured_Ch3a_there = True
    except:
        u_structured_Ch3a_there = False
    try:
        u_structured_Ch3b = ncfile.variables['u_structured_Ch3b'][:,:]
        u_structured_Ch3b_there = True
    except:
        u_structured_Ch3b_there = False
    try:
        u_structured_Ch4 = ncfile.variables['u_structured_Ch4'][:,:]
        u_structured_Ch4_there = True
    except:
        u_structured_Ch4_there = False
    try:
        u_structured_Ch5 = ncfile.variables['u_structured_Ch5'][:,:]
        u_structured_Ch5_there = True
    except:
        u_structured_Ch5_there = False

    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(6, 6), sharex=True, sharey=False)
    if u_structured_Ch1_there:
        axes[0, 0].plot(u_structured_Ch1,'-',linewidth=0.2,marker='.',markersize=2,color='k')
    if u_structured_Ch2_there:
        axes[0, 1].plot(u_structured_Ch2,'-',linewidth=0.2,marker='.',markersize=2,color='k')
    if u_structured_Ch3a_there:
        axes[1, 0].plot(u_structured_Ch3a,'-',linewidth=0.2,marker='.',markersize=2,color='k')
    if u_structured_Ch3b_there:
        axes[1, 1].plot(u_structured_Ch3b,'-',linewidth=0.2,marker='.',markersize=2,color='k')
    if u_structured_Ch4_there:
        axes[2, 0].plot(u_structured_Ch4,'-',linewidth=0.2,marker='.',markersize=2,color='k')
    if u_structured_Ch5_there:
        axes[2, 1].plot(u_structured_Ch5,'-',linewidth=0.2,marker='.',markersize=2,color='k')
    axes[0, 0].set_title('Ch1_structured')
    axes[0, 1].set_title('Ch2_structured')
    axes[1, 0].set_title('Ch3a_structured')
    axes[1, 1].set_title('Ch3b_structured')
    axes[2, 0].set_title('Ch4_structured')
    axes[2, 1].set_title('Ch5_structured')

    fig.subplots_adjust(hspace=0.4)
    #plt.show()
    plt.savefig('chx_structured.png')
    
if __name__ == "__main__":

    parser = OptionParser('usage: %filename varname')
    (options, args) = parser.parse_args()
    filename = args[0]
    ncfile = Dataset(filename,'r')
    nargs = len(args)

    if nargs !=1 and nargs !=2:
        parser.error('incorrect number of input arguments')

    if nargs == 1:
        plot_chx(ncfile,filename)
        plot_chx_independent(ncfile,filename)
        plot_chx_structured(ncfile,filename)
    elif nargs == 2:
        varname = args[1]
        plot_var(ncfile,filename,varname)   





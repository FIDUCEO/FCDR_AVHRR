#!/usr/bin/env python2.7                                                                    

# * Copyright (C) 2017 M.Taylor University of Reading                                       
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

# call as: python2.7 quickplot_dat.py filename

# ==================================================
# PLOT 5-point stats for each orbit in sensor series
# ==================================================
# Version 0.1
# 14 December, 2017
# michael.taylor AT reading DOT ac DOT uk
# ==================================================

from  optparse import OptionParser 
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

def plot_stats(dataset,filename):

#   'year', 'month', 'day', 'hour', 'minute', [00-04]
#   'channel','ndata','minval','maxval','x',  [05-09]
#   'x2','mean','stdev','var','q_25',         [10-14]
#   'median','q_75','nan_fraction')           [15-17]

    v_x = []
    v_min = []
    v_q1 = []
    v_q2 = []
    v_q3 = []
    v_max = []
    for line in dataset:      
        p = line.split()
        v_x.append(float(p[9]))
        v_min.append(float(p[7]))
        v_q1.append(float(p[14]))
        v_q2.append(float(p[15]))
        v_q3.append(float(p[16]))
        v_max.append(float(p[8]))
        
    x = range(0,len(v_x))
    cutoff = -10000

    y_x = np.array(v_x,dtype=float)
    bd = (y_x < cutoff)
    y_x[bd] = np.nan

    y_min = np.array(v_min,dtype=float)
    bd = (y_min < cutoff)
    y_min[bd] = np.nan

    y_q1 = np.array(v_q1,dtype=float)
    bd = (y_q1 < cutoff)
    y_q1[bd] = np.nan

    y_q2 = np.array(v_q2,dtype=float)
    bd = (y_q2 < cutoff)
    y_q2[bd] = np.nan

    y_q3 = np.array(v_q3,dtype=float)
    bd = (y_q3 < cutoff)
    y_q3[bd] = np.nan

    y_max = np.array(v_max,dtype=float)
    bd = (y_max < cutoff)    
    y_max[bd] = np.nan

    colormap = plt.cm.gist_ncar 
#    labels = ['max','q3','q2','q1','min']
    labels = ['max','q2','min']
   
    fig, axes = plt.subplots()
    plt.plot(x,y_max)
#    plt.plot(x,y_q3)
    plt.plot(x,y_q2)
#    plt.plot(x,y_q1)
    plt.plot(x,y_min)
    plt.xlabel('Orbit')
    plt.ylabel('Value')
    plt.legend(labels, ncol=4, loc='upper center', columnspacing=1.0, labelspacing=0.0, handletextpad=0.0, handlelength=1.5, fancybox=True, shadow=True)
    plt.title(filename)
    plt.savefig('stats_'+filename+'.png')

if __name__ == "__main__":

    parser = OptionParser('usage: %filename')
    (options, args) = parser.parse_args()
    filename = args[0]

    f = open(filename, 'r') 
    dataset = f.readlines() 
    f.close() 

    plot_stats(dataset,filename)  


#!/usr/bin/env python2.7

# call as: python2.7 quickplot_ascii.py filename, column

# =======================================
# PLOT column of an array
# =======================================
# Version 1.0
# 14 December, 2017
# michael.taylor AT reading DOT ac DOT uk
# =======================================

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
#    plt.legend(labels, ncol=4, loc='upper center', bbox_to_anchor=[0.5, 1.1], columnspacing=1.0, labelspacing=0.0, handletextpad=0.0, handlelength=1.5, fancybox=True, shadow=True)
    plt.title(filename)
#    plt.show()
    plt.savefig('stats_'+filename+'.png')

if __name__ == "__main__":

#    parser = OptionParser('usage: %filename column')
    parser = OptionParser('usage: %filename')
    (options, args) = parser.parse_args()
    filename = args[0]
#    column = int(args[1])

    f = open(filename, 'r') 
    dataset = f.readlines() 
    f.close() 

    plot_stats(dataset,filename)  

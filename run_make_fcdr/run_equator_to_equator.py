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


import numpy as np
import datetime 
import calendar
import os
from  optparse import OptionParser
import subprocess

def run_equator_to_equator(name):

    if name == 'NOAA06':
        dirname = 'AVHRR06_G'
    elif name == 'NOAA07':
        dirname = 'AVHRR07_G'
    elif name == 'NOAA08':
        dirname = 'AVHRR08_G'
    elif name == 'NOAA09':
        dirname = 'AVHRR09_G'
    elif name == 'NOAA10':
        dirname = 'AVHRR10_G'
    elif name == 'NOAA11':
        dirname = 'AVHRR11_G'
    elif name == 'NOAA12':
        dirname = 'AVHRR12_G'
    elif name == 'NOAA14':
        dirname = 'AVHRR14_G'
    elif name == 'NOAA15':
        dirname = 'AVHRR15_G'
    elif name == 'NOAA16':
        dirname = 'AVHRR16_G'
    elif name == 'NOAA17':
        dirname = 'AVHRR17_G'
    elif name == 'NOAA18':
        dirname = 'AVHRR18_G'
    elif name == 'NOAA19':
        dirname = 'AVHRR19_G'
    elif name == 'METOPA':
        dirname = 'AVHRRMTA_G'

    for year in range(1978,2017):
        for month in range(1,13):
            maxday = calendar.monthrange(year,month)[1]
            for day in range(1,maxday+1):
                directory = '/group_workspaces/cems2/esacci_sst/input/avhrr/l1b/{0}/v1/{1:04d}/{2:02d}/{3:02d}'.format(dirname,year,month,day)
                if os.path.isdir(directory):
                    command='python2.7 equator_to_equator.py {0} {1:04d} {2:02d} {3:02d} Y'.format(name,year,month,day)
                    subprocess.call(command,shell=True)

if __name__ == "__main__":
    
    parser = OptionParser("usage: %prog instr_name")
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")

    run_equator_to_equator(args[0])


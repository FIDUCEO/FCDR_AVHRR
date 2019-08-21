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
import os.path
import sys
import uuid
import glob
import stat
from  optparse import OptionParser
import subprocess

def __get_avhrr_dirname(name):

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

    return dirname

def make_shell_command(dirname,year,month,day,file_nc,i,directory):
         
    currentdir = os.getcwd()
    outdir = '{0}/{1}/{2:04d}/{3:02d}/{4:02d}'.\
        format(currentdir,dirname,year,month,day)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    os.chdir(outdir)
    try:
        os.symlink('/group_workspaces/cems2/fiduceo/Users/mtaylor/FCDR/make_fcdr_co\
de/calc_stats.py','calc_stats.py')
    except:
        pass

    # job submission script file 
    job_file = 'run.{0:06d}.sh'.format(i)
    job_log = 'run.{0:06d}.log'.format(i)
    with open(job_file,'w') as fp:
        job_stem = []
        base = os.path.basename(file_nc)
        str_tuple = os.path.splitext(base)
        job_stem = str_tuple[0]+'.'+str(uuid.uuid4())
        job_str = 'python2.7 calc_stats.py {0} {1:04d} {2:02d} {3:02d} {4}\n'.format(dirname,year,month,day,file_nc,job_stem)
        fp.write(job_str)    

    os.chmod(job_file,stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)
    job_name='./'+job_file
    job = ['bsub','-q', 'short-serial','-W', '01:00','-oo', job_log, job_name]
    subprocess.call(job)
    os.chdir(currentdir)                    

def run_calc_stats(dirname):

    for year in range(1978,2017):
        for month in range(1,13):
            maxday = calendar.monthrange(year,month)[1]
            for day in range(1,maxday+1):
                directory = '/group_workspaces/cems2/fiduceo/Users/mtaylor/FCDR/make_fcdr_code/{0}/{1:04d}/{2:02d}/{3:02d}'.format(dirname,year,month,day)
                if os.path.isdir(directory):       
#                    filelist = os.listdir(directory)
                    nclist = os.path.join(directory,'*.nc')
                    filelist = glob.glob(nclist)
                    for i in range(len(filelist)):
                        file_nc_all = str(filelist[i])
                        file_nc = os.path.basename(file_nc_all)
                        make_shell_command(dirname,year,month,day,file_nc,i,directory)       
 
if __name__ == "__main__":
    
    parser = OptionParser("usage: %prog instr_name")
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")

    dirname = __get_avhrr_dirname(args[0])
    run_calc_stats(dirname)




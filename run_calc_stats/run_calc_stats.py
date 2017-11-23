import numpy as np
import datetime 
import calendar
import os
import os.path
import sys
from  optparse import OptionParser
import subprocess

def run_calc_stats(name):

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
                directory = '/group_workspaces/cems2/fiduceo/Users/mtaylor/FCDR/make_fcdr_code/{0}/{1:04d}/{2:02d}/{3:02d}'.format(dirname,year,month,day)
                if os.path.isdir(directory):
                    directory_stats = '/group_workspaces/cems2/fiduceo/Users/mtaylor/FCDR/make_fcdr_code/stats/{0}/{1:04d}/{2:02d}/{3:02d}'.format(dirname,year,month,day)
                    if not os.path.isdir(directory_stats):
                        os.makedirs(directory_stats)

                        filelist = os.listdir(directory)
                        print filelist
                        for i in range(len(filelist)):
                            file_nc = filelist[i]
                            command='python2.7 calc_stats.py {0} {1:04d} {2:02d} {3:02d} {4}'.format(name,year,month,day,file_nc)
                            subprocess.call(command,shell=True)

if __name__ == "__main__":
    
    parser = OptionParser("usage: %prog instr_name")
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")

    run_calc_stats(args[0])


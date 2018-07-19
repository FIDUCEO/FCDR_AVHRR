# * Copyright (C) 2017 J.Mittaz University of Reading
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
# Code to try and go equator to equator for the AVHRR
# and generate the FIDUCEO FCDR (pre-beta easy)
#
# Written by J.Mittaz UoR
# --------------------------------
# Modified M.Taylor UoR: 
# MT: 24-10-2017: import local ECTs to implement additional logical processing in orbit merging
# MT: 24-10-2017: add logic for case where end of orbit overlaps start of ECT window
# MT: 24-10-2017: add logic for case where start of orbit overlaps end of ECT window
# MT: 24-10-2017: add more logic for case where end of orbit overlaps start of ECT window
# MT: 15-11-2017: convert if to elif to fix repetition of orbit files in filelist



import numpy as np
import ephem
import datetime
import os
import glob
import uuid
import stat
import subprocess
import gzip
import bz2
from  optparse import OptionParser

# Get AVHRR type from filename
def get_avhrr_dir_name(avhrr_name):

    if 'TIROSN' == avhrr_name:
        avhrr_dir_name = 'AVHRRTN_G'
    elif 'NOAA06' == avhrr_name:
        avhrr_dir_name = 'AVHRR06_G'
    elif 'NOAA07' == avhrr_name:
        avhrr_dir_name = 'AVHRR07_G'
    elif 'NOAA08' == avhrr_name:
        avhrr_dir_name = 'AVHRR08_G'
    elif 'NOAA09' == avhrr_name:
        avhrr_dir_name = 'AVHRR09_G'
    elif 'NOAA10' == avhrr_name:
        avhrr_dir_name = 'AVHRR10_G'
    elif 'NOAA11' == avhrr_name:
        avhrr_dir_name = 'AVHRR11_G'
    elif 'NOAA12' == avhrr_name:
        avhrr_dir_name = 'AVHRR12_G'
    elif 'NOAA14' == avhrr_name:
        avhrr_dir_name = 'AVHRR14_G'
    elif 'NOAA15' == avhrr_name:
        avhrr_dir_name = 'AVHRR15_G'
    elif 'NOAA16' == avhrr_name:
        avhrr_dir_name = 'AVHRR16_G'
    elif 'NOAA17' == avhrr_name:
        avhrr_dir_name = 'AVHRR17_G'
    elif 'NOAA18' == avhrr_name:
        avhrr_dir_name = 'AVHRR18_G'
    elif 'NOAA19' == avhrr_name:
        avhrr_dir_name = 'AVHRR19_G'
    elif 'METOPA' == avhrr_name:
        avhrr_dir_name = 'AVHRRMTA_G'
    elif 'METOPB' == avhrr_name:
        avhrr_dir_name = 'AVHRRMTB_G'
    else:
        print avhrr_name
        raise Exception,"Cannot find avhrr_name"
    return avhrr_dir_name

# Get information from filename
class get_file_data(object):

    # Get AVHRR type from filename
    def get_avhrr_type(self,infilename):

        filesp = os.path.split(infilename)
        filename = filesp[1]

        # Extract instrument part of filename
        avhrr_name=filename[9:11]

        if 'TN' == avhrr_name:
            self.avhrr_dir_name = 'AVHRRTN_G'
            self.instr = 'TIROSN'
        elif 'NA' == avhrr_name:
            self.avhrr_dir_name = 'AVHRR06_G'
            self.instr = 'NOAA06'
        elif 'NC' == avhrr_name:
            self.avhrr_dir_name = 'AVHRR07_G'
            self.instr = 'NOAA07'
        elif 'NE' == avhrr_name:
            self.avhrr_dir_name = 'AVHRR08_G'
            self.instr = 'NOAA08'
        elif 'NF' == avhrr_name:
            self.avhrr_dir_name = 'AVHRR09_G'
            self.instr = 'NOAA09'
        elif 'NG' == avhrr_name:
            self.avhrr_dir_name = 'AVHRR10_G'
            self.instr = 'NOAA10'
        elif 'NH' == avhrr_name:
            self.avhrr_dir_name = 'AVHRR11_G'
            self.instr = 'NOAA11'
        elif 'ND' == avhrr_name:
            self.avhrr_dir_name = 'AVHRR12_G'
            self.instr = 'NOAA12'
        elif 'NJ' == avhrr_name:
            self.avhrr_dir_name = 'AVHRR14_G'
            self.instr = 'NOAA14'
        elif 'NK' == avhrr_name:
            self.avhrr_dir_name = 'AVHRR15_G'
            self.instr = 'NOAA15'
        elif 'NL' == avhrr_name:
            self.avhrr_dir_name = 'AVHRR16_G'
            self.instr = 'NOAA16'
        elif 'NM' == avhrr_name:
            self.avhrr_dir_name = 'AVHRR17_G'
            self.instr = 'NOAA17'
        elif 'NN' == avhrr_name:
            self.avhrr_dir_name = 'AVHRR18_G'
            self.instr = 'NOAA18'
        elif 'NP' == avhrr_name:
            self.avhrr_dir_name = 'AVHRR19_G'
            self.instr = 'NOAA19'
        elif 'M2' == avhrr_name:
            self.avhrr_dir_name = 'AVHRRMTA_G'
            self.instr = 'METOPA'
        elif 'M1' == avhrr_name:
            self.avhrr_dir_name = 'AVHRRMTB_G'
            self.instr = 'METOPB'
        else:
            print 'get_avhrr_type'
            print 'avhrr_name = ',avhrr_name
            raise Exception,"Cannot find avhrr_name"

    # Get start/end date from filename
    def get_avhrr_time(self,filename_in):

        file_tuple = os.path.split(filename_in)
        filename = file_tuple[1]

        # Get year
        year = int(filename[13:15])
        if year < 50:
            year = 2000+year
        else:
            year = 1900+year

        # Day of year
        dayno = int(filename[15:18])

        # Start time
        start_hour = int(filename[20:22])
        start_minute = int(filename[22:24])

        # End time
        end_hour = int(filename[26:28])
        end_minute = int(filename[28:30])

        # Get start time in datetime format
        start_string = '{0:04d}-{1:03d}:{2:02d}:{3:02d}'.\
            format(year,dayno,start_hour,start_minute)
        self.start_time = datetime.datetime.strptime(start_string,\
                                                         "%Y-%j:%H:%M")

        # Get end time
        # Have to deal with crossing day boundary
        start_time = start_hour + (start_minute/60.)
        end_time = end_hour + (end_minute/60.)

        # If end_time < start_time then crossed day boundary
        if end_time < start_time:
            # Get start time (year/day)
            start_date_string = '{0:04d}-{1:03d}'.format(year,dayno)
            start_date = datetime.datetime.strptime(start_date_string,"%Y-%j")
            # Add a day to datetime
            end_date = start_date + datetime.timedelta(1)
            yday = end_date.timetuple().tm_yday
            end_string = '{0:04d}-{1:03d}:{2:02d}:{3:02d}'.\
                format(end_date.year,yday,\
                           end_hour,end_minute)
        else:
            # If within same day just get string
            end_string = '{0:04d}-{1:03d}:{2:02d}:{3:02d}'.\
                format(year,dayno,end_hour,end_minute)

        # End time datetime
        self.end_time = datetime.datetime.strptime(end_string,"%Y-%j:%H:%M")

        # Get middle time - have to use timedelta for arithmetic
        self.time = self.start_time + (self.end_time - self.start_time)/2

    # Get information from filename
    def get_file_information(self,filename):

        # filemame is of structure
        #
        # NSS.GHRR.NL.D06001.S0000.E0043.B2719899.GC.g
        #          
        # Get AVHRR type
        self.get_avhrr_type(filename)
        
        # Get start/end time
        self.get_avhrr_time(filename)

    # Init for class
    def __init__(self,filename):
        self.get_file_information(filename)

# Deal with blacklist data
class blacklist(object):
    # Read in one day of blacklist
    def get_blacklist_day(self,instrument,year,month,day):
        # Get filename of blacklist
        mtaylor_path = "/group_workspaces/cems2/fiduceo/Users/mtaylor/avhrr_l1b/"
        dir_in = '{0}/listes_orbites/liste_orbites_day/{1}/{2:04d}/{3:02d}/'.\
            format(mtaylor_path,instrument,year,month)
        file_in = '{0}/blackliste_{1}_{2:04d}_{3:02d}_{4:02d}.data'.\
            format(dir_in,instrument,year,month,day)

        # Read blacklist file
        if os.path.isfile(file_in)==True:
            try:
                blacklist,redundant,reason= np.loadtxt(file_in_2,unpack=True,dtype=np.str,usecols=(0,2,3))
            except:
                blacklist=['None']
                redundant = -1
                reason = -1
        else:
            blacklist=['None']
            redundant = -1
            reason = -1

        return blacklist,redundant,reason

    # Read in blacklist based on a time
    # add_day to get blacklist over two days if needed
    def get_blacklist(self,instrument,year,month,day,add_day):

        blacklist,redundant,reason = \
            self.get_blacklist_day(instrument,year,month,day)

        if redundant != -1 and reason != -1:
            # Move on a day
            if add_day:
                date = datetime.datetime(year,month,day) + \
                    datetime.timedelta(1)
                blacklist2,redundant2,reason2 = \
                    self.get_blacklist_day(instrument,date.year,date.month,\
                                               date.day)
                if redundant2 != -1 and reason2 != -1:
                    blacklist = np.append(blacklist,blacklist2)
                    redundant = np.append(redundant,redundant2)
                    reason = np.append(reason,reason2)
        self.blacklist = blacklist
        self.redundant = redundant
        self.reason = reason
                
    # Check file against blacklist
    def check_blacklist(self,filename):

        bad_orbit = False
        # Check is we have an array or just a string
        if isinstance(self.blacklist,np.ndarray):
            # is the file in blacklist?
            if any(x == filename for x in self.blacklist):
                b = np.where(self.blacklist==filename)
                if self.reason[b]=="too_small" or \
                        self.reason[b]=="too_long" or \
                        self.reason[b]=='bad_l1c_quality' or \
                        self.reason[b]=="ground_station_duplicate" or \
                        self.reason[b]=='along_track_too_long' or \
                        self.redundant[b]=='1':
                    bad_orbit = True
        else:
            # Just a single string
            if self.blacklist == filename:
                if self.reason=="too_small" or \
                        self.reason=="too_long" or \
                        self.reason=='bad_l1c_quality' or \
                        self.reason=="ground_station_duplicate" or \
                        self.reason=='along_track_too_long' or \
                        self.redundant=='1':
                    bad_orbit = True

        # return value
        return bad_orbit

    def __init__(self,instrument,year,month,day,add_day=False):

        self.get_blacklist(instrument,year,month,day,add_day)

# Get equator crossing times
class tle_data(object):

    def read_tle(self,instr,year,month,day):

        time = datetime.datetime(year,month,day)
        self.name = instr
        if 'NOAA06' == instr:
            filename = '/group_workspaces/cems2/fiduceo/Users/mtaylor/FCDR/make_fcdr_code/merge_code/TLE/noaa-06.txt'
        elif 'NOAA07' == instr:
            filename = '/group_workspaces/cems2/fiduceo/Users/mtaylor/FCDR/make_fcdr_code/merge_code/TLE/noaa-07.txt'
        elif 'NOAA08' == instr:
            filename = '/group_workspaces/cems2/fiduceo/Users/mtaylor/FCDR/make_fcdr_code/merge_code/TLE/noaa-08.txt'
        elif 'NOAA09' == instr:
            filename = '/group_workspaces/cems2/fiduceo/Users/mtaylor/FCDR/make_fcdr_code/merge_code/TLE/noaa-09.txt'
        elif 'NOAA10' == instr:
            filename = '/group_workspaces/cems2/fiduceo/Users/mtaylor/FCDR/make_fcdr_code/merge_code/TLE/noaa-10.txt'
        elif 'NOAA11' == instr:
            filename = '/group_workspaces/cems2/fiduceo/Users/mtaylor/FCDR/make_fcdr_code/merge_code/TLE/noaa-11.txt'
        elif 'NOAA12' == instr:
            filename = '/group_workspaces/cems2/fiduceo/Users/mtaylor/FCDR/make_fcdr_code/merge_code/TLE/noaa-12.txt'
        elif 'NOAA14' == instr:
            filename = '/group_workspaces/cems2/fiduceo/Users/mtaylor/FCDR/make_fcdr_code/merge_code/TLE/noaa-14.txt'
        elif 'NOAA15' == instr:
            filename = '/group_workspaces/cems2/fiduceo/Users/mtaylor/FCDR/make_fcdr_code/merge_code/TLE/noaa-15.txt'
        elif 'NOAA16' == instr:
            filename = '/group_workspaces/cems2/fiduceo/Users/mtaylor/FCDR/make_fcdr_code/merge_code/TLE/noaa-16.txt'
        elif 'NOAA17' == instr:
            filename = '/group_workspaces/cems2/fiduceo/Users/mtaylor/FCDR/make_fcdr_code/merge_code/TLE/noaa-17.txt'
        elif 'NOAA18' == instr:
            filename = '/group_workspaces/cems2/fiduceo/Users/mtaylor/FCDR/make_fcdr_code/merge_code/TLE/noaa-18.txt'
        elif 'NOAA19' == instr:
            filename = '/group_workspaces/cems2/fiduceo/Users/mtaylor/FCDR/make_fcdr_code/merge_code/TLE/noaa-19.txt'
        elif 'METOPA' == instr:
            filename = '/group_workspaces/cems2/fiduceo/Users/mtaylor/FCDR/make_fcdr_code/merge_code/TLE/metop-a.txt'
        elif 'METOPB' == instr:
            filename = '/group_workspaces/cems2/fiduceo/Users/mtaylor/FCDR/make_fcdr_code/merge_code/TLE/metop-b.txt'
        else:
            print 'instr = ',instr
            raise Exception,"Cannot find instr"

        if not os.path.exists(filename):
            print 'instr = ',instr
            raise Exception,"TLE file not found : "+filename

        # Read file and get closest pair
        min_time = 1e30
        with open(filename,'r') as fp:
            lines = fp.readlines()
            k=0            
            while k < len(lines):
                line1 = lines[k].strip()
                line2 = lines[k+1].strip()
                year = int(line1[17:20])
                daynum = float(line1[20:32])
                day = int(daynum)
                temp = (daynum-day)*24.
                hour = int(temp)
                temp = (temp-hour)*60.
                minute = int(temp)
                temp = (temp-minute)*60.
                second = int(temp)
                timestr = '{0:02d}-{1:03d}:{2:02d}:{3:02d}:{4:02d}'.\
                    format(year,day,hour,minute,second)
                tletime = datetime.datetime.strptime(timestr,"%y-%j:%H:%M:%S")
                difftime = (time-tletime)
                if abs(difftime.total_seconds()) < min_time:
                    stored_line1 = line1
                    stored_line2 = line2
                    min_time = difftime.total_seconds()
                k=k+2

        return stored_line1,stored_line2

    def get_ascending_descending_type(self):

        if 'TIROSN' == self.name:
            ascending=True
        elif 'NOAA06' == self.name:
            ascending=False
        elif 'NOAA07' == self.name:
            ascending=True
        elif 'NOAA08' == self.name:
            ascending=False
        elif 'NOAA09' == self.name:
            ascending=True
        elif 'NOAA10' == self.name:
            ascending=False
        elif 'NOAA11' == self.name:
            ascending=True
        elif 'NOAA12' == self.name:
            ascending=False
        elif 'NOAA14' == self.name:
            ascending=True
        elif 'NOAA15' == self.name:
            ascending=False
        elif 'NOAA16' == self.name:
            ascending=True
        elif 'NOAA17' == self.name:
            ascending=False
        elif 'NOAA18' == self.name:
            ascending=True
        elif 'NOAA19' == self.name:
            ascending=True
        elif 'METOPA' == self.name:
            ascending=False
        elif 'METOPB' == self.name:
            ascending=False
        else:
            raise Exception,"Cannot recognise self.name (equ crossing type"

        return ascending

    # Get equator crossing 1 and 2 in each direction
    def get_nearest_time(self,line1,line2,year,month,day):
        
        # Get TLE record
        tle_rec = ephem.readtle(self.name, line1, line2)

        # Start with day boundary
        # First work out if ascending or not
        time = datetime.datetime(year,month,day)
        time_add_day = time + datetime.timedelta(seconds=86400)
        # Now get times of equator crossing - depends on satellite for daytime case (needed for solar contamination)
        ascending=self.get_ascending_descending_type()

        # Now Loop round a complete day to get times of equator crossing
        # including first one in next day
        ok = False
        i=0
        time_arr = []
        while True:
            intime = time + datetime.timedelta(seconds=i*60)
            if intime >= time_add_day:
                break
            tle_rec.compute(intime)
            if 0 == i:
                start_lat = tle_rec.sublat
            else:
                if ascending:
                    if tle_rec.sublat >= 0. and start_lat <= 0.:
                        ok = True
                        time_arr.append(intime)
                else:
                    if tle_rec.sublat <= 0. and start_lat >= 0.:
                        ok = True
                        time_arr.append(intime)
            start_lat = tle_rec.sublat
            i=i+1

        # Get first cross over day boundary
        ok = False
        i=0
        while True:
            intime = time_add_day + datetime.timedelta(seconds=i*60)
            tle_rec.compute(intime)
            if 0 == i:
                start_lat = tle_rec.sublat
            else:
                if not ascending:
                    if tle_rec.sublat >= 0. and start_lat <= 0.:
                        ok = True
                        time_arr.append(intime)
                        break
                else:
                    if tle_rec.sublat <= 0. and start_lat >= 0.:
                        ok = True
                        time_arr.append(intime)
                        break
            start_lat = tle_rec.sublat
            i=i+1

        if not ok:
            raise Exception,"ERROR: Cannot find equator from TLE"
        return time_arr

    def __init__(self,instr,year,month,day):

        time = datetime.datetime(year,month,day)
        # Read TLE file and find closest line pair
        line1,line2 = self.read_tle(instr,year,month,day)
        # Get nearest equator crosstime time to time value
        self.times = self.get_nearest_time(line1,line2,year,month,day)

# Get list of files that contain low/high time
def __get_filelist(low_time,high_time,avhrr_dir_name):

    # MT: 24-10-2017: Import local ECTs to implement additional logical processing
    ect1_time = low_time + datetime.timedelta(seconds=5400) 
    ect2_time = high_time - datetime.timedelta(seconds=5400)

    # Make directory name
    dir_name_list = '/group_workspaces/cems2/esacci_sst/input/avhrr/l1b/{0}/v1/{1:04d}/{2:02d}/{3:02d}/NSS.*'.\
        format(avhrr_dir_name,low_time.year,\
                   low_time.month,low_time.day)
    
    # Get list of all files
    filelist = glob.glob(dir_name_list)

    # If low_time.day == high_time.day then everything is in the same day
    # If not then need second day
    add_day=False
    if low_time.day != high_time.day:
        dir_name_list = '/group_workspaces/cems2/esacci_sst/input/avhrr/l1b/{0}/v1/{1:04d}/{2:02d}/{3:02d}/NSS.*'.\
            format(avhrr_dir_name,high_time.year,\
                       high_time.month,high_time.day)

        filelist2 = glob.glob(dir_name_list)
        filelist = filelist+filelist2
        add_day=True

    # Check list against low_time/high_time to select required files
    stored_file = []
    ok = False
    for i in range(len(filelist)):
        # Get file data
        infile = get_file_data(filelist[i])
        # Check if there is overlap of file with low/high time        
        if infile.start_time <= low_time and infile.end_time >= low_time:
# MT: 24-10-2017: add logic for case where end of orbit overlaps start of ECT window
            if infile.end_time >= ect1_time: 
                ok = True
                stored_file.append(filelist[i])
        elif infile.start_time <= high_time and infile.end_time >= high_time:
# MT: 24-10-2017: add logic for case where start of orbit overlaps end of ECT window
            if infile.start_time <= ect2_time: 
                ok = True
                stored_file.append(filelist[i])
        elif infile.start_time >= low_time and infile.end_time <= high_time:
# MT: 24-10-2017: add logic for case where end of orbit overlaps start of ECT window
            if infile.start_time <= ect1_time and infile.end_time <= ect2_time: 
                ok = True
                stored_file.append(filelist[i])
# MT: 24-10-2017: end of orbit overlaps start of ECT window 
# MT: 15-11-2017: convert if to elif to fix repetition of orbit files in filelist
            elif infile.start_time >= ect1_time and infile.end_time >= ect2_time: 
                ok = True
                stored_file.append(filelist[i]) 
# MT: 24-10-2017: start of orbit overlaps end of ECT window
# MT: 15-11-2017: convert if to elif to fix repetition of orbit files in filelist
            elif infile.start_time >= ect1_time and infile.start_time <= ect2_time: 
                ok = True
                stored_file.append(filelist[i])

    return ok,add_day,stored_file

def __find_avhrr_list_good(filelist,year,month,day,instr):

    filename_list=[]
    for i in range(len(filelist)):
        tuple_path = os.path.split(filelist[i])
        filename_list.append(tuple_path[1])

    new_filelist=[]
    new_accepted=[]
    ok=False
    black_list = blacklist(instr,year,month,day)
    for i in range(len(filename_list)):
        if not black_list.check_blacklist(filename_list[i]):
            ok=True
            new_filelist.append(filelist[i])
            new_accepted.append(filename_list[i])

    # Return status and list
    return ok,new_filelist,new_accepted

# Makes dir structure and shell script for submission
# This is for a single run
def make_shell_command(filelist,instr,avhrr_dir_name,year,month,day,i,\
                           equ_time1,equ_time2,split_single,spawn_job,\
                           test=False,gbcs_l1c_args='N',\
                           walton_only=False,keep_temp=False,\
                           write_fcdr=True,walton_cal=False):
    
    curr_dir = os.getcwd()
    # Check to see if we've already made this directory
    outdir = '{0}/{1}/{2:04d}/{3:02d}/{4:02d}'.\
        format(curr_dir,avhrr_dir_name,year,\
                   month,day)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    os.chdir(outdir)
    # Make link to make_fcdr.exe
    try:
        os.symlink('/group_workspaces/cems2/fiduceo/Users/jmittaz/FCDR/Mike/FCDR_AVHRR/make_fcdr.exe','make_fcdr.exe')
    except:
        pass
    try:
        os.symlink('/group_workspaces/cems2/fiduceo/Users/jmittaz/FCDR/Mike/FCDR_AVHRR/write_easy_fcdr_from_netcdf.py','write_easy_fcdr_from_netcdf.py')
    except:
        pass
    try:
        os.symlink('/group_workspaces/cems2/fiduceo/Users/jmittaz/FCDR/Mike/FCDR_AVHRR/write_l1c_data.py','write_l1c_data.py')
    except:
        pass
    try:
        os.symlink('/group_workspaces/cems2/fiduceo/Users/jmittaz/FCDR/Mike/FCDR_AVHRR/write_easy_fcdr.sh','write_easy_fcdr.sh')
    except:
        pass
    try:
        os.symlink('/home/users/jpdmittaz/Python/jpdm/lib/python2.7/site-packages/pygac/gac_run.py','gac_run.py')
    except:
        pass
    # Write script files
    outfile = 'run.{0:06d}.sh'.format(i)
    file_log = 'run.{0:06d}.log'.format(i)
    with open(outfile,'w') as fp:
        fp.write('. /home/users/jpdmittaz/Python/jpdm/bin/activate\n')
        fp.write('export PYGAC_CONFIG_FILE=/home/users/jpdmittaz/pygac.cfg\n')
        outfile_stem = []
        for j in range(len(filelist)):
            # Convert listed data via pygac
            base = os.path.basename(filelist[j])
            str_tuple = os.path.splitext(base)
            out_file_stem = str_tuple[0]+'.'+str(uuid.uuid4())
            # Deal with different compressions
            gzip_file = False
            bzip_file = False
            if str_tuple[1] == '.gz':
                newstr = 'cp -f {0} {1}.gz\n'.format(filelist[j],\
                                                         out_file_stem)
                fp.write(newstr)
                newstr = 'gunzip -f {0}.gz\n'.format(out_file_stem)
                fp.write(newstr)
                gzip_file = True
            elif str_tuple[1] == '.bz2':
                newstr = 'cp -f {0} {1}.bz2\n'.format(filelist[j],\
                                                          out_file_stem)
                fp.write(newstr)
                newstr = 'bunzip2 -f {0}.bz2\n'.format(out_file_stem)
                fp.write(newstr)
                bzip_file = True
            else:
                newstr = 'cp -f {0} {1}\n'.format(filelist[j],\
                                                      out_file_stem)
                fp.write(newstr)
            # Check to see if we have an archive header and remove if
            # necessary            
            if gzip_file:
                with gzip.open(filelist[j],'r') as fphead:
                    header = fphead.read(512)
            elif bzip_file:
                with bz2.BZ2File(filelist[j],'r') as fphead:
                    header = fphead.read(512)
            else:
                with open(filelist[j],'r') as fphead:
                    header = fphead.read(512)
            #
            # Look for ALL in Lat/Lon start/stop fields of archive
            # header
            #
            if "ALL" == header[75:78] and "ALL" == header[78:81] and \
                    "ALL" == header[81:84] and "ALL" == header[85:88]:
                newstr='echo "Removing NOAA CLASS Header"\n'
                fp.write(newstr)
                #
                # Remove leading 512 bytes from file
                #
                newstr='dd skip=1 ibs=512 if={0} of={0}.tmp\n'.\
                    format(out_file_stem)
                fp.write(newstr)
                newstr='mv -f {0}.tmp {0}\n'.\
                    format(out_file_stem)
                fp.write(newstr)
            # Make temporary directory to run pygac in
            # This is so we can find the output filename
            tempDir = uuid.uuid4()
            newstr = 'mkdir -p {0}\n'.format(tempDir)
            fp.write(newstr)
            newstr = 'cd {0}\n'.format(tempDir)
            fp.write(newstr)
            newstr = 'ln -s ../{0} .\n'.format(out_file_stem)
            fp.write(newstr)
            fp.write('ln -s ../gac_run.py\n')
            newstr = 'python2.7 gac_run.py {0} 0 0\n'.\
                format(out_file_stem)
            fp.write(newstr)
            # Find first pygac file name (_avhrr_ case)
            newstr='pygac{0:1d}=`ls ECC_GAC_avhrr*.h5`\n'.format(j+1)
            fp.write(newstr)
            # Find second pygac file name (_avhrr_ case)
            newstr='pygac{0:1d}_2=`ls ECC_GAC_qualflags*.h5`\n'.format(j+1)
            fp.write(newstr)
            # Find third pygac file name (_avhrr_ case)
            newstr='pygac{0:1d}_3=`ls ECC_GAC_sunsatangles*.h5`\n'.format(j+1)
            fp.write(newstr)
            # copy pygac files back up            
            fp.write('mv -f ECC_GAC_*.h5 ..\n')
            # cd up a level
            fp.write('cd ..\n')
            # remove temporary directory
            newstr='rm -rf {0}\n'.format(tempDir)
            fp.write(newstr)
            # Get output filename
            outfile_stem.append(out_file_stem)
        # Write merge command with all files                    
        newstr = './make_fcdr.exe '+str(uuid.uuid4())+' '+\
            instr+' '+gbcs_l1c_args+' '
        # Add equator crossing time estimates
        newstr = newstr + '{0:04d} {1:02d} {2:02d} {3:02d} {4:02d} '.\
            format(equ_time1.year,equ_time1.month,equ_time1.day,\
                       equ_time1.hour,equ_time1.minute)
        newstr = newstr + '{0:04d} {1:02d} {2:02d} {3:02d} {4:02d} '.\
            format(equ_time2.year,equ_time2.month,equ_time2.day,\
                       equ_time2.hour,equ_time2.minute)
        if split_single:
            newstr = newstr + ' Y'
        else:
            newstr = newstr + ' N'
        if walton_cal:
            if walton_only:
                newstr = newstr + ' Y'
            else:
                newstr = newstr + ' N'
        else:
            newstr = newstr + ' F'
        if keep_temp:
            newstr = newstr + ' Y'
        else:
            newstr = newstr + ' N'
        if write_fcdr:
            newstr = newstr + ' Y'
        else:
            newstr = newstr + ' N'
        newstr = newstr + ' {0:2d}'.format(len(filelist))
        for j in range(len(filelist)):
            newstr = newstr+' '+outfile_stem[j]
        if len(filelist) == 1:
            newstr = newstr+' ${pygac1}'
        elif len(filelist) == 2:
            newstr = newstr+' ${pygac1} ${pygac2}'
        elif len(filelist) == 3:
            newstr = newstr+' ${pygac1} ${pygac2} ${pygac3}'
        elif len(filelist) == 4:
            newstr = newstr+' ${pygac1} ${pygac2} ${pygac3} ${pygac4}'
        elif len(filelist) == 5:
            newstr = newstr+' ${pygac1} ${pygac2} ${pygac3} ${pygac4} ${pygac5}'
        newstr = newstr+'\n'
        fp.write(newstr)

        for j in range(len(filelist)):
# MT: 12-04-2018: keep temp files for uncertainty component analysis
            newstr = 'rm -f '+outfile_stem[j]+'*\n'
            fp.write(newstr)
#        fp.write('rm -f *.h5')
#        fp.write('rm -f make_fcdr.exe')
#        fp.write('rm -f gac_run.py')
#        fp.write('rm -f write_easy_fcdr_from_netcdf.py')
            if 0 == j:
                fp.write('rm -f ${pygac1}\n')
                fp.write('rm -f ${pygac1_2}\n')
                fp.write('rm -f ${pygac1_3}\n')
            elif 1 == j:
                fp.write('rm -f ${pygac2}\n')
                fp.write('rm -f ${pygac2_2}\n')
                fp.write('rm -f ${pygac2_3}\n')
            elif 2 == j:
                fp.write('rm -f ${pygac3}\n')
                fp.write('rm -f ${pygac3_2}\n')
                fp.write('rm -f ${pygac3_3}\n')
            elif 3 == j:
                fp.write('rm -f ${pygac4}\n')
                fp.write('rm -f ${pygac4_2}\n')
                fp.write('rm -f ${pygac4_3}\n')
            elif 4 == j:
                fp.write('rm -f ${pygac5}\n')
                fp.write('rm -f ${pygac5_2}\n')
                fp.write('rm -f ${pygac5_3}\n')

    # submit jobs
    os.chmod(outfile,stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)
    job_name='./'+outfile
    job = ['bsub','-q', 'short-serial','-W', '01:00','-oo', file_log, job_name]
    # Actually submit jobs
    if spawn_job:
        subprocess.call(job)
    os.chdir(curr_dir)

# Write all shell command scripts for complete day
def write_commands(instr,year,month,day,split_single,spawn_job,\
                       timestep=60,test=False,\
                       gbcs_l1c_args='N',\
                       walton_only=False,\
                       keep_temp=False,\
                       write_fcdr=True.\
                       walton_cal=False):
    
    # Get equator crossing times in the day
    t = tle_data(instr,year,month,day)
    if not t.times:
        raise Exception,"No TLE equator crossing times found"

    avhrr_dir_name = get_avhrr_dir_name(instr)

    # Loop round crossing times
    nwrites=0
    for eqtr in range(len(t.times)-1):
        # start and end time of equator crossing and 1/2 hour timewindow
        stime = t.times[eqtr] - datetime.timedelta(seconds=5400)
        etime = t.times[eqtr+1] + datetime.timedelta(seconds=5400)
        # Get files that contain stime/etime
        ok,add_day,list_of_files = \
            __get_filelist(stime,etime,avhrr_dir_name)
        
        if ok:
            # Apply blacklist
            ok2,filelist,accepted = __find_avhrr_list_good(list_of_files,\
                                                               year,month,day,\
                                                               instr)
            make_shell_command(filelist,instr,avhrr_dir_name,year,month,day,\
                                   eqtr,t.times[eqtr],t.times[eqtr+1],\
                                   split_single,spawn_job,test=test,\
                                   gbcs_l1c_args=gbcs_l1c_args,\
                                   walton_only=walton_only,keep_temp=keep_temp,\
                                   write_fcdr=write_fcdr,walton_cal=walton_cal)
            nwrites=nwrites+1

    print 'Number of command files : ',nwrites

if __name__ == "__main__":

    parser = OptionParser("usage: %prog instr year month day split_single_file GBCS_L1C(Y/N/C) Spawn_Jobs(Y/N) walton_only(Y/N/F) keep_temp(Y/N) write_fcdr{Y/N) (test=Y/N)")
    (options, args) = parser.parse_args()
    if len(args) != 10 and len(args) != 11:
        parser.error("incorrect number of arguments")
    year = int(args[1])
    month = int(args[2])
    day = int(args[3])
    split_single_file = args[4]
    gbcs_l1c_args = args[5]
    if args[6] == 'Y':
        spawn_job=True
    else:
        spawn_job=False
    walton_cal=True
    if args[7] == 'Y':
        walton_only=True
    else:
        walton_only=False
        if args[7] == 'F':
            walcon_cal=False
    if args[8] == 'Y':
        keep_temp=True
    else:
        keep_temp=False
    if args[9] == 'Y':
        write_fcdr=True
    else:
        write_fcdr=False
    if len(args) == 11:
        if args[10] == 'Y':
            test=True
        else:
            test=False
    else:
        test=False

    if 'Y' == split_single_file or 'y' == split_single_file:
        split_single=True
    else:
        split_single=False

    write_commands(args[0],year,month,day,split_single,spawn_job,test=test,\
                       gbcs_l1c_args=gbcs_l1c_args,walton_only=walton_only,\
                       keep_temp=keep_temp,write_fcdr=write_fcdr,\
                       walton_cal=walton_cal)


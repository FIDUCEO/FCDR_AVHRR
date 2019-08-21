from __future__ import print_function
import equator_to_equator_ensemble as eq
import os
import glob
import subprocess
import datetime

def find_harmonisation(instr,year,month,day,harmonisation_dir):
    '''Find relevant harmonisation files for given date/sensor'''

    if 'TIROSN' != instr and \
            'NOAA06' != instr and \
            'NOAA07' != instr and \
            'NOAA08' != instr and \
            'NOAA09' != instr and \
            'NOAA10' != instr and \
            'NOAA11' != instr and \
            'NOAA12' != instr and \
            'NOAA14' != instr and \
            'NOAA15' != instr and \
            'NOAA16' != instr and \
            'NOAA17' != instr and \
            'NOAA18' != instr and \
            'NOAA19' != instr and \
            'METOPA' != instr and \
            'METOPB' != instr:
        print(instr)
        raise Exception("Incorrect avhrr_name")

    #
    # Find Harmonisation files
    # Based on sensor name
    #
    harmonisation_files = \
        glob.glob('{0}/FIDUCEO_Harmonisation_Data_{1}_*.nc'.\
                      format(harmonisation_dir,instr))

    if None == harmonisation_files:
        raise Exception("No matching harmonisation files found (instrument name)")
    
    #
    # Get times/channel from filenames
    #
    start = []
    end = []
    channel = []
    for filename_with_dir in harmonisation_files:
        filename = os.path.basename(filename_with_dir)
        syear = int(filename[37:41])
        sday = int(filename[41:44])
        start_string = '{0:04d}-{1:03d}'.format(syear,sday)
        eyear = int(filename[45:49])
        eday = int(filename[49:52])
        end_string = '{0:04d}-{1:03d}'.format(eyear,eday)
        start.append(datetime.datetime.strptime(start_string,'%Y-%j'))
        end.append(datetime.datetime.strptime(end_string,'%Y-%j'))
        channel.append(int(filename[34:36]))

    #
    # For each channel find associated file in time range
    #
    time_string = '{0:04d}-{1:02d}-{2:02d}'.format(year,month,day)
    intime = datetime.datetime.strptime(time_string,'%Y-%m-%H')
    ch3b_there = False
    ch4_there = False
    ch5_there = False
    for i in range(len(channel)):
        if start[i] <= intime and end[i] >= intime:
            if channel[i] == 37:
                ch3b_there = True
                ch3b_name = harmonisation_files[i]
            elif channel[i] == 11:
                ch4_there = True
                ch4_name = harmonisation_files[i]
            elif channel[i] == 12:
                ch5_there = True
                ch5_name = harmonisation_files[i]

    if not ch3b_there:
        raise Exception('CH3B Harmonisation file not present for this date')
    if not ch4_there:
        raise Exception('CH4 Harmonisation file not present for this date')
    if not ch5_there:
        raise Exception('CH5 Harmonisation file not present for this date')

    #
    # Find MC delta parameter file(s)
    #
    MC_files = \
        glob.glob('{0}/FIDUCEO_Harmonisation_Data_MC_{1}*.nc'.\
                      format(harmonisation_dir,instr))

    if None == MC_files:
        raise Exception("No matching harmonisation files found (instrument name)")
    #
    # Find those MC files in time zone
    #
    start = []
    end = []
    for filename_with_dir in MC_files:
        filename = os.path.basename(filename_with_dir)
        syear = int(filename[37:41])
        sday = int(filename[41:44])
        start_string = '{0:04d}-{1:03d}'.format(syear,sday)
        eyear = int(filename[45:49])
        eday = int(filename[49:52])
        end_string = '{0:04d}-{1:03d}'.format(eyear,eday)
        start.append(datetime.datetime.strptime(start_string,'%Y-%j'))
        end.append(datetime.datetime.strptime(end_string,'%Y-%j'))

    #
    # For each channel find associated file in time range
    #
    MC_there = False
    for i in range(len(start)):
        if start[i] <= intime and end[i] >= intime:
            MC_there = True
            MC_name = MC_files[i]

    if not MC_there:
        raise Exception('MC Harmonisation file not present for this date')
    #
    # Return files
    #
    return ch3b_name,ch4_name,ch5_name,MC_name

def run(instr,year,month,day,in_harmonisation_dir,oldHarm=False):
    '''Run an individual sensor/year/month/day to generate baseline and 
    Ensemble FCDR files'''

    #
    # Expand harmonisation_dir to full path
    #
    harmonisation_dir = os.path.abspath(in_harmonisation_dir)

    #
    # Get appropriate files for Harmonisation
    #
    harm_37,harm_11,harm_12,fiduceo_mc_harm = \
        find_harmonisation(instr,year,month,day,harmonisation_dir)

    #
    # Write run.*.sh scripts
    #
    eq.main_python(instr,year,month,day,"Y N N M N Y N",fiduceo_mc_harm)

    dirname = '{0}/{1:04d}/{2:02d}/{3:02d}'.\
        format(eq.get_avhrr_dir_name(instr),year,month,day)
    dirname = os.path.abspath(dirname)

    #
    # Process all run.*.sh files in one go as part of this process
    #
    os.chdir(dirname)
    filelist = glob.glob('run.*.sh')
    #
    # Write master scripts including harmonisation files
    #
    if filelist != None:
        with open('master_scripts.sh','w') as fp:
            #
            # If new calibration set evironment parameters
            #
            if not oldHarm:
                fp.write('export FIDUCEO_Harmonisation_Data_37={0}\n'.\
                             format(harm_37))
                fp.write('export FIDUCEO_Harmonisation_Data_11={0}\n'.\
                             format(harm_11))
                fp.write('export FIDUCEO_Harmonisation_Data_12={0}\n'.\
                             format(harm_12))
            fp.write('cd {0}\n'.format(dirname))
            for filename in filelist:
                fp.write('. ./{0}\n'.format(filename))

            if not oldHarm:
                fp.write('unset FIDUCEO_Harmonisation_Data_37\n'.\
                             format(harm_37))
                fp.write('unset FIDUCEO_Harmonisation_Data_11\n'.\
                             format(harm_11))
                fp.write('unset FIDUCEO_Harmonisation_Data_12\n'.\
                             format(harm_12))
        job = ['.', 'master_scripts.sh', '>', 'master_scripts.log', '2>&1']
#        subprocess.call(job)

    os.chdir('../../../..')
    return dirname

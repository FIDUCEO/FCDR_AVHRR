
#!/usr/bin/env python

#La difference par rapport a l'autre version est le repertoire de sortie, ici group workspace fiduceo

import os
import os.path
import optparse
import sys
import subprocess
import cemsjob
import numpy as np
import func
from func import *

__version__="1.0"
instrument=sys.argv[1]
if sys.argv[2]=="all":
	year_start, year_end=year(instrument)
	m_beg=1
	m_end=12
else:
	year_start=int(sys.argv[2])
	year_end=int(sys.argv[3])
	m_beg=int(sys.argv[4])
	m_end=int(sys.argv[5])
path_github="/home/users/mdesmons/gbcs/src/AVHRR/"
path="/group_workspaces/cems2/fiduceo/Users/mdesmons/avhrr_l1b/"
path_scratch="/work/scratch/mdesmons/"
path_out=path_scratch+"temp/test_"+instrument+"/"
#if not os.path.isdir(dir_out):
#        os.mkdir (dir_out)

for i_year in range(year_start,year_end+1):
       year= str("%04d" % i_year)
       dir_out=path_out+year+"/"
       if not os.path.isdir(dir_out):
               #print "dir_out doesn exist"
               os.mkdir(dir_out)
       for i_month in range(m_beg, m_end+1):
              month =  str("%02d" % i_month)
              print "year, month", i_year, i_month 
              dir_in=path+"listes_orbites/liste_orbites_day/"+instrument+"/"+year+"/"+month+"/"
              dir_out=path_out+year+"/"+month+"/"
              #print "dir_out", dir_out
              if not os.path.isdir(dir_out):
                      #print "dir_out doesn exist"
                      os.mkdir(dir_out)
              for i_day in range(1, 32):
                     day =  str("%02d" % i_day)
                     dir_out=path_out+year+"/"+month+"/"+day+"/"
                     #print "dir_out", dir_out
                     if not os.path.isdir(dir_out):
                             #print "dir_out doesn exist"
                             os.mkdir(dir_out)
                     file_in=dir_in+"liste_orbites_"+instrument+"_"+year+"_"+month+"_"+day+".data"
                     file_in_2=dir_in+"blackliste_"+instrument+"_"+year+"_"+month+"_"+day+".data"
                     print file_in_2
                     if os.path.isfile(file_in)==True: 
                            f = open(file_in, 'r')
                            n = 0
                            for line in f:
                                   n += 1
                            print "nb lignes dans file_in",n     
                            if n>1:    
                                    if os.path.isfile(file_in_2)==True: 
                                            g = open(file_in_2, 'r')
                                            ng = 0
                                            for line in g:
                                                    ng += 1
                                            print "nb lignes dans file_in_2",ng     
#----------------------------------On lit la liste des orbites  
                                    liste_orbites= np.loadtxt(file_in,unpack=True, dtype=np.str)
                                    for i in range(len(liste_orbites)):
                                            print "i", i, "/", n
					    file_out=dir_out+"test_"+liste_orbites[i][84:99]+".out"
                                            file_err=dir_out+"test_"+liste_orbites[i][84:99]+".err"
                                            if os.path.isfile(file_out):
                                                os.remove (file_out)
                                            if os.path.isfile(file_err):
                                                os.remove (file_err)
                                            job = ['bsub',
                                            '-q', 'short-serial',
                                            "-W", "00:40",
                                            '-o', file_out,
                                            '-e', file_err,
                                            path_github+"check_37_chan_sortie_par_orbite_with_channel_3_4_5.FCDR.essai_other_output.exe",
                                            instrument, str(i_year), str(i_month), str(i_day), "easy_fcdr", str(liste_orbites[i])]
#-------------------------------------------Si la blackliste n'est pas vide, on la lit
                                            if ng>2:        
                                            	blackliste,redundant,reason= np.loadtxt(file_in_2,unpack=True,dtype=np.str,usecols=(0,2,3))
                                                a=liste_orbites[i][75:120]
                                                if any(x==a for x in blackliste): 
							b=np.where(blackliste==a)
                                                        print blackliste[b],reason[b]
                                                        if reason[b]=="too_small" \
							or reason[b]=="too_long" \
							or reason[b]=='bad_l1c_quality'\
							or reason[b]=="ground_station_duplicate" \
							or reason[b]=='along_track_too_long'\
							or redundant[b]=='1':
                                                        	print "black", blackliste[b], reason[b] 
                                                        else:
								print job
								subprocess.call(job)	
	                                        else:
							print job
							subprocess.call(job)
                                            else:
						print job
                                                subprocess.call(job)
# 2 usages posssibles:
# all : si on veut faire tourner le script sur toute la duree de l'instrument
# year_start, year_end, m_start, m_end: si on veut faire tourner le script sur une duree plus specifique

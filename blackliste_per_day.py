import os
import sys
import sqlite3
import func
import re
import numpy as np

instrument=sys.argv[1]

year_start,year_end=func.year(instrument)
print year_start,year_end
ins,ins2=func.instrument_to_ins(instrument)

dir_instrument="/group_workspaces/cems2/fiduceo/Users/mdesmons/avhrr_l1b/listes_orbites/liste_orbites_day/"+instrument+"/"
if not os.path.isdir(dir_instrument):
    os.mkdir(dir_instrument) 

######################################
sqlite_file = '/home/users/mdesmons/gbcs/src/AVHRR/AVHRR_GAC_archive_v2_201603_post_overlap.sqlite3'    # name of the sqlite database file
table_name = 'vw_std'   # name of the table to be queried
column_1 = 'filename'
column_2 = 'blacklist'
column_3 = 'satellite_name'
column_4 = 'start_time_l1b'
column_5 ='end_time_l1b'
column_6="missing_scanlines"
column_7="redundant"
column_8="blacklist_reason"

conn=sqlite3.connect(sqlite_file)
c = conn.cursor()            

query = """
     SELECT
              {fl},{cn},{dn},{ms}
            FROM
              {tn}
            WHERE
              ({cn}=1 OR
              {dn}=1) AND
              {sn}=? AND
              {st} >? AND
              {et} <?
""".format(fl=column_1,ms=column_8,tn=table_name, cn=column_2,dn=column_7,sn=column_3, st=column_4, et=column_5)
#####################################

for y in range(year_start,year_end+1):
    print y
    dir_year=dir_instrument+str(y)+"/"
    if not os.path.isdir(dir_year):
        os.mkdir(dir_instrument)
    for m in range(1,13):
        output_dir=dir_year+"%2.2d"%m+"/"
        if not os.path.isdir(output_dir):
            os.mkdir(dir_month)
        for d in range(1,32):
            output_file_1=output_dir+"blackliste_"+instrument+"_"+str(y)+"_"+str("%2.2d"%m)+"_"+str("%2.2d"%d)+".data"
            if os.path.isfile(output_file_1):
                os.remove(output_file_1)
                print "rm"

            file_to_write=open(output_file_1,"w")
            start=str(y)+"-"+str("%2.2d"%m)+"-"+str("%2.2d"%d)+" 00:00"
            end=str(y)+"-"+str("%2.2d"%m)+"-"+str("%2.2d"%d)+" 23:59"
            print ins2,start,end
           
            c.execute(query,(ins2,start,end,))
            all_rows = c.fetchall()
            #print all_rows
            for i in range(len(all_rows)):
                a=all_rows[i]
                b_1 = str(all_rows[i]).replace("(u'", ' ')
                b_2 = b_1.replace("')", ' ')
                b_3 = b_2.replace("'", ' ')
                b_4 = b_3.replace("u", ' ')
                b_5 = b_4.replace(" ", '')
                b_6 = b_5.replace(",", ' ')
                file_to_write.write("\n"+b_6)
                #print b_6
               
          

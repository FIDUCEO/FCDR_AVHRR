find AVHRR06_G -name '*.log' -print -exec grep CRITICAL {} \; > /tmp/logs_06.log
find AVHRR07_G -name '*.log' -print -exec grep CRITICAL {} \; > /tmp/logs_07.log
find AVHRR08_G -name '*.log' -print -exec grep CRITICAL {} \; > /tmp/logs_08.log
find AVHRR09_G -name '*.log' -print -exec grep CRITICAL {} \; > /tmp/logs_09.log
find AVHRR10_G -name '*.log' -print -exec grep CRITICAL {} \; > /tmp/logs_10.log
find AVHRR11_G -name '*.log' -print -exec grep CRITICAL {} \; > /tmp/logs_11.log
find AVHRR12_G -name '*.log' -print -exec grep CRITICAL {} \; > /tmp/logs_12.log
find AVHRR14_G -name '*.log' -print -exec grep CRITICAL {} \; > /tmp/logs_14.log
find AVHRR15_G -name '*.log' -print -exec grep CRITICAL {} \; > /tmp/logs_15.log
find AVHRR16_G -name '*.log' -print -exec grep CRITICAL {} \; > /tmp/logs_16.log
find AVHRR17_G -name '*.log' -print -exec grep CRITICAL {} \; > /tmp/logs_17.log
find AVHRR18_G -name '*.log' -print -exec grep CRITICAL {} \; > /tmp/logs_18.log
find AVHRR19_G -name '*.log' -print -exec grep CRITICAL {} \; > /tmp/logs_19.log
find AVHRRMTA_G -name '*.log' -print -exec grep CRITICAL {} \; > /tmp/logs_MTA.log

#grep -o 'CRITICAL : Cannot find equator where it should be (first)' logs_06.log | wc -l
#grep -o 'CRITICAL : No overlap between file1 and file2' logs_06.log | wc -l
#grep -o 'CRITICAL : No overlap between file2 and file3' logs_06.log | wc -l
#grep -o 'CRITICAL : No overlap between file3 and file4' logs_06.log | wc -l
#grep -o 'CRITICAL : No overlap between file4 and file5' logs_06.log | wc -l
#grep -o 'CRITICAL : Cannot find any good Tinstr data' logs_06.log | wc -l
#grep -o 'CRITICAL : USAGE: ./extract_l1b_data.exe uuid outfile eq_year1 eq_month1 eq_day1 eq_hour1 eq_min1 eq_year2 eq_month2 eq_day2 eq_hour2 eq_min2 split_single file1 (file2) (file3) (file4) (file5)' logs_06.log | wc –l
#grep -o 'CRITICAL : No data available' logs_06.log | wc -l
#grep -o 'CRITICAL : Cannot find time index' logs_06.log | wc –l



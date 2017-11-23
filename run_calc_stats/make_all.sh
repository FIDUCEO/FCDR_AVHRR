rm -f AVHRR*_G_*.dat
echo ./make_stats.sh AVHRR06_G > make1.sh
bsub -q short-serial -W24:00 -oo make1.log < make1.sh
echo ./make_stats.sh AVHRR07_G > make2.sh
bsub -q short-serial -W24:00 -oo make2.log < make2.sh
echo ./make_stats.sh AVHRR08_G > make3.sh
bsub -q short-serial -W24:00 -oo make3.log < make3.sh
echo ./make_stats.sh AVHRR09_G > make4.sh
bsub -q short-serial -W24:00 -oo make4.log < make4.sh
echo ./make_stats.sh AVHRR10_G > make5.sh
bsub -q short-serial -W24:00 -oo make5.log < make5.sh
echo ./make_stats.sh AVHRR11_G > make6.sh
bsub -q short-serial -W24:00 -oo make6.log < make6.sh
echo ./make_stats.sh AVHRR12_G > make7.sh
bsub -q short-serial -W24:00 -oo make7.log < make7.sh
echo ./make_stats.sh AVHRR14_G > make8.sh
bsub -q short-serial -W24:00 -oo make8.log < make8.sh
echo ./make_stats.sh AVHRR15_G > make9.sh
bsub -q short-serial -W24:00 -oo make9.log < make9.sh
echo ./make_stats.sh AVHRR16_G > make10.sh
bsub -q short-serial -W24:00 -oo make10.log < make10.sh
echo ./make_stats.sh AVHRR17_G > make11.sh
bsub -q short-serial -W24:00 -oo make11.log < make11.sh
echo ./make_stats.sh AVHRR18_G > make12.sh
bsub -q short-serial -W24:00 -oo make12.log < make12.sh
echo ./make_stats.sh AVHRR19_G > make13.sh
bsub -q short-serial -W24:00 -oo make13.log < make13.sh
echo ./make_stats.sh AVHRRMTA_G > make14.sh
bsub -q short-serial -W24:00 -oo make14.log < make14.sh

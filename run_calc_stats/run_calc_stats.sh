#!/bin/sh

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

echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.1.sh
echo python2.7 run_calc_stats.py NOAA06 >> run.1.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.2.sh
echo python2.7 run_calc_stats.py NOAA07 >> run.2.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.3.sh
echo python2.7 run_calc_stats.py NOAA08 >> run.3.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.4.sh
echo python2.7 run_calc_stats.py NOAA09 >> run.4.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.5.sh
echo python2.7 run_calc_stats.py NOAA10 >> run.5.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.6.sh
echo python2.7 run_calc_stats.py NOAA11 >> run.6.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.7.sh
echo python2.7 run_calc_stats.py NOAA12 >> run.7.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.8.sh
echo python2.7 run_calc_stats.py NOAA14 >> run.8.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.9.sh
echo python2.7 run_calc_stats.py NOAA15 >> run.9.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.10.sh
echo python2.7 run_calc_stats.py NOAA16 >> run.10.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.11.sh
echo python2.7 run_calc_stats.py NOAA17 >> run.11.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.12.sh
echo python2.7 run_calc_stats.py NOAA18 >> run.12.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.13.sh
echo python2.7 run_calc_stats.py NOAA19 >> run.13.sh
echo '. /home/users/mtaylor/pyenv/pygac/bin/activate' > run.14.sh
echo python2.7 run_calc_stats.py METOPA >> run.14.sh
bsub -q short-serial -W24:00 -oo run.1.log < run.1.sh
bsub -q short-serial -W24:00 -oo run.2.log < run.2.sh
bsub -q short-serial -W24:00 -oo run.3.log < run.3.sh
bsub -q short-serial -W24:00 -oo run.4.log < run.4.sh
bsub -q short-serial -W24:00 -oo run.5.log < run.5.sh
bsub -q short-serial -W24:00 -oo run.6.log < run.6.sh
bsub -q short-serial -W24:00 -oo run.7.log < run.7.sh
bsub -q short-serial -W24:00 -oo run.8.log < run.8.sh
bsub -q short-serial -W24:00 -oo run.9.log < run.9.sh
bsub -q short-serial -W24:00 -oo run.10.log < run.10.sh
bsub -q short-serial -W24:00 -oo run.11.log < run.11.sh
bsub -q short-serial -W24:00 -oo run.12.log < run.12.sh
bsub -q short-serial -W24:00 -oo run.13.log < run.13.sh
bsub -q short-serial -W24:00 -oo run.14.log < run.14.sh


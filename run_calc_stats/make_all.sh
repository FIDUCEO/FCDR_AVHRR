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

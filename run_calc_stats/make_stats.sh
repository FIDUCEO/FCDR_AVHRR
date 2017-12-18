#!/bin/bash

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

function get_channel {

    temp_name=`uuid`
    search_name=${1}'_*_'${2}'_'${3}'.data'
    find ${1} -name "${search_name}" > ${temp_name}.list
    rm -f ${1}_${2}_${3}.dat
    while read line
    do
	cat ${line} >> ${1}_${2}_${3}.dat
    done < ${temp_name}.list
    rm -f ${temp_name}.list

}

if [ $# -ne 1 ]
then
    echo "USAGE: AVHRR"
    exit -1
fi

get_channel ${1} independent Ch1
get_channel ${1} independent Ch2
get_channel ${1} independent Ch3a
get_channel ${1} independent Ch3b
get_channel ${1} independent Ch4
get_channel ${1} independent Ch5

get_channel ${1} structured Ch1
get_channel ${1} structured Ch2
get_channel ${1} structured Ch3a
get_channel ${1} structured Ch3b
get_channel ${1} structured Ch4
get_channel ${1} structured Ch5

get_channel ${1} measurement Ch1
get_channel ${1} measurement Ch2
get_channel ${1} measurement Ch3a
get_channel ${1} measurement Ch3b
get_channel ${1} measurement Ch4
get_channel ${1} measurement Ch5



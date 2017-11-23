#!/bin/bash

function get_channel {

    temp_name=`uuid`
    search_name='stats_*_'${2}'_'${3}'.data'
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



#!/bin/bash

if [ $# -ne 1 ] && [ $# -ne 2 ]
then
    echo "USAGE: ./write_easy_fcdr.sh temp_file (optional: output_filename)"
    exit -1
fi

# Define where Anaconda environment is with Gerrits CURUC code
#anaconda_env='/group_workspaces/cems2/fiduceo/Users/jmittaz/Anaconda/bin/'
anaconda_env='/gws/nopw/j04/fiduceo/Users/jmittaz/Anaconda/etc/profile.d/'

# Find if running jpdm virtual environment
python_there=0
anaconda_there=0
if [ -z $VIRTUAL_ENV ] 
then
    # Check to see if jpdm version of anaconda is running
    if conda info > /dev/null 
    then
	python_loc=`which python`
	if echo $python_loc | grep Users/jmittaz/Anaconda
	then
	    # jpdm version is running
	    anaconda_there=1
	else
	    # other conda - deactivate
	    conda deactivate
	fi
    fi
else
    python_loc=$VIRTUAL_ENV
    if echo $VIRTUAL_ENV | grep jpdm 
    then
	python_there=1	
    fi
fi

# If there then deactivate to start Gerrits environment
if [ 1 -eq ${python_there} ]
then
    deactivate
fi

# Start anaconda/Gerrits environment for CURUC
#export PATH=/group_workspaces/cems2/fiduceo/Users/jmittaz/Anaconda/bin:$PATH
#source /group_workspaces/cems2/fiduceo/Users/jmittaz/Anaconda/bin/activate py36_gh
# start setup of anaconda
if [ 0 -eq ${anaconda_there} ]
then
    echo '...Setting up Gerrits python3 enviroment for CURUC...'
    . /gws/nopw/j04/fiduceo/Users/jmittaz/Anaconda/etc/profile.d/conda.sh
    conda activate py36_gh
fi

# Run python writer
if [ $# -eq 1 ]
then
    python3 write_easy_fcdr_from_netcdf.py ${1}
else
    python3 write_easy_fcdr_from_netcdf.py ${1} ${2}
fi

# Exit anaconda/Gerrits environment for CURUC
#source /group_workspaces/cems2/fiduceo/Users/jmittaz/Anaconda/bin/deactivate
conda deactivate

# Turn back old environment (python 2)
if [ 1 -eq ${python_there} ]
then
    . ${python_loc}/bin/activate
fi

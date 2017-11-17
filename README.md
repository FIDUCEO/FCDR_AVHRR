
 * Copyright (C) 2017 J.Mittaz University of Reading
 * This code was developed for the EC project: Fidelity and Uncertainty in   
 * Climate Data Records from Earth Observations (FIDUCEO).
 * Grant Agreement: 638822
 * <Version> Reviewed and approved by <name, instituton>, <date>
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 * 
 * A copy of the GNU General Public License should have been supplied along
 * with this program; if not, see http://www.gnu.org/licenses/

# FCDR_AVHRR
Development code for AVHRR FCDR Uncertainties

Written by J.Mittaz, University of Reading 31/08/2017
Revised by M.Taylor, University of Reading 17/11/2017

Code here will create the pre-beta FIDUCEO FCDR. The data will be 
equator-to-equator data with all fields defined as in the specification. Note 
that the Fortran code is dependent on having the UoR GBCS (Generalised Bayesian
 Cloud Screening) installed from the ESA CCI SST project. 

The files in the directory are

equator-to-equator.py : top level script to run the process. Accepts 
instrument name and year month day and will create all FCDR files for that day.

write_easy_fcdr_from_netcdf.py : python code to convert temporary file created
by write_fcdf.exe into FIDUCEO netcdf format using Tom Block's writer available
from FIDUCEO/FCDRTools/.

Makefile: Makefile set up to make .exe file on CEMS. The Makefile assumes that 
the GBCS is installed as a GBCS directory in the source directory.

write_fcdr.f90 : Top level program to make FCDR data from an input list of 
AVHRR Level 1B files.

combine_orbits.f90 : Code to merge AVHRR L1B files rogether into a continuous
stream replacing bad scan lines from one with good scan lines from another file.
Then cuts out the equator-to-equator data and recalibrates it.

fiduceo_uncertainties.f90 : Code to create uncertainties for easy FCDR. 
Also specifies sensor specific channel correlation matrices and spatial
correlation scale. Writes temporary output netcdf file which is then converted using
write_easy_fcdr_from_netcdf.py within the code.

quickplot.py : Code to automate quick time series plot of reflectance and IR channel data (chx.png) as well as the independent uncertainties (chx_independent.png) and structured uncertainties (chx_structured.png).

run_instrument.py : script to run equator_to_equator.py over a sensor series

run.sh : shell script to spawn run_instrument.py from within a python environment over the archive

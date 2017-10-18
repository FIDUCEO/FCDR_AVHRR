# * Copyright (C) 2017 J.Mittaz University of Reading
# * This code was developed for the EC project “Fidelity and Uncertainty in   
# * Climate Data Records from Earth Observations (FIDUCEO)”. 
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

FC = gfortran
FL = gfortran
#FC = ifort
#FL = ifort
NCDFPATH = /group_workspaces/cems2/esacci_sst/software/common.gfortran
#NCDFPATH = /group_workspaces/cems2/esacci_sst/software/common
#FC_FLAGS = -I/group_workspaces/cems2/esacci_sst/software/common/include -I. -fpp -free -fPIC -assume byterecl -fpconstant -std03 -check all -debug all -traceback -g -O2
#F_FLAGS = -fpp -fPIC -assume byterecl -fpconstant -check all -debug all -traceback -g -O2
#FC_FLAGS = -Wall -g -fbounds-check -I. -I$(NCDFPATH)/include -cpp 
#FC_FLAGS = -g -fbounds-check -I. -I$(NCDFPATH)/include -cpp 
FC_FLAGS = -O -I. -I$(NCDFPATH)/include -cpp 
C_FLAGS = -I GBCS/thirdparty/epr_api/src -IGBCS/src/EPR_API
LIBRARIES = $(NCDFPATH)/lib/libnetcdff.a $(NCDFPATH)/lib/libnetcdf.a $(NCDFPATH)/lib/libhdf5hl_fortran.a $(NCDFPATH)/lib/libhdf5_hl.a $(NCDFPATH)/lib/libhdf5_fortran.a $(NCDFPATH)/lib/libhdf5.a -ldl -lz -LGBCS/thirdparty/epr_api/obj/Linux_x86_64 -lepr
#LIBRARIES += -L/usr/lib/gcc/x86_64-redhat-linux/4.1.1 -lgfortran -lpng 
LIBRARIES += -L/usr/lib/gcc/x86_64-redhat-linux/4.1.1 -lpng 

GBCSOBJECTS = GbcsKinds.o GbcsBaseTypes.o GbcsConstants.o GbcsStringUtil.o \
	GbcsDateTime.o GbcsPath.o GbcsPixel.o GbcsErrorHandler.o GbcsPDFLookup.o \
	GbcsConfFile.o GbcsLUTShift.o GbcsTypes.o GbcsAtmPhysics.o GbcsGeopotential.o GbcsProfileGenerator.o \
	GbcsTimeSpace.o GbcsSystemTools.o GbcsInterpolators.o GbcsForecastModel.o \
	GbcsMatrixOps.o GbcsPDFLoaders.o \
	GbcsScattPhysics.o GbcsCleanUp.o \
	GbcsPixelLoaders.o AVHRR_Filter_Data.o NOAA_LoadAVHRRLevel1B.o GbcsImageUtil.o \
	epr_wrapper_c.o epr_wrapper.o ARC_ATSRVarNEdT.o ARC_ATSR1.o ARC_L1bCorrection.o ARC_LoadImagery.o

OBJECTS = $(GBCSOBJECTS) fiduceo_uncertainties.o combine_orbits.o
OBJECTS_ALL1 =  $(OBJECTS) write_fcdr.o
OBJECTS_ALL2 =  $(OBJECTS) write_fcdr_single.o

vpath %.f90 GBCS/src/GbcsMod_Globals GBCS/src/GbcsMod_SystemTools \
	GBCS/src/GbcsMod_ErrorHandler GBCS/src/GbcsMod_CleanUp \
	GBCS/src/GbcsMod_Misc GBCS/src/GbcsMod_Bayesian \
	GENERIC IASI_L1 GBCS/src/GbcsMod_Profiles \
	GBCS/src/GbcsMod_ProcessImagery GBCS/src/GbcsMod_Physics \
	GBCS/src/GbcsMod_Maths GBCS/src/GbcsMod_Initialize \
	GBCS/src/GbcsMod_DataReaders GBCS/src/GbcsMod_RTM \
	GBCS/src/AVHRR GBCS/src/ARC GBCS/thirdparty/epr_api/src GBCS/src/EPR_API
vpath %.c GBCS/src/EPR_API
vpath %.f .

all: $(OBJECTS_ALL1)
	$(FC) -o write_fcdr.exe write_fcdr.o $(OBJECTS) $(LIBRARIES)

.PHONY: clean
clean:
	rm -rf *.o *.mod *.exe

include GBCS/src/make.rules

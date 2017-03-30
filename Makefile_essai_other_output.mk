FL = gfortran
FC = gfortran
NCDFPATH=/glusterfs/surft/users/pf906247/mySoft

FC_FLAGS = -g -fbounds-check -I /usr/include -I. -I$(NCDFPATH)/include -cpp -I$(CRTMPATH)/include

LIBRARIES =  -L /usr/lib -lnetcdff -lnetcdf -L/usr/lib/gcc/x86_64-redhat-linux/4.1.1  -lgfortran -lpng 

FFTROOT = /apps/libs/fftw/gnu/3.3.2
INCLUDE = $(FFTROOT)/include
LIB = $(FFTROOT)/lib
FLAG = -lfftw3


#module_fft_with_filtering.o : module_fft_with_filtering.f90
#	$(FC)  -c module_fft_with_filtering.f90 -I $(INCLUDE) -L $(LIB)  $(FLAG) -o $@


GBCSOBJECTS = GbcsKinds.o GbcsBaseTypes.o GbcsConstants.o GbcsStringUtil.o \
	GbcsDateTime.o GbcsPath.o GbcsPixel.o GbcsErrorHandler.o GbcsPDFLookup.o \
	GbcsConfFile.o GbcsLUTShift.o GbcsTypes.o GbcsAtmPhysics.o GbcsGeopotential.o GbcsProfileGenerator.o \
	GbcsTimeSpace.o GbcsSystemTools.o GbcsInterpolators.o GbcsForecastModel.o \
	GbcsMatrixOps.o GbcsPDFLoaders.o \
	GbcsScattPhysics.o GbcsCleanUp.o \
	GbcsPixelLoaders.o \
	module_type.v3.o \
	module_convert.o \
	module_functions.o \
	module_functions_parse.o \
	module_read_data_header.o \
	module_filtering_thresholds.o \
	module_filtering_thresholds_prt.o \
	module_average_per_scanline.o \
	module_average_per_orbite.v3.o \
	module_sensitivities.o \
        module_urict.o\
	module_radiance_uncertainties.v2.o \
	module_rescale.o \
	module_easy_FCDR.o\
	module_netcdf.o \
	NOAA_LoadAVHRRLevel1B.v3.o
	#NOAA_LoadAVHRRLevel1B_sortie_par_orbite_with_channel_3_4_5.and_filtering.FCDR.v3.o
 
OBJECTS = $(GBCSOBJECTS) prt_utils.o  check_37_chan_sortie_par_orbite_with_channel_3_4_5.FCDR.essai_other_output.v3.o 
#module_fft_with_filtering.o check_37_chan_sortie_par_orbite_with_channel_3_4_5.v12.o 

vpath %.f90 /home/users/mdesmons/gbcs/src/GbcsMod_Globals /home/users/mdesmons/gbcs/src/GbcsMod_SystemTools \
	/home/users/mdesmons/gbcs/src/GbcsMod_ErrorHandler /home/users/mdesmons/gbcs/src/GbcsMod_CleanUp \
	/home/users/mdesmons/gbcs/src/NOAA /home/users/mdesmons/gbcs/src/GbcsMod_Misc /home/users/mdesmons/gbcs/src/GbcsMod_Bayesian \
	GENERIC IASI_L1 /home/users/mdesmons/gbcs/src/GbcsMod_Profiles \
	/home/users/mdesmons/gbcs/src/GbcsMod_ProcessImagery /home/users/mdesmons/gbcs/src/GbcsMod_Physics \
	/home/users/mdesmons/gbcs/src/GbcsMod_Maths /home/users/mdesmons/gbcs/src/GbcsMod_Initialize \
	/home/users/mdesmons/gbcs/src/GbcsMod_DataReaders /home/users/mdesmons/gbcs/src/GbcsMod_RTM  \
	/home/users/mdesmons/gbcs/src/AVHRR/


all: $(OBJECTS)
	$(FC) -o check_37_chan_sortie_par_orbite_with_channel_3_4_5.FCDR.essai_other_output.exe $(OBJECTS)  \
	$(LIBRARIES) -I $(INCLUDE) -L $(LIB) $(FLAG)

.PHONY: clean
clean:
	rm -rf *.o *.mod *.exe

include /home/users/mdesmons/gbcs/src/make.rules

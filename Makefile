 OBJS = module_array_utilities.o \
	module_init_utilities.o \
	module_mars24.o \
	module_model_constants.o \
	module_m_time.o \
	module_nrutils.o \
	module_planet_utilities.o \
	module_ra_chapman.o \
	module_ra_mars_common.o \
	module_ra_kdm.o \
	module_ra_mars_uv.o \
	module_ra_mars_wbm.o \
	module_ra_utils.o \
	module_ra_valverde.o \
	module_read_soundings.o \
	module_setup_tpz.o \
	module_wrf_error.o \
	module_diag_common.o \
	module_radiation_driver.o \
	namelist.o


NATIVE_RWORDSIZE = 4
RWORDSIZE = $(NATIVE_RWORDSIZE)
IWORDSIZE = 4
DWORDSIZE = 8
LWORDSIZE = 4
MAX_DOMAINS = 1

arch=-DIWORDSIZE=$(IWORDSIZE) -DDWORDSIZE=$(DWORDSIZE) -DRWORDSIZE=$(RWORDSIZE) -DLWORDSIZE=$(LWORDSIZE)
defines=-DWRF_MARS -DWRF_PLANET -DMAX_DOMAINS_F=$(MAX_DOMAINS) -g -O0 

# NCINC=-I/nasa/netcdf/4.4.1.1_serial/include
# NCLIB=-L/nasa/netcdf/4.4.1.1_serial/lib

NCINC=-I/usr/include
NCLIB=-L/usr/lib
compile=-ffree-form -ffree-line-length-0 -cpp $(NCINC) $(defines) $(arch) 
link=$(NCLIB)
# compiler=/opt/intel/oneapi/compiler/latest/bin/ifort -diag-disable=10448
compiler=gfortran

AR      = ar

main: main.f90 $(OBJS)
	$(compiler) $(compile) $(link) $^ -o main

%.o : %.f90
	$(compiler) -c $(compile) $^

module_planet_utilities.o : module_wrf_error.o module_nrutils.o module_model_constants.o module_mars24.o

module_ra_mars_common.o : module_wrf_error.o module_mars24.o module_model_constants.o module_nrutils.o

module_ra_kdm.o : module_ra_utils.o module_ra_valverde.o module_ra_mars_common.o module_ra_mars_uv.o module_ra_mars_wbm.o module_setup_tpz.o module_array_utilities.o module_diag_common.o module_planet_utilities.o

module_ra_mars_uv.o : module_ra_chapman.o module_state_description.o module_model_constants.o module_wrf_error.o module_ra_mars_common.o

module_setup_tpz.o : module_init_utilities.o module_nrutils.o module_read_soundings.o

module_read_soundings.o : module_m_time.o module_wrf_error.o module_model_constants.o module_mars24.o module_nrutils.o

module_nrutils.o : module_array_utilities.o

module_ra_valverde.o : module_array_utilities.o

module_radiation_driver.o : module_wrf_error.o module_model_constants.o

module_diag_common.o : module_wrf_error.o

main.o : module_radiation_driver.o module_ra_kdm.o module_planet_utilities.o namelist.o

clean:
	rm -f *.o *.mod

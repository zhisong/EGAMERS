# EGAMERS Makefile

#****** Please choose the following options *******
# COMPILER SELECTIONS
# for Intel, use COMPILER = ifort
# for gfortran, use COMPILER = gfortran
# !! To use gfortran, one needs to install the GNU FORTRAN Compiler
# !! To use ifort, one needs to install the INTEL FORTRAN Compiler (license required)
# !! We haven't tested on other fortran compilers

# MKL
# compile with MKL, use MKL = yes
# compile with LAPACK, use MKL = no
# !! EGAMERS needs matrix operation subroutines to run
# !! One can use either LAPACK+BLAS(GNU) or MKL(INTEL), depending on the availability
# !! To use MKL, one needs to install the INTEL Math Kernal Library (license required)
# !! The code will fail on newest verion of MKL, please use <2016 version

# MPI
# compile without MPI (can use only one cpu, good for laptop), use MPI = no
# compile with MPI (can use many cpus and cluster), use MPI = yes
# !! Message Passing Interface (MPI) is a communications protocol for parallel computation
# !! MPI is highly recommended for PIC mode. Linear solver can work resonablely fast without MPI.
# !! To use MPI, one need to install a MPI library
# !! Have tested OpenMPI 1.10.2 and Intel MPI 2018. 

# NETCDF
# compile without NETCDF (output binary unformatted data), use NETCDF = no
# compile with serial NETCDF (output NETCDF format), use NETCDF = serial
# compile with parallel NETCDF, use NETCDF = parallel
# !! NETCDF is a stucture for output data, can be easily read by Python or Matlab
# !! NOTE: NETCDF/binary will only be used for PIC output, other modes will output text files

COMPILER = gfortran
MKL = yes
MPI = yes
DEBUG = no
NETCDF = serial
# The path to NETCDF, please adjust
NETCDF_INCLUDE = /usr/include/
#NETCDF_INCLUDE = $(NETCDFHOME)/include
NETCDF_LIB = /usr/lib/x86_64-linux-gnu/
#NETCDF_LIB = $(NETCDFHOME)/lib

#****************************************************

# the name of the executable
TARGET = ../egamers.exe

# the source files
OBJS =  paras_phy.o paras_num.o interfaces.o matrix_module.o sgtsl.o cubic.o spline_module.o \
        profile.o distribution_fun.o sintable.o hermite.o landau_integral.o radial_grid.o \
        orbit_classify.o orbit.o orbit_integral.o mpi.o random.o trap_grid.o trap_matrix.o \
        mhd_matrix.o field.o eigen.o continuum.o particles.o test_particles.o pic.o \
        diagnostics.o nl.o io.o egamers.o

ifeq ($(COMPILER), ifort)
  FC = mpif90
  FFLAGS = -r8 -cpp
  ifeq ($(DEBUG), yes)
    FFLAGS += -O0 -g
  else
    FFLAGS += -O2
  endif
endif

ifeq ($(COMPILER), gfortran)
  FC = mpif90
  FFLAGS = -fdefault-real-8 -cpp
  ifeq ($(DEBUG), yes)
    FFLAGS += -O0 -g
  else
    FFLAGS += -O3
  endif

endif

ifeq ($(MKL), yes)
  LIB = -mkl
  MACRO += -DMKL
else
  LIB = -llapack -lblas
endif

ifeq ($(MPI), yes)
  MACRO += -DMPI
endif

ifeq ($(NETCDF), serial)
  MACRO += -DNC
  LIB += -lnetcdff -L$(NETCDF_LIB)
  INCLUDE += -I$(NETCDF_INCLUDE)
endif
ifeq ($(NETCDF), parallel)
  MACRO += -DNC -DNCPAR
  LIB += -lnetcdff -L$(NETCDF_LIB)
  INCLUDE += -I$(NETCDF_INCLUDE)
endif


all :	$(OBJS)
	$(FC) ${FFLAGS} ${MACRO} -o $(TARGET) $(MODULES) $(OBJS) $(LIB)

clean :
	rm *.o
	rm *.mod
	rm $(TARGET)

%.o : %.f90
	$(FC) ${FFLAGS} ${MACRO} ${INCLUDE} -c $<

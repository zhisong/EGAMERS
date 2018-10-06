# COMPILER SELECTIONS
# for Intel, use COMPILER = ifort
# for gfortran, use COMPILER = gfortran
COMPILER = gfortran

# MKL
# compile with MKL, use MKL = yes
# compile with LAPACK, use MKL = no
MKL = no

# MPI
# compile without MPI, use MPI = no
# compile with MPI, use MPI = yes
MPI = yes

# the name of the executable
TARGET = ../egamers.exe

# the source files
OBJS =  paras_phy.o paras_num.o interfaces.o matrix_module.o sgtsl.o cubic.o spline_module.o profile.o distribution_fun.o sintable.o hermite.o landau_integral.o radial_grid.o orbit_classify.o orbit.o orbit_integral.o mpi.o random.o trap_grid.o trap_matrix.o mhd_matrix.o field.o eigen.o particles.o pic.o nl.o io.o egamers.o

ifeq ($(COMPILER), ifort)
  FC = mpiifort
  FFLAGS = -r8 -O2 -cpp
endif

ifeq ($(COMPILER), gfortran)
  FC = mpif90
  FFLAGS = -fdefault-real-8 -O3 -cpp
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

all :	$(OBJS)
	$(FC) ${FFLAGS} ${MACRO} -o $(TARGET) $(MODULES) $(OBJS) $(LIB)

clean :
	rm *.o
	rm *.mod
	rm $(TARGET)

%.o : %.f90
	$(FC) ${FFLAGS} ${MACRO} -c $<

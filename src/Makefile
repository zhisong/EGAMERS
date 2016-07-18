FC = ifort
LIB = -lblas -llapack
TARGET = ../egamers.exe
OBJS =  paras_phy.o paras_num.o interfaces.o matrix_module.o sgtsl.o cubic.o spline_module.o profile.o distribution_fun.o sintable.o hermite.o landau_integral.o radial_grid.o orbit.o orbit_integral.o orbit_classify.o trap_grid.o trap_matrix.o mhd_matrix.o eigen.o io.o egamers.o
FFLAGS = -r8

all :   $(OBJS)
	$(FC) ${FFLAGS} -o $(TARGET) $(MODULES) $(OBJS) $(LIB)

clean :
	rm *.o
	rm *.mod
	rm $(TARGET)

%.o : %.f90
	$(FC) ${FFLAGS} -c $<

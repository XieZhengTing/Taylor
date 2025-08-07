FC              =    nvfortran
OPLEVEL         =   -O3
LDFLAGS := -lnetcdff -lexoIIv2for32 -lexodus -lnetcdf
EXDIR           =   /usr/local
EXDIRLIB        =   $(EXDIR)/lib
EXDIRINC        =   $(EXDIR)/include
NETCDFDIR       =   /usr/local
NETCDFDIRLIB    =   $(NETCDFDIR)/lib
NETCDFDIRINC    =   $(NETCDFDIR)/include


FFLAGS  =   $(OPLEVEL)  $(EXSWITCH) -L $(EXDIRLIB) -I $(EXDIRINC) -I $(NETCDFDIRINC) -L $(NETCDFDIRLIB) -acc -gpu=cc75 -Minfo=accel -mp
SRCS1=$(wildcard *.f)
SRCS2=$(wildcard *.F)
SRCS3=$(wildcard *.F90)
SRCS4=$(wildcard *.f90)
OBJS1=$(SRCS1:.f=.o)
OBJS2=$(SRCS2:.F=.o)
OBJS3=$(SRCS3:.F90=.o)
OBJS4=$(SRCS4:.f90=.o)
OBJS1 := $(filter-out FINT_FUNCTIONS.o MODEL.o CONTROL.o GPU_ERROR.o, $(OBJS1))
OBJECTS = FINT_FUNCTIONS.o MODEL.o CONTROL.o GPU_ERROR.o $(OBJS4) $(OBJS3) $(OBJS2) $(OBJS1)
OBJECTS2 = $(OBJS4) $(OBJS3) $(OBJS2) $(OBJS1)
EXECUTABLE = MEGA
all:    $(EXECUTABLE)
CONTROL.o CONTROL.mod : CONTROL.F90
	$(FC) $(FFLAGS) -c -o $@ $<

FINT_FUNCTIONS.o FINT_FUNCTIONS.mod : FINT_FUNCTIONS.F90
	$(FC) $(FFLAGS) -c -o $@ $<

MODEL.o MODEL.mod : MODEL.F90
	$(FC) $(FFLAGS) -c -o $@ $<
GPU_ERROR.o GPU_ERROR.mod : GPU_ERROR.F90
	$(FC) $(FFLAGS) -c -o $@ $<
%.o: %.f90 FINT_FUNCTIONS.mod GPU_ERROR.mod
	$(FC) $(FFLAGS) -c $< -o $@
%.o: %.F90 FINT_FUNCTIONS.mod GPU_ERROR.mod
	$(FC) $(FFLAGS) -c $< -o $@
%.o: %.F FINT_FUNCTIONS.mod
	$(FC) $(FFLAGS)	-c $< -o $@
%.o: %.f FINT_FUNCTIONS.mod
	$(FC) $(FFLAGS)	-c $< -o $@
$(EXECUTABLE): $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS2)  -o $(EXECUTABLE) $(LDFLAGS)
clean:
	$(RM) $(EXECUTABLE)
	$(RM) $(OBJECTS)
	$(RM) control.mod fint_functions.mod model.mod gpu_error.mod
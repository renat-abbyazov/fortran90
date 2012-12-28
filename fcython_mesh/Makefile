# Modify these for your Fortran compiler:

# GFortran
F90 = gfortran
F90FLAGS = -Wall -std=f2008 -Wextra -Wimplicit-interface -fPIC
# Debug flags:
F90FLAGS += -g -fcheck=all
# Release flags:
#F90FLAGS += -O3 -march=native -ffast-math -funroll-loops

# Intel ifort
#F90 = ifort-12.0.191
#F90FLAGS = -stand f95 -warn all
# Debug flags:
#F90FLAGS += -check all
# Release flags:
#F90FLAGS += -xHOST -O3 -no-prec-div -static

#GCC
CC = gcc
CFLAGS = -Wall -fPIC -I$(SPKG_LOCAL)/include/python2.6 -I$(SPKG_LOCAL)/lib/python2.6/site-packages/numpy/core/include

fcython.so: fmesh_wrapper.o fcython.o fmesh.o utils.o types.o
	$(F90) -shared -o fcython.so $^ -lpython2.6

fcython.c: fcython.pyx
	cython fcython.pyx

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

%.o: %.c
	$(CC) $(CFLAGS) -c $<

fmesh_wrapper.o: fmesh.o
fmesh.o: types.o utils.o
utils.o: types.o

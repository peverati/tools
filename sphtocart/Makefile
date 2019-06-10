#! /usr/bin/make

obj = sph2real
all : $(obj)

FCOMPL = gfortran
FFLAGC = -O2
#-Wpedantic -g -pg -Wunused -fbacktrace -fcheck=bounds,mem,pointer,do,array-temps -Wall 
CCOMPL = gcc
CFLAGC =
opt    = -O2
SYSDEP =
ldflag =

.SUFFIXES:
.SUFFIXES: .f90 .o .F90

.F90.o:
	$(FCOMPL) -c $(FFLAGC) $(opt) -o $@ $<

.f90.o:
	$(FCOMPL) -c $(FFLAGC) $(opt) -o $@ $<

OBJECTS = mod_io.o mod_param.o sph2real.o

# Clean objects

clean:
	@rm -f core $(SYSDEP) $(VECTOR) $(OBJECTS) *.mod

veryclean:
	@rm -f core $(SYSDEP) $(VECTOR) $(OBJECTS) $(obj) *.mod

$(obj): $(SYSDEP) $(VECTOR) $(OBJECTS)
	$(FCOMPL) $(FFLAGC) -o $(obj) $(OBJECTS) $(SYSDEP) $(ldflag)


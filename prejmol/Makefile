#! /usr/bin/make

obj = prejmol
all : $(obj)

FCOMPL = gfortran
FFLAGC = 
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

OBJECTS = mod_io.o mod_param.o mod_wfn.o prejmol.o

# Clean objects

clean:
	@rm -f core $(SYSDEP) $(VECTOR) $(OBJECTS) *.mod

veryclean:
	@rm -f core $(SYSDEP) $(VECTOR) $(OBJECTS) *.mod $(obj)

$(obj): $(SYSDEP) $(VECTOR) $(OBJECTS)
	$(FCOMPL) $(FFLAGC) -o $(obj) $(OBJECTS) $(SYSDEP) $(ldflag)

prejmol.o : mod_param.mod
prejmol.o : mod_wfn.mod
prejmol.o : mod_io.mod
mod_wfn.o : mod_io.mod

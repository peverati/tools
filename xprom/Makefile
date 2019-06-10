#! /usr/bin/make

obj = xprom
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

OBJECTS = mod_io.o mod_datatm.o mod_xprom.o xprom.o \
					pair.o

# Clean objects

clean:
	@rm -f core $(SYSDEP) $(VECTOR) $(OBJECTS) *.mod

veryclean:
	@rm -f core $(SYSDEP) $(VECTOR) $(OBJECTS) $(obj) *.mod

$(obj): $(SYSDEP) $(VECTOR) $(OBJECTS)
	$(FCOMPL) $(FFLAGC) -o $(obj) $(OBJECTS) $(SYSDEP) $(ldflag)

mod_xprom.o mod_xprom.mod : mod_datatm.mod
xprom.o xprom.mod pair.o pair.mod : mod_xprom.mod
xprom.o xprom.mod pair.o pair.mod mod_xprom.o mod_xprom.mod : mod_io.mod

#! /usr/bin/make

include make.inc

all: integ
exe = integ
obj = mod_prec.o mod_io.o mod_memory.o mod_param.o \
			mod_futils.o mod_wfn.o mod_gto.o integ.o 
BLAS = -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lm -lpthread 

.SUFFIXES: .o .f90

%.o: %.f90
	$(FCOMPL) -c $(FFLAGC) $(OPT) $(FDEBUG) -o $@ $<

%.o: %.F90
	$(FCOMPL) -c $(FFLAGC) $(OPT) $(FDEBUG) -o $@ $<

integ: $(obj)
	$(FCOMPL) $(LDFLAG) -o $(exe) $(obj) libcint.a $(BLAS)

clean:
	@rm -f $(obj) *.mod $(exe)

.PHONY: all integ clean

# DEPS

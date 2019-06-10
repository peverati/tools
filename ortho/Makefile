#! /usr/bin/make

include make.inc

all: ortho
exe = ortho
obj = mod_prec.o mod_io.o mod_memory.o mod_param.o \
			mod_futils.o mod_wfn.o mod_gto.o ortho.o 

.SUFFIXES: .o .f90

%.o: %.f90
	$(FCOMPL) -c $(FFLAGC) $(OPT) $(FDEBUG) -o $@ $<

%.o: %.F90
	$(FCOMPL) -c $(FFLAGC) $(OPT) $(FDEBUG) -o $@ $<

ortho: $(obj)
	$(FCOMPL) $(LDFLAG) -o $(exe) $(obj) libcint.a -lblas -llapack

clean:
	@rm -f $(obj) *.mod $(exe)

.PHONY: all ortho clean

# DEPS
mod_gto.o mod_gto.mod : mod_futils.mod
ortho.o ortho.mod : mod_gto.mod
mod_wfn.o mod_wfn.mod mod_gto.o mod_gto.mod : mod_memory.mod
ortho.o ortho.mod mod_gto.o mod_gto.mod : mod_param.mod
ortho.o ortho.mod mod_wfn.o mod_wfn.mod mod_param.o mod_param.mod mod_memory.o mod_memory.mod mod_gto.o mod_gto.mod mod_futils.o mod_futils.mod : mod_prec.mod
ortho.o ortho.mod mod_gto.o mod_gto.mod : mod_wfn.mod
ortho.o ortho.mod mod_wfn.o mod_wfn.mod mod_memory.o mod_memory.mod mod_gto.o mod_gto.mod mod_futils.o mod_futils.mod : mod_io.mod

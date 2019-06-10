#! /usr/bin/make

include make.inc

all: evalwfn
exe = evalwfn
obj = mod_prec.o mod_io.o mod_memory.o mod_param.o mod_wfn.o mod_fields.o mod_plot.o evalwfn.o

.SUFFIXES: .o .f90

%.o: %.f90
	$(FCOMPL) -c $(FFLAGC) $(OPT) $(FDEBUG) -o $@ $<

%.o: %.F90
	$(FCOMPL) -c $(FFLAGC) $(OPT) $(FDEBUG) -o $@ $<

evalwfn: $(obj)
	$(FCOMPL) $(LDFLAG) -o $(exe) $(obj)

clean:
	@rm -f $(obj) *.mod $(exe)

.PHONY: all evalwfn clean

# DEPS
mod_plot.o mod_plot.mod : mod_fields.mod
mod_wfn.o mod_wfn.mod mod_plot.o mod_plot.mod : mod_memory.mod
mod_wfn.o mod_wfn.mod mod_plot.o mod_plot.mod evalwfn.o evalwfn.mod : mod_param.mod
evalwfn.o evalwfn.mod : mod_plot.mod
mod_wfn.o mod_wfn.mod mod_plot.o mod_plot.mod mod_param.o mod_param.mod mod_memory.o mod_memory.mod mod_fields.o mod_fields.mod evalwfn.o evalwfn.mod : mod_prec.mod
mod_plot.o mod_plot.mod mod_fields.o mod_fields.mod evalwfn.o evalwfn.mod : mod_wfn.mod
mod_wfn.o mod_wfn.mod mod_plot.o mod_plot.mod mod_memory.o mod_memory.mod evalwfn.o evalwfn.mod : mod_io.mod


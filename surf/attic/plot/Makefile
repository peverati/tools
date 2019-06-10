#! /usr/bin/make

include make.inc

all: plot
exe = plot
obj = mod_prec.o mod_io.o mod_memory.o mod_param.o mod_plot.o \
      plot.o

.SUFFIXES: .o .f90

%.o: %.f90
	$(FCOMPL) -c $(FFLAGC) $(OPT) $(FDEBUG) -o $@ $<

%.o: %.f
	$(FCOMPL) -c $(FFLAGC) $(OPT) $(FDEBUG) -o $@ $<

%.o: %.F90
	$(FCOMPL) -c $(FFLAGC) $(OPT) $(FDEBUG) -o $@ $<

plot: $(obj)
	$(FCOMPL) $(LDFLAG) -o $(exe) $(obj)

clean:
	@rm -f $(obj) *.mod $(exe)

.PHONY: all plot clean

# DEPS

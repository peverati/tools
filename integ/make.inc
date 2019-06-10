FCOMPL = gfortran
#FDEBUG = -Wpedantic -g -pg -Wunused -fbacktrace -fcheck=bounds,mem,pointer,do,array-temps -Wall
FFLAGC = -m64 -pipe -falign-commons -fopenmp
CCOMPL = gcc
CFLAGC = -O2 -m64 -mtune=native -Wunused
OPT    = -O3 -ftree-loop-optimize \
-ftree-vectorize -ftree-loop-distribution -mtune=native -march=native
LDFLAG = -m64 -fopenmp
AR = ar rc
RANLIB = ranlib

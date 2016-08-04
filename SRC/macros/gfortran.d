# ------------------------------------------------------ #
#	architecture Makefile programme CHE (cold-hot-plot)
# ------------------------------------------------------ #

COMPIL = gfortran
COMPILMPI = mpif90 
OPTIONC = -O3 -cpp
OPTIOND = -fno-range-check # pb in mt19937ar avec gfortran

FFLAGS = -fbounds-check

# -------------------- Options ------------------------- #

# test :

#FFLAGS = -fdefault-real-8 -O0 -g -fbounds-check -Wall -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -ftrapv -fimplicit-none -fno-automatic -ffree-form -Wconversion -Wunderflow -Wunreachable-code -g3 -fstack-protector-all

# ------------------------------------------------------ #
# cpp -> permet de définir _FILE_DIR_ en fonction du compilateur
# ------------------------------------------------------ #



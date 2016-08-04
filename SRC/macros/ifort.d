# ------------------------------------------------------ #
#	architecture Makefile programme CHE (cold-hot-plot)
# ------------------------------------------------------ #


# mpimodfile=-I/U.../MPICH2/mpich-3.1.2/mpich-install/include/ ###  sudo find . -name mpi.mod

COMPIL = ifort
COMPILMPI = mpif90 -f90=ifort
OPTIONC = -O3 -fpp -diag-disable 8291 -diag-disable 8290 -assume byterec

# FFLAGS = -check bounds  

# ------------------------------------------------------ #
# fpp -> permet de définir _FILE_DIR_ en fonction du compilateur
# ------------------------------------------------------ #
# -diag-disable 8291 -diag-disable 8290  -> format d'écriture dans des fichiers des variables output
# remark #8291: Recommended relationship between field width 'W' and the number of fractional digits 'D' in this edit descriptor is 'W>=D+7'
# remark #8290: Recommended relationship between field width 'W' and the number of fractional digits 'D' in this edit descriptor is 'W>=D+3'
# ------------------------------------------------------ #
# assume byterec -> force to use byte units in access='direct' file
# ------------------------------------------------------ #

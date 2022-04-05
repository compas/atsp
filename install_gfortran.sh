#!/bin/bash
#########################################################################

export FC="gfortran"                       # Fortran compiler
export FC_FLAGS="-O2 -fno-automatic -fcray-pointer"
export FC_LD="-static-libgfortran"         # Serial linker flags
export ATSP="${PWD}"
export LAPACK_LIB="-llapack -lblas"
export LAPACK_DIR="/usr/lib"
export FC_MALLOC="LINUX"

#########################################################################

export FC_MPI="mpif90"                     # MPI compiler
export FC_MPIFLAGS="-O2 -fno-automatic -fcray-pointer"       # Parallel code compiler flags
export FC_MPILD="-O"                         # Parallel linker flags
export MPI_TMP="/home/nemouchi/tmp_mpi"

#########################################################################

export CPP="g++"                           # C++ compiler
export CPP_FLAGS="-O3"                     # C++ compiler flags
export CPP_LD="-static"                    # C++ linker

#########################################################################

cd ${ATSP}/src
make $* 2>&1 | tee ${ATSP}/make.err
grep "severe" ../make.err

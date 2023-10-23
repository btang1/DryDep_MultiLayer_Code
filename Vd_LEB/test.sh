#!/bin/bash

set -ex

rm -f *.mod *.o *.exe

# Compile modules
gfortran -Wall -c ACCESS_Constants.f90
gfortran -Wall -c ACCESS_Modules.f90
gfortran -Wall -c DryDep.f90

gfortran -Wall -o test_get_coeffs.exe test_get_coeffs.f90 DryDep.o ACCESS_Modules.o ACCESS_Constants.o

./test_get_coeffs.exe

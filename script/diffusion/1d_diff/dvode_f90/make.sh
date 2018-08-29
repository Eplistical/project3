#!/bin/bash

gfortran -c -fPIC *.f90
gfortran -shared -fPIC -o libdvode_f90.so *.o
rm *.o

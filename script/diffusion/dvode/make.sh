#!/bin/bash

gfortran -c -fPIC *.f90
gfortran -shared -fPIC -o libdvode.so *.o
rm *.o

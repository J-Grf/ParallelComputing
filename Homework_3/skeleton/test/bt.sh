#!/bin/bash

cd ../build
module load DEVELOP; module load intel/19.0; export CXX=`which icpc`; export C=`which icc`; make;

cd ../test
./2d_Unsteady_OpenMP settings.coarse.in
./2d_Unsteady_OpenMP settings.coarse.in
./2d_Unsteady_OpenMP settings.coarse.in

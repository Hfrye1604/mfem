#!/bin/sh

#g++  -O3 -std=c++11 -I../../  NSnonlininteg_modified.cpp -L../../ -lmfem -lrt
g++  -g -O3 -std=c++11 -I../../ perturbed_drift_diffusion_1d.cpp -o pdd1d -L../../ -lmfem 


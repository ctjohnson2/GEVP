#!/bin/bash



g++ -O2 -fopenmp -I/home/chris/code/hadspec/scattering_devel/itpp/itpp-4.2/ -c -o fitter.o fitter.cc -L/home/chris/code/hadspec/install/itpp/lib -litpp -llapack -lcblas -lf77blas -latlas -lgsl

g++ -O2 -fopenmp -I/home/chris/code/hadspec/scattering_devel/itpp/itpp-4.2/ -o fitter fitter.o -L/home/chris/code/hadspec/install/itpp/lib -litpp -llapack -lcblas -lf77blas -latlas -lgsl

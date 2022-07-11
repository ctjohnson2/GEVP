#!/bin/zsh

set itpp_path="-I/Users/christopherjohnson/code/hadspec/scattering_devel/itpp/itpp-4.2"
set itpp_lib="/Users/christopherjohnson/install/itpp/lib"

set LD_FLAGS="-L$itpp_lib -litpp -llapack -lcblas -lf77blas -latlas -lgsl"
CXX="/opt/homebrew/Cellar/gcc/11.2.0_3/bin/g++-11"
CXX_FLAGS="-fopenmp -O2"
#g++ -O2 -fopenmp -I/$itpp_path -o fitter.o -c fitter.cc -L$itpp_lib -litpp -llapack -lcblas -lf77blas -latlas -lgsl

$CXX -O2 -fopenmp -I/opt/homebrew/Cellar/gsl/2.7.1/include -I/Users/christopherjohnson/code/hadspec/scattering_devel/itpp/itpp-4.2 -L/Users/christopherjohnson/install/itpp/lib -litpp -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgsl -o fitter fitter.cc

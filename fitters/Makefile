## simple make file, make sure you change paths for compilers and libraries
INSTALL_DIR=$(HOME)/projects/GEVP/fitters/
CXX="/opt/homebrew/Cellar/gcc/11.2.0_3/bin/g++-11"
CXXFLAGS= -fopenmp -O3 -I/opt/homebrew/Cellar/gsl/2.7.1/include -I/Users/christopherjohnson/code/hadspec/scattering_devel/itpp/itpp-4.2 \

LDFLAGS=-L/Users/christopherjohnson/install/itpp/lib -L/opt/homebrew/Cellar/gsl/2.7.1/lib
LIBS=-litpp -lgsl

SRCS = fitter.cc 

OBJS = $(SRCS:.cc=.o)

fitter: $(OBJS)
	$(CXX)  $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ $(SRCS)


install: fitter
	mkdir -p $(INSTALL_DIR)
	cp fitter $(INSTALL_DIR)

uninstall:
	rm -f $(INSTALL_DIR)/fitter
	
clean:
	rm -f *.o; rm -f fitter;

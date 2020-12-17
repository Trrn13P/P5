#CPPflags=  clang++ -Xpreprocessor -fopenmp -std=c++14
LIB = -larmadillo -llapack -lblas
CPPflags= g++ -std=c++14

# n, tsteps, saved_tsteps, alpha, type
#make sure alpha is a decimal point number, 0.25, not 1/4 for ex.
CONDITIONS = 10 100 5 0.25 1d

all: compile execute

compile:
	${CPPflags} main.cpp diffusion.cpp tridiag_solver.cpp -o ./main.out ${LIB}

execute:
	./main.out ${CONDITIONS}

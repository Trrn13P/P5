#CPPflags=  clang++ -Xpreprocessor -fopenmp -std=c++14
LIB = -larmadillo -llapack -lblas
CPPflags= g++ -std=c++14


all: compile execute

compile:
	${CPPflags} main.cpp diffusion.cpp tridiag_solver.cpp -o ./main.out ${LIB}

execute:
	./main.out

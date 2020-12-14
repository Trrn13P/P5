#include "diffusion.hpp"

#include <iostream>
#include <fstream>
using namespace std;

  int main(int argc, char const *argv[]) {
    int n = 100;
    int tsteps = 10;
    float dx = 1./(n+1);
    float alpha = 1./4;
    float dt = alpha*dx*dx;


    diffusion *test;
    test = new diffusion(n, tsteps, dx, dt);
    test -> backward_euler();
    delete test;


    /*
    ofstream outfile("./test.txt");

    float a,b,c;
    a=b=c=1;

    tridiag_solver *test;
    test = new tridiag_solver(n,a,b,c);
    for(int t=0;t<tsteps;t++){
      //test -> forward_solver(1-2*alpha,alpha);
      //test -> backward_solver();
      //test -> forward_euler();
    }

    //test -> crank_nicolson(n,tsteps,delta_x,alpha);

    outfile.close();
    delete test;

    test = new matrix_solver;
    test -> backward_euler(n,tsteps,delta_x,alpha);
    delete test;

    test = new matrix_solver;
    test -> forward_euler(n,tsteps,delta_x,alpha);
    delete test;
    return 0;
    */

  }

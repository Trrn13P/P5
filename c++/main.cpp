#include "diffusion.hpp"

#include <iostream>
#include <fstream>
using namespace std;

  int main(int argc, char const *argv[]) {
    //Setting up parameters
    int n = 10;
    int tsteps = 1000;
    float dx = 1./(n+1);
    float alpha = 1./4;
    float dt = alpha*dx*dx;

    //Solving equation with set parameters
    diffusion *test;
    test = new diffusion(n, tsteps, dx, dt);
    test -> backward_euler();
    delete test;
  }

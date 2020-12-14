#include "diffusion.hpp"

#include <iostream>
#include <fstream>
using namespace std;

  int main(int argc, char const *argv[]) {
    int n = 10;
    int tsteps = 100;
    float dx = 1./(n+1);
    float alpha = 1./4;
    float dt = alpha*dx*dx;


    diffusion *test;
    test = new diffusion(n, tsteps, dx, dt);
    test -> backward_euler();
    delete test;
  }

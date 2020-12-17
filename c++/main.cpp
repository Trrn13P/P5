#include "diffusion.hpp"

#include <iostream>
#include <fstream>
//#include <string>
using namespace std;

  int main(int argc, char const *argv[]) {
    if (argc <= 5) {
      cout << "Bad Usage: " << argv[0] << endl;
      exit(1);
    }
    int n, tsteps, saved_tsteps;
    float alpha;
    string type;
    if (argc > 1) {
      //Setting up parameters from makefile
      n = atoi(argv[1]);
      tsteps = atoi(argv[2]);
      //number of saved timesteps, 0 for saving all. 1 for saving every second
      saved_tsteps = atoi(argv[3]);
      alpha = atof(argv[4]);
      //type can be 1d or 2d
      type = argv[5];
    }

    float dx = 1./(n+1);
    float dt = alpha*dx*dx;

    //Solving equation with set parameters
    diffusion *solve;
    solve = new diffusion(n ,tsteps ,saved_tsteps, dx, dt);
    if(type=="1d"){
      solve -> backward_euler();
      solve -> forward_euler();
      solve -> crank_nicolson();
    }
    if(type=="2d"){
      solve -> forward_2d();
    }
    delete solve;
  }

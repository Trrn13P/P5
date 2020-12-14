#include <armadillo>
using namespace arma;

#include <fstream>
using namespace std;

#ifndef DIFFUSION_HPP
#define DIFFUSION_HPP

class diffusion {
private:
  //Setting up ints, floats
  float a ,b ,c;
  int n, tsteps;
  float dx, dt, alpha;

  //Setting up vectors
  vec u, unew;

  void Initialize(int n_, int tsteps_, float dx_, float dt_){
    //Initialing stepsize, numbers ect.
    dx = dx_; dt = dt_;
    alpha = dt*1./dx*1./dx;

    n = n_;
    tsteps = tsteps_;

    //creating zero array for u and unew
    u = unew = zeros(n+1);
  }

public:
  float func(float x);
  void backward_euler();
  void forward_euler();
  void crank_nicolson();

  void writetofile(ofstream &outfile ,vec u_);

  //Initializing
  diffusion(int n_, int tsteps_, float dx_, float dt_){
    Initialize(n_, tsteps_, dx_, dt_);
  }
};
#endif

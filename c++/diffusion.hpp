#include <armadillo>
using namespace arma;

#ifndef DIFFUSION_HPP
#define DIFFUSION_HPP

class diffusion {
private:
  float a ,b ,c;
  //Setting up integers
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

    u = unew = zeros(n+1);

    //Initializing x
    //x_vec = linspace(0,1,n+1);
  }


public:
  float func(float x);
  void backward_euler();

  //Initializing
  diffusion(int n_, int tsteps_, float dx_, float dt_){
    Initialize(n_, tsteps_, dx_, dt_);
  }
};
#endif

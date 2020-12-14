#include <armadillo>
using namespace arma;

#ifndef TRIDIAG_SOLVER_HPP
#define TRIDIAG_SOLVER_HPP

class tridiag_solver {
  private:
    vec b_tilde ,a_tilde , g_tilde;

    float a, b;


    int n;
    void Initialize(int n_, float a_, float b_){
      n = n_;
      a = a_;
      b = b_;

      //setting up the tilde vectors
      b_tilde = a_tilde = g_tilde = zeros(n+1);
    }

public:
  void forward_solver(vec &v_vec, vec &g_vec);
  void backward_solver(vec &v_vec, vec &g_vec);

  //Initializing
  tridiag_solver(int n_, float a_, float b_){
    Initialize(n_, a_, b_);
  }
};

#endif

#include <armadillo>
using namespace arma;

#ifndef TRIDIAG_SOLVER_HPP
#define TRIDIAG_SOLVER_HPP

class tridiag_solver {
  private:
    vec b_tilde ,c_tilde , g_tilde;

    float a, b, c;


    int n;
    void Initialize(int n_, float a_, float b_, float c_){
      n = n_;
      a = a_;
      b = b_;
      c = c_;

      //setting up the tilde vectors
      b_tilde = c_tilde = g_tilde = zeros(n);
    }

public:
  void forward_solver(vec &v_vec, vec &g_vec);
  void backward_solver(vec &v_vec, vec &g_vec);

  //Initializing
  tridiag_solver(int n_, float a_, float b_, float c_){
    Initialize(n_, a_, b_, c_);
  }
};

#endif

#include "tridiag_solver.hpp"

void tridiag_solver::forward_solver(vec &v_vec, vec &g_vec){
  //Setting inital values for the tilde vectors
  b_tilde(0) = b; c_tilde(0) = c; g_tilde(0) = g_vec(0);

  //Solving the forward algorythm
  for(int i=1;i<n;i++){
    c_tilde(i) = b_tilde(i-1);
    b_tilde(i) = -2*b_tilde(i-1)-c_tilde(i-1);
    g_tilde(i) = -g_vec(i)*b_tilde(i-1)-g_tilde(i-1);
  }
}

void tridiag_solver::backward_solver(vec &v_vec, vec &g_vec){
  //Setting up the endpoint
  v_vec(n) = g_tilde(n-1)*1./b_tilde(n-1);

  //Solving the backward algorythm
  int j;
  for(int i=2;i<n+1;i++){
    j = (n-i);
    v_vec(j+1) = (g_tilde(j) - c_tilde(j) *v_vec(j+2))*1./b_tilde(j);
  }
}

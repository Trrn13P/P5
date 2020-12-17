#include "tridiag_solver.hpp"

/*
THIS FORWARD SOLVER IS SPECIALIZED FOR EQUAL DIAGONAL ELEMENTS a=c
*/

void tridiag_solver::forward_solver(vec &v_vec, vec &g_vec){
  //Setting inital values for the tilde vectors
  b_tilde(0) = b; a_tilde(0) = a; g_tilde(0) = g_vec(0);


  //precalculating b/a
  float b_over_a = b*1./a;
  //Solving the forward algorythm, specialized
  for(int i=1;i<n+1;i++){
    a_tilde(i) = b_tilde(i-1);
    b_tilde(i) = b_over_a*b_tilde(i-1)-a_tilde(i-1);
    g_tilde(i) = g_vec(i)*b_tilde(i-1)*1./a    - g_tilde(i-1);
  }
}

void tridiag_solver::backward_solver(vec &v_vec, vec &g_vec){
  //Setting up the endpoint
  v_vec(n) = g_tilde(n)*1./b_tilde(n);

  //Backward algorythm
  for(int j=n-1;j>=0;j--){
    v_vec(j) = (g_tilde(j) - a_tilde(j) *v_vec(j+1))*1./b_tilde(j);
  }
}

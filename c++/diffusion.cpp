#include "tridiag_solver.hpp"
#include "diffusion.hpp"

#include <fstream>
using namespace std;

//one line pr timestep
void diffusion::writetofile(ofstream &outfile, vec u){
  outfile << u.t();
}

//g(x) boundary function
float diffusion::func(float x){
  return 1;
}


void diffusion::forward_euler(){
  ofstream outfile("../textfiles/forward_euler.txt");
  a = alpha;
  b = (1-2*alpha);

  //Boundary conditions g(x)
  for (int i = 1; i < n+1; i++) {
    u(i) = func(dx*(i));
  }
  // Boundary conditions a(t) b(t) (zero here)
  unew(n) = u(n) = u(0) = unew(0) = 0;

  //writing initial conditions
  writetofile(outfile,u);

  //going trough timesteps
  for (int t = 1; t <= tsteps; t++) {
    for (int i = 1; i < n; i++) {
    // Discretized diff eq
      unew(i) = alpha * u(i-1) + (1 - 2*alpha) * u(i) + alpha * u(i+1);
    }
    //setting boundary and updating u
    unew(0) = unew(n) = 0;
    u = unew;

    //print statement
    writetofile(outfile,u);
  }
  outfile.close();
}





void diffusion::backward_euler(){
  ofstream outfile("../textfiles/backward_euler.txt");
  a = -alpha;
  b = (1+2*alpha);

  //Boundary conditions g(x)
  for (int i = 1; i < n+1; i++) {
    u(i) = func(dx*(i));
  }
  // Boundary conditions a(t) b(t) (zero here)
  unew(n) = u(n) = u(0) = unew(0) = 0;

  //Setting up class pointer
  tridiag_solver *solver;
  solver = new tridiag_solver(n,a,b);

  //writing initial conditions
  writetofile(outfile,u);

  //going trough timesteps
  for (int t = 1; t <= tsteps; t++) {
    //forward-backward solve
    solver->forward_solver(unew,u);
    solver->backward_solver(unew,u);

    //setting boundary and updating u
    unew(0) = unew(n) = 0;
    u = unew;

    //print statement
    writetofile(outfile,u);
  }
  outfile.close();
  delete solver;
}

void diffusion::crank_nicolson(){
  ofstream outfile("../textfiles/crank_nicolson.txt");
  a = - alpha;
  b = 2 + 2*alpha;

  vec r(n+1); // Right side of matrix equation Au=r


  //Boundary conditions g(x)
  for (int i = 1; i < n+1; i++) {
    u(i) = func(dx*(i));
  }
  // Boundary conditions a(t) b(t) (zero here)
  unew(n) = u(n) = u(0) = unew(0) = 0;

  //writing initial conditions
  writetofile(outfile,u);

  //Setting up class pointer
  tridiag_solver *solver;
  solver = new tridiag_solver(n,a,b);

  // Time iteration
  for (int t = 1; t <= tsteps; t++) {
    // Calculate r for use in tridag, right hand side of the Crank Nicolson method
    for (int i = 1; i < n; i++) {
      r(i) = alpha*u(i-1) + (2 - 2*alpha)*u(i) + alpha*u(i+1);
    }
    r(0) = 0;
    r(n) = 0;
  // Then solve the tridiagonal matrix
  solver -> forward_solver(u,r);
  solver -> backward_solver(u,r);
  r = u;
  u(0) = 0;
  u(n) = 0;
  // Eventual print statements etc
  writetofile(outfile,u);
  }
  outfile.close();
  delete solver;
}

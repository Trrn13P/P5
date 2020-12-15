#include "tridiag_solver.hpp"
#include "diffusion.hpp"

#include <fstream>
using namespace std;

//writing u(t_j) to file if this function is called
void diffusion::writetofile(ofstream &outfile, vec u){
  outfile << u.t();
}

//g(x) boundary function, set to 1 for easy comparrison with analytic solution
float diffusion::func(float x){
  return 1;
}


void diffusion::forward_euler(){
  //opening outfile and creating first lines for reading in python
  ofstream outfile("../textfiles/forward_euler.txt");
  outfile << "dt=" << dt << " saved_tsteps=" << saved_tsteps  <<" alpha=" << alpha <<" dx=" << dx << endl;
  outfile << "u(t):, each new line is a new timestep * (saved_tsteps+1)" << endl;
  //setting diagonal elements
  a = alpha;
  b = (1-2*alpha);

  //Boundary conditions g(x)
  for (int i = 1; i < n+1; i++) {
    u(i) = func(dx*(i));
  }
  // Boundary conditions a(t) b(t) (zero here)
  unew(n) = u(n) = u(0) = unew(0) = 0;

  //writing initial conditions to file
  writetofile(outfile,u);

  //k is a counter for knowing when to print to file
  int k = 0;
  //going trough timesteps
  for (int t = 1; t <= tsteps; t++) {
    //going trough xsteps
    for (int i = 1; i < n; i++) {
    // Discretized diff eq
      unew(i) = alpha * u(i-1) + (1 - 2*alpha) * u(i) + alpha * u(i+1);
    }
    //setting boundary and updating u
    unew(0) = unew(n) = 0;
    u = unew;

    //printing to file
    if (k==saved_tsteps){
      writetofile(outfile,u);
      k=0;
    }
    else{
      k+=1;
    }
  }
  //closing outfile
  outfile.close();
}

void diffusion::backward_euler(){
  //opening outfile and creating first lines for reading in python
  ofstream outfile("../textfiles/backward_euler.txt");
  outfile << "dt=" << dt << " saved_tsteps=" << saved_tsteps  <<" alpha="<< alpha <<" dx=" << dx << endl;
  outfile << "u(t):, each new line is a new timestep * saved_tsteps" << endl;
  //setting diagonal elements
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

  int k = 0;
  //going trough timesteps
  for (int t = 1; t <= tsteps; t++) {
    //forward-backward solving
    solver->forward_solver(unew,u);
    solver->backward_solver(unew,u);

    //setting boundary and updating u
    unew(0) = unew(n) = 0;
    u = unew;

    //printing to file
    if (k==saved_tsteps){
      writetofile(outfile,u);
      k=0;
    }
    else{
      k+=1;
    }
  }
  //closing outfile
  outfile.close();
  delete solver;
}

void diffusion::crank_nicolson(){
  //opening outfile and creating first lines for reading in python
  ofstream outfile("../textfiles/crank_nicolson.txt");
  outfile << "dt=" << dt << " saved_tsteps=" << saved_tsteps  <<" alpha=" << alpha <<" dx=" << dx << endl;
  outfile << "u(t):, each new line is a new timestep * saved_tsteps" << endl;
  //setting diagonal elements
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

  int k = 0;
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
  //setting boundary and updating r
  r = u;
  u(0) = 0;
  u(n) = 0;

  //printing to file
  if (k==saved_tsteps){
    writetofile(outfile,u);
    k=0;
  }
  else{
    k+=1;
  }


  }
  //closing outfile
  outfile.close();
  delete solver;
}

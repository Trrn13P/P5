#include "tridiag_solver.hpp"
#include "diffusion.hpp"

#include <fstream>
#include <iostream>
using namespace std;

//writing u(t_j) to file if this function is called
void diffusion::writetofile(ofstream &outfile, vec u){
  outfile << u.t();
}

void diffusion::writetofile2d(ofstream &outfile, mat u_){
  outfile << u_ << endl;
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

void diffusion::forward_2d(){
  //opening outfile and creating first lines for reading in python
  ofstream outfile("../textfiles/2d_euler.txt");
  outfile  << "dt=" << dt << " saved_tsteps=" << saved_tsteps
  <<" alpha=" << alpha <<" dx=" << dx <<" n=" << n<< endl;
  outfile << "u(t):, each block of lines is at a time dt * (saved_tsteps+1) *current_timestep" << endl;

  //setting up matrixes
  mat u_, unew_;
  u_ = unew_ = mat(n+1,n+1);

  //setting boundary
  for(int i=0;i<n+1;i++){
    u_(0,i) = 0;
    u_(n,i) = 0;
    u_(i,0) = 0;
    u_(i,n) = 0;
  }
  //Setting initial conditions
  for(int i=1;i<n;i++){
    for(int j=1;j<n;j++){
      u_(i,j) = func(1);
    }
  }
  unew_ = u_;
  writetofile2d(outfile,u_);
  // Time iteration
  int k = 0;
  for (int t = 1; t <= tsteps; t++) {
    //Forward algorythm
    for(int i=1;i<n;i++){
      for(int j=1;j<n;j++){
        unew_(i,j) = u_(i,j) + alpha*(u_(i+1,j)+u_(i-1,j)+u_(i,j+1)+u_(i,j-1)-4*u_(i,j));
      }
    }

    //setting boundary and updating u
    for(int i=0;i<n+1;i++){
      unew_(0,i) = 0;
      unew_(n,i) = 0;
      unew_(i,0) = 0;
      unew_(i,n) = 0;
    }

    u_ = unew_;

    //printing to file
    if (k==saved_tsteps){
      writetofile2d(outfile,u);
      k=0;
    }
    else{
      k+=1;
    }
  }
  outfile.close();
}

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


void diffusion::backward_euler(){
  ofstream outfile("../textfiles/backward_euler.txt");
  a = c = -alpha;
  b = (1+2*alpha);

  //Boundary conditions g(x)
  for (int i = 1; i < n+1; i++) {
    unew(i) = u(i) = func(dx*(i));
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

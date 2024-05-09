#include "banded.hpp"
#include "pdma.hpp"

#include <fstream>

std::vector<std::vector<double> > iterate_step(penta_thomas_solver &solver, std::vector<std::vector<double>> &pos, std::vector<std::vector<double>> &vel, double dt, double m, double alpha, double f_bend, double L, int n_fixed){
    int n_balls = pos[0].size();
    std::vector<double> rhs(n_balls);
    std::vector<double> sol(n_balls);
    std::vector<std::vector<double> > new_pos;
    new_pos.resize(3);

  for(int j = 0; j < 3; j++){

    for(int i=0; i<n_fixed; i++){
        rhs[i] = pos[j][i];
    }
    for(int i=n_fixed; i<n_balls; i++){
        rhs[i] = pos[j][i] * (m/dt/dt + alpha/dt) + vel[j][i]*m/dt;
    }
    if(j == 2){ rhs[n_balls-1] += f_bend; }
    sol = solver.solve(rhs);
    new_pos[j] = sol;
  }
  return new_pos;
};

int main(){


double time = 1.0; // Seconds
double dt = 1e-2;
int n_iter = time/dt;

double real_len = 10.0; // microns
int n_balls = 20;
double m = 3e-4 * real_len / (double)n_balls; // picograms

double L = real_len/(float)n_balls;

double  k =  1.0e16; // Spring constant
double  alpha = 5.0e8;// Viscosity

int n_fixed = 2;

double f_bend = 1.0e9; // Converts to picoNewtons

penta_repeating mat(n_balls, 2, 2);
//banded_general_penta mat(n_balls);

mat.fill_from_elemBanded(elementary_matrices::SpringBanded);
mat.multiply_all(k);
double correc = m/dt/dt + alpha/dt;
mat.add_identity_factor(-correc);
mat.multiply_all(-1.0);

for(int i=0; i<n_fixed; i++){
    mat.set_identity_row(i);
}

std::cout<<"Matrix is :"<<std::endl;
mat.print();
std::cout<<"Or :"<<std::endl;
mat.expanded_print();

penta_thomas_solver solver;
solver.update_array(mat);
std::cout<<"Verifying LU decomposition"<<std::endl;
bool err = solver.verify_stored_array(mat);
if(err) std::cout<<"Problem in decmposition"<<std::endl;

std::vector<std::vector<double> > pos, old_pos, vel;

pos.resize(3);
vel.resize(3);
for(auto & p : pos) p.resize(n_balls);
for(auto & p : vel) p.resize(n_balls);

for(int i=1; i<n_balls; i++){
    pos[0][i] = pos[0][i-1] + L;
}

for(int t = 0; t<n_iter;t++){
    old_pos = pos;
    pos = iterate_step(solver, pos, vel, dt, m, alpha, f_bend, L, n_fixed);
    for(int j = 0; j<3; j++){
      for(int i=0; i<n_balls; i++){
        vel[j][i] = (pos[j][i] - old_pos[j][i])/dt;
      }
    }
}

std::ofstream outfile;
outfile.open("output.csv");
for(int i=0; i<n_balls; i++) outfile<< pos[0][i]<<","<<pos[1][i]<<","<<pos[2][i]<<"\n";
std::cout<<std::endl;
outfile.close();

}

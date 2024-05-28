#include "banded.hpp"
#include "pdma.hpp"

int main(){

const int len = 8;

//banded_general_penta mat(len);
penta_repeating mat(len, 2, 2);

mat.fill_from_elemBanded(elementary_matrices::SpringBanded);
mat.add_identity_factor(0.1);
std::cout<<"Matrix is :"<<std::endl;
mat.print();
std::cout<<"Matrix in full is :"<<std::endl;
mat.pretty_print();

penta_thomas_solver solver;
solver.update_array(mat);
std::cout<<"Verifying LU decomposition"<<std::endl;
bool err = solver.verify_stored_array(mat);
if(err) std::cout<<"Problem in decomposition"<<std::endl;

std::vector<double> rhs, sol;
rhs.resize(len);
sol.resize(len);
for(int i=1; i<len; i++){
    sol[i] = sol[i-1] + 0.1;
}
long n_iter = 10;
for(int t = 0; t<n_iter;t++){
    // Do a bunch of solves, checking the result is within the expected error
    mat.clear();
    mat.fill_from_elemBanded(elementary_matrices::SpringBanded);
    mat.add_identity_factor(0.1);
    solver.update_array(mat);
    for(int i=0; i<len; i++){
      rhs[i] = sol[i] * 0.1;
    }
    sol = solver.solve(rhs);
    auto b = matvecmult(mat, sol);
    double max_err = 0.0;
    for(int i=0; i<len; i++){
      if(std::abs(b[i] - rhs[i]) > max_err){
        max_err = std::abs(b[i] - rhs[i]);
      }
    }
    // Check ABSOLUTE error
    std::cout<<"Iteration "<<t<<", error "<<max_err<<std::endl;
}

  // Try out the nudge resizing
  std::cout<<"Growing matrix \n";
  int new_len = len + 2;
  mat.change_len(new_len);
  mat.pretty_print();
  solver.nudge_array_len(mat);
  err = solver.verify_stored_array(mat);
  if(err) std::cout<<"Problem in decomposition "<<std::endl;

  rhs.resize(new_len);
  for(int i=1; i<len; i++){
    rhs[i] =  rhs[i-1] + 0.01; // Same original rhs
  }
  sol = solver.solve(rhs);
  auto b = matvecmult(mat, sol);
  double max_err = 0.0;
  for(int i=0; i<len; i++){
    if(std::abs(b[i] - rhs[i]) > max_err){
      max_err = std::abs(b[i] - rhs[i]);
    }
  }
  // Check ABSOLUTE error
  std::cout<<"Resized array, error "<<max_err<<std::endl;

}

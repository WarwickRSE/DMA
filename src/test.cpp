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
if(err) std::cout<<"Problem in decmposition"<<std::endl;

std::vector<double> rhs, sol;
rhs.resize(len);
for(int i=1; i<len; i++){
    rhs[i] = rhs[i-1] + 0.1;
}
long n_iter = 1;
for(int t = 0; t<n_iter;t++){
    mat.clear();
    mat.fill_from_elemBanded(elementary_matrices::SpringBanded);
    mat.add_identity_factor(0.1);
    solver.update_array(mat);
    rhs[0] += 0.1 * float(t)/float(n_iter);
    sol = solver.solve(rhs);
    std::cout<<"Print to trick the optimiser: "<<sol[0]<<std::endl;
}

for(int i=0; i<len; i++) std::cout<<rhs[i]<<",";
std::cout<<std::endl;

for(int i=0; i<len; i++) std::cout<<sol[i]<<",";
std::cout<<std::endl;


  // Try out the nudge resizing
  std::cout<<"Growing matrix \n";
  int new_len = len - 3;
  mat.change_len(new_len);
  mat.pretty_print();
  solver.nudge_array_len(mat);
  err = solver.verify_stored_array(mat);
  if(err) std::cout<<"Problem in decomposition "<<std::endl;

  rhs.resize(new_len);
  for(int i=1; i<new_len; i++){
    rhs[i] = rhs[i-1] + 0.1;
  }
  sol = solver.solve(rhs);
  std::cout<<sol[0]<<", "<<sol[new_len-1]<<std::endl;

}

#include "banded.hpp"
#include "pdma.hpp"

int main(){

const int len = 200;

banded_general_penta mat(len);

mat.fill_from_elemBanded(elementary_matrices::SpringBanded);
mat.add_identity_factor(0.1);
mat.print();

reference_penta_thomas solver;
solver.update_array(mat);

solver.verify_stored_array(mat);

std::vector<double> rhs, sol;
rhs.resize(len);
for(int i=1; i<len; i++){
    rhs[i] = rhs[i-1] + 0.1;
}
long n_iter = 100000;
for(int t = 0; t<n_iter;t++){
    mat.clear();
    mat.fill_from_elemBanded(elementary_matrices::SpringBanded);
    mat.add_identity_factor(0.1);
    solver.update_array(mat);
    rhs[0] += 0.1 * float(t)/float(n_iter);
    sol = solver.solve(rhs);
    std::cout<<sol[0]<<std::endl;
}

for(int i=0; i<len; i++) std::cout<<rhs[i]<<",";
std::cout<<std::endl;

for(int i=0; i<len; i++) std::cout<<sol[i]<<",";
std::cout<<std::endl;

}

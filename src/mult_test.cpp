
#include "banded.hpp"
#include <iostream>
#include <cassert>

int main(){

  const int len = 6;

  // Bandwidth 1
  banded_general<1> mat1_1(len), mat1_2(len);
  //Fill in
  for(int i = 0; i < len; i++){
    mat1_1.set(0, i, i);
  }
  mat1_2 = mat1_1;
  mat1_2.multiply_all(-1.0);

  auto mat1_3 = matmul(mat1_1, mat1_2);
  std::cout<<"Multiplying two one-diagonal matrices\n";
  mat1_3.print();
  std::cout<<"______________________\n";
  std::cout<<"Result has bandwidth "<<diag_bandw(mat1_3)<<std::endl;
  std::cout<<"______________________\n";
 
  // Bandwidth 3
  banded_general<3> mat3_1(len), mat3_2(len);
  // Fill in
  for(int i = 0; i < len; i++){
      mat3_1.set(1, i, i+1);
      if(i < len-1){
          mat3_1.set(0, i, i+1);
          mat3_1.set(2, i, i+1);
      }
  }
  mat3_2 = mat3_1;
  mat3_2.multiply_all(-1.0);

  auto mat3_3 = matmul(mat3_1, mat3_2);
  std::cout<<"Multiplying two tridiagonal matrices\n";
  mat3_3.print();
  std::cout<<"______________________\n";
  std::cout<<"Result has bandwidth "<<diag_bandw(mat3_3)<<std::endl;
  std::cout<<"______________________\n";

  banded_general<5> mat5_1(len), mat5_2(len);

    // Fill in
    for(int i = 0; i < len; i++){
        mat5_1.set(2, i, i+1);
        if(i < len-1){
            mat5_1.set(1, i, i+1);
            mat5_1.set(3, i, i+1);
        }
        if(i < len-2){
            mat5_1.set(0, i, i+1);
            mat5_1.set(4, i, i+1);
        }
    }
    mat5_2 = mat5_1;
    mat5_2.multiply_all(-1.0);

    std::cout<<"Checking equality of matrices\n";
    assert(mat5_1 != mat5_2);
    assert(mat5_1 == mat5_1);
    auto mat5_1_copy = mat5_1.to_general();
    assert(mat5_1 == mat5_1_copy);
    std::cout<<"..... done\n";

    auto mat5_3 = matmul(mat5_1, mat5_2);
    std::cout<<"Multiplying two penta-diagonal matrices\n";
    mat5_3.print();
    std::cout<<"______________________\n";
    std::cout<<"Result has bandwidth "<<diag_bandw(mat5_3)<<std::endl;
    std::cout<<"______________________\n";

    //Checking against the expanded multiplication
    auto mat_e1 = mat5_1.to_general();
    auto mat_e2 = mat5_2.to_general();

    auto mat_e3 = matmul(mat_e1, mat_e2);
    std::cout<<"Checking banded and non banded multiplication\n";
    assert(mat_e3 == mat5_3);
    std::cout<<"..... done\n";

    //Looking at repeating matrices
    penta_repeating mat_repeating(len, 2, 2);
    mat_repeating.fill_from_elemBanded(elementary_matrices::SpringBanded);

    auto mat_repeating_expanded = mat_repeating.to_full();
    std::cout<<"Creating a banded matrix from a repeating one\n";
    mat_repeating.print();
    std::cout<<"______________________\n";
    mat_repeating_expanded.print();
    assert(mat_repeating == mat_repeating_expanded);
    std::cout<<"..... done\n";

    // Conversion to real matrix
    auto mat_repeating_general = mat_repeating.to_general();
    std::cout<<"Now creating full square matrix\n";
    mat_repeating_general.print();
    std::cout<<"______________________\n";
    mat_repeating.expanded_print();
    std::cout<<"______________________\n";
    assert(mat_repeating_general == mat_repeating_expanded);
    std::cout<<"..... done\n";
 
    // Matmul on repeating
    std::cout<<"Verifying multiplication with repeating matrices\n";
    auto mat_repeating_2 = matmul(mat_repeating, mat_repeating);
    mat_repeating_2.print();
    std::cout<<"______________________\n";
    auto mat_expanded_2 = matmul(mat_repeating_general, mat_repeating_general);
    mat_expanded_2.print();
    assert(mat_repeating_2 == mat_expanded_2);
    std::cout<<"..... done\n";

}
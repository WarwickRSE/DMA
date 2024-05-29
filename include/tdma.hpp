#ifndef TDMA_HPP
#define TDMA_HPP
#include <algorithm>
#include <iostream>
#include <cassert>

#include "banded.hpp"

const double zero_thresh = 1e-14;

class thomas_solver{
/**
 * @brief Solves A x = b where A is a tri-diagonal matrix
 * 
 * Optimised for the case where matrix is small enough to hold 3xN elements in memory, and where it either does not change size, or changes slowly, as re-calculation is minimised; and where a large number of solves are done for a given matrix. C.f. thomas_onthefly for case where calculation on the fly is superior.
 * 
 */

  std::vector<std::vector<double> > lu_values;
  size_t len;
  const int w=0, b_st=1, c_st=2;
  const int a=0, b=1, c=2;
  public:

  // Update stored values for a new array
  // Will completely wipe out any previous values
  template < typename T >
  void update_array(T const & new_arr){

    len = new_arr.len;

    lu_values.clear();
    lu_values.resize(3); // cprime and dprime (Wiki)
    for(int i = 0; i<3; i++) lu_values[i].resize(len);

    lu_values[b_st][0] = new_arr.get(b, 0);
    for(size_t i = 1; i < len; i++){
        lu_values[w][i] = new_arr.get(a, i-1)/lu_values[b_st][i-1];
        lu_values[b_st][i] = new_arr.get(b, i) - lu_values[w][i] * new_arr.get(c, i-1);
        lu_values[c_st][i-1] = new_arr.get(c, i-1);
    }
  }

  void nudge_array_len(tri_repeating const &new_arr){
    // MUST pass the same array as the one originally set, EXCEPT the length is changed
    // Short-cuts certain logic for faster update on resize, esp growing/shrinking by amount << len
    // Only recalculates the last ftr + delta_len rows (for growth) or ftr + 1 for shrinkage
    // MAY cause memory re-alocation internally, especially on growing
    // Only really useful for a repeating type matrix, because otherwise where do the new values come from?

    size_t old_len = len;
    len = new_arr.len;
    long delta_len = len - old_len;
    bool grow = delta_len > 0;
    size_t start_of_recalc = len - new_arr.ftr;
    if(grow){
        start_of_recalc -= delta_len;
    }else{
        start_of_recalc -= 1;
    }

    for(size_t i = 0; i < 3; i++){
        lu_values[i].resize(len);
    }

    for(size_t i = start_of_recalc; i < len; i++){
        lu_values[w][i] = new_arr.get(a, i-1)/lu_values[b_st][i-1];
        lu_values[b_st][i] = new_arr.get(b, i) - lu_values[w][i] * new_arr.get(c, i-1);
        lu_values[c_st][i-1] = new_arr.get(c, i-1);
    }

  }

template < typename T >
bool verify_stored_array(T& expected){
    banded_general_tri L(len), U(len);

    if(len != expected.len){
        std::cout<<"Canot verify, length incorrect\n";
        return true;
    }

    // Construct them
    for(size_t i = 0; i< len; i++){
        L.set(1, i, 1.0);
        if(i < len-1){
          L.set(0, i, lu_values[w][i+1]);
        }
    }

    for(size_t i = 0; i < len; i++){
        U.set(1, i, lu_values[b_st][i]);
        if(i < len-1){
          U.set(2, i, lu_values[c_st][i]);
        }
    }

    // Multiply them
    // Obtain a non-banded matrix in general. Best to check for absence of non-zeros off the band to be sure
    auto res = matmul(L, U);

    bool err = false;
    double max_err = 0.0;
    double non_zero_err = 0.0;
    for(size_t i = 0; i<len; i++){
      for(size_t j = 0; j<len; j++){
        if(std::abs((long)(i - j)) < 2){
          if(expected.get_real(i,j) != 0.0){// Exact equality correct, just protecting from divide-by-zero
            double diff = std::abs((res.get(i, j) - expected.get_real(i, j))/expected.get_real(i, j));
            if(diff > zero_thresh && diff > max_err){
              max_err = diff;
              err = true;
            }
          }
        }else{
          double val = std::abs(res.get(i, j));
          if(val > zero_thresh && val > non_zero_err){
            non_zero_err = val;
            err = true;
          }
        }
      }
    }
    if(err) res.print();
    std::cout<<"Max diff from expected: "<<max_err<<std::endl;
    std::cout<<"Max val of element which should be zero: "<<non_zero_err<<std::endl;
    assert(!err);
    return err;
}

std::vector<double> solve(const std::vector<double> & rhs)const{

#ifdef DEBUG
  assert(len == rhs.size());
#endif

  std::vector<double> rho, psi;
  rho.resize(len);
  psi.resize(len);

  // Forward sweep
  rho[0] = rhs[0];
  for(size_t i=1; i<len; i++){
    rho[i] = rhs[i] - lu_values[w][i] * rho[i-1];
  }

  // Reverse
  psi[len-1] = rho[len-1]/lu_values[b_st][len-1];
  for(long i = len-2; i > -1; i--){
    psi[i] = (rho[i] - lu_values[c_st][i] * psi[i+1])/lu_values[b_st][i];
  }

  return psi;
}

};

class thomas_onthefly{
/**
 * @brief Solves A x = b where A is a penta-diagonal matrix
 * 
 * Solves the system on the fly, without storing intermediate values. This is useful when the matrix is large and the solves are infrequent, or when the matrix changes size frequently, so that storing intermediate values is not efficient 
 */

  public:

  template < typename T >
  std::vector<double> solve(const T & A, const std::vector<double> & rhs)const{

    const int a=0, b=1, c=2;

    auto len = A.len;
#ifdef DEBUG
  assert(len == rhs.size());
#endif

    std::vector<double> rho, psi, b_pr;
    b_pr.resize(len);
    rho.resize(len);
    psi.resize(len);

    // Forward sweep and store new b for later
    double w;
    rho[0] = rhs[0];
    b_pr[0] = A.get(b, 0);
    for(size_t i=1; i<len; i++){
      w = A.get(a, i-1)/b_pr[i-1];
      b_pr[i] = A.get(b, i) - w * A.get(c, i-1);
      rho[i] = rhs[i] -  w * rho[i-1];
    }

    // Reverse
    psi[len-1] = rho[len-1]/b_pr[len-1];
    for(long i = len-2; i > -1; i--){
      psi[i] = (rho[i] - A.get(c, i) * psi[i+1])/b_pr[i];
    }

    return psi;
  }
};

#endif

#ifndef PDMA_HPP
#define PDMA_HPP
#include <algorithm>
#include <iostream>
#include <cassert>

#include "banded.hpp"

const double zero_thresh = 1e-14;

class penta_thomas_solver{
/**
 * @brief Solves A x = b where A is a penta-diagonal matrix
 * 
 * Optimised for the case where matrix is small enough to hold 5xN elements in memory, and where it either does not change size, or changes slowly, as re-calculation is minimised.
 * 
 */

  std::vector<std::vector<double> > lu_values;
  size_t len;
  const int x=0, y=1, z=2, c_st=3, e_st=4;
  const int a=2, b=3, c=4, d=1, e=0;
  public:

  // Update stored values for a new array
  // Will completely wipe out any previous values
  template < typename T >
  void update_array(T const & new_arr){

    len = new_arr.len;

    lu_values.clear();
    lu_values.resize(5); // x, y, z
    for(int i = 0; i<5; i++) lu_values[i].resize(len);

    // Initial values
    lu_values[x][0] = new_arr.get(a, 0);
    lu_values[y][0] = new_arr.get(b, 0);
    lu_values[z][1] = new_arr.get(d, 0) / lu_values[x][0];
    lu_values[y][1] = new_arr.get(b, 1) - lu_values[z][1] * new_arr.get(c,0);
    lu_values[x][1] = new_arr.get(a, 1) - lu_values[y][0] * lu_values[z][1];

    for(size_t i = 2; i < len; i++){
        lu_values[z][i] = ( new_arr.get(d,i-1) - (new_arr.get(e, i-2) * lu_values[y][i-2]/lu_values[x][i-2]))/lu_values[x][i-1];
        if(i < len-1){
          lu_values[y][i] = new_arr.get(b, i) -  lu_values[z][i]* new_arr.get(c, i-1);
        }
        lu_values[x][i] = new_arr.get(a, i) - (lu_values[y][i-1] * lu_values[z][i]) - (new_arr.get(e, i-2) * new_arr.get(c, i-2)/lu_values[x][i-2]);
    }

    // Renumber for simpler maths
    for(size_t i = 0; i < len-2; i++){
        lu_values[c_st][i] = new_arr.get(c, i);
        lu_values[e_st][i+2] = new_arr.get(e, i)/lu_values[x][i];
    }

#ifdef DEBUG
    // ALL x's must be non-zero
    assert(std::all_of(lu_values[x].begin(), lu_values[x].end(), [](double d){return std::abs(d)>zero_thresh;}));
#endif

  }

  void nudge_array_len(penta_repeating const &new_arr){
    // MUST pass the same array as the one originally set, EXCEPT the length is changed
    // Short-cuts certain logic for faster update on resize, esp growing/shrinking by amount << len
    // Only recalculates the last ftr + delta_len rows (for growth) or ftr + 1 for shrinkage
    // MAY cause memory re-alocation internally, especially on growing

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

    for(size_t i = 0; i < 5; i++){
        lu_values[i].resize(len);
    }

    for(size_t i = start_of_recalc; i < len; i++){
        lu_values[z][i] = ( new_arr.get(d,i-1) - (new_arr.get(e, i-2) * lu_values[y][i-2]/lu_values[x][i-2]))/lu_values[x][i-1];
        if(i < len-1){
          lu_values[y][i] = new_arr.get(b, i) -  lu_values[z][i]* new_arr.get(c, i-1);
        }
        lu_values[x][i] = new_arr.get(a, i) - (lu_values[y][i-1] * lu_values[z][i]) - (new_arr.get(e, i-2) * new_arr.get(c, i-2)/lu_values[x][i-2]);
    }

    // Renumber for simpler maths
    for(size_t i = start_of_recalc; i < len-2; i++){
        lu_values[c_st][i] = new_arr.get(c, i);
        lu_values[e_st][i+2] = new_arr.get(e, i)/lu_values[x][i];
    }

  }

template < typename T >
bool verify_stored_array(T expected){
    banded_general_penta L(len), U(len);

    if(len != expected.len){
        std::cout<<"Canot verify, length incorrect\n";
        return true;
    }

    // Construct them
    for(size_t i = 0; i< len; i++){
        L.set(2, i, 1.0);
        if(i < len-1){
          L.set(1, i, lu_values[z][i+1]);
        }
        if(i < len-2){
          L.set(0, i, lu_values[e_st][i+2]);
        }
    }

    for(size_t i = 0; i < len; i++){
        U.set(2, i, lu_values[x][i]);
        if(i < len-1){
          U.set(3, i, lu_values[y][i]);
        }
        if(i < len-2){
          U.set(4, i, lu_values[c_st][i]);
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
        if(std::abs((long)(i - j)) < 3){
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

std::vector<double> solve(const std::vector<double> rhs){

#ifdef DEBUG
  assert(len == rhs.size());
#endif

  std::vector<double> rho, psi;
  rho.resize(len);
  psi.resize(len);

  // Forward sweep
  rho[0] = rhs[0];
  rho[1] = rhs[1] - lu_values[z][1] * rho[0];
  for(size_t i =2; i<len; i++){
    rho[i] = rhs[i] - lu_values[z][i] * rho[i-1] - lu_values[e_st][i]*rho[i-2];
  }

  // Reverse
  psi[len-1] = rho[len-1]/lu_values[x][len-1];
  psi[len-2] = (rho[len-2] - lu_values[y][len-2]*psi[len-1])/lu_values[x][len-2];
  for(long i = len-3; i > -1; i--){
    psi[i] = (rho[i] - lu_values[y][i] * psi[i+1] - lu_values[c_st][i] * psi[i+2])/lu_values[x][i];
  }

  return psi;
}

};
#endif

#ifndef PDMA_HPP
#define PDMA_HPP
#include <algorithm>
#include <iostream>
#include <cassert>

#include "banded.hpp"

const double zero_thresh = 1e-14;

class penta_thomas_solver{

  std::vector<std::vector<double> > lu_values;
  int len;
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

    for(int i = 2; i < len; i++){
        lu_values[z][i] = ( new_arr.get(d,i-1) - (new_arr.get(e, i-2) * lu_values[y][i-2]/lu_values[x][i-2]))/lu_values[x][i-1];
        if(i < len-1){
          lu_values[y][i] = new_arr.get(b, i) -  lu_values[z][i]* new_arr.get(c, i-1);
        }
        lu_values[x][i] = new_arr.get(a, i) - (lu_values[y][i-1] * lu_values[z][i]) - (new_arr.get(e, i-2) * new_arr.get(c, i-2)/lu_values[x][i-2]);
    }

    // Renumber for simpler maths
    for(int i = 0; i < len-2; i++){
        lu_values[c_st][i] = new_arr.get(c, i);
        lu_values[e_st][i+2] = new_arr.get(e, i);
    }

#ifdef DEBUG
    // ALL x's must be non-zero
    assert(std::all_of(lu_values[x].begin(), lu_values[x].end(), [](int i){return std::abs(i)>zero_thresh}));
#endif

  }

  void nudge_array_len(penta_repeating const &new_arr){
    // MUST pass the same array as the one originally set, EXCEPT the length is changed
    // Short-cuts certain logic for faster update on resize, esp growing/shrinking by amount << len
    // Only recalculates the last ftr + delta_len rows (for growth) or ftr + 1 for shrinkage
    // MAY cause memory re-alocation internally, especially on growing

    int old_len = len;
    len = new_arr.len;
    int delta_len = len - old_len;
    bool grow = delta_len > 0;
    int start_of_recalc = len - new_arr.ftr;
    if(grow){
        start_of_recalc -= delta_len;
    }else{
        start_of_recalc -= 1;
    }

    for(int i = 0; i < 5; i++){
        lu_values[i].resize(len);
    }

    for(int i = start_of_recalc; i < len; i++){
        lu_values[z][i] = ( new_arr.get(d,i-1) - (new_arr.get(e, i-2) * lu_values[y][i-2]/lu_values[x][i-2]))/lu_values[x][i-1];
        if(i < len-1){
          lu_values[y][i] = new_arr.get(b, i) -  lu_values[z][i]* new_arr.get(c, i-1);
        }
        lu_values[x][i] = new_arr.get(a, i) - (lu_values[y][i-1] * lu_values[z][i]) - (new_arr.get(e, i-2) * new_arr.get(c, i-2)/lu_values[x][i-2]);
    }

    // Renumber for simpler maths
    for(int i = start_of_recalc; i < len-2; i++){
        lu_values[c_st][i] = new_arr.get(c, i);
        lu_values[e_st][i+2] = new_arr.get(e, i);
    }

  }

template < typename T >
bool verify_stored_array(T expected){
    banded_general_penta L(len), U(len), res(len);

    if(len != expected.len){
        std::cout<<"Canot verify, length incorrect\n";
        return true;
    }

    // Construct them
    for(int i = 0; i< len; i++){
        L.set(2, i, 1.0);
        if(i < len-1){
          L.set(1, i, lu_values[z][i+1]);
        }
        if(i < len-2){
          L.set(0, i, lu_values[e_st][i+2]/ lu_values[x][i]);
        }
    }

    for(int i = 0; i< len; i++){
        U.set(2, i, lu_values[x][i]);
        if(i < len-1){
          U.set(3, i, lu_values[y][i]);
        }
        if(i < len-2){
          U.set(4, i, lu_values[c_st][i]);
        }
    }

    const int a=2, b=3, c=4, d=1, e=0;

    // Multiply them
    res.set(a, 0, L.get(a, 0)*U.get(a, 0));
    res.set(a, 1, L.get(d, 0)*U.get(b, 0) + L.get(a, 1)*U.get(a, 1));

    res.set(d, 0, L.get(d, 0)*U.get(a, 0));
    res.set(b, 0, L.get(a, 0)*U.get(b, 0));

    for(int i = 0; i < len; i++){
        if(i > 1){
          res.set(a, i, L.get(e, i-2)*U.get(c, i-2) + L.get(d, i-1)*U.get(b, i-1) + L.get(a, i)*U.get(a, i));
        }
        if(i > 0 && i < len-1){
          res.set(b, i, L.get(d, i-1)*U.get(c, i-1) + L.get(a, i)*U.get(b, i));
          res.set(d, i, L.get(e, i-1)*U.get(b, i-1) + L.get(d, i)*U.get(a, i));
        }
        if(i < len-2){
          res.set(c, i, L.get(a, i)*U.get(c, i));
          res.set(e, i, L.get(e, i)*U.get(a, i));
        }
    }

    bool err = false;
    double max_err = 0.0;
    for(int i = 0; i < 5; i++){
        int l = res.len - std::abs(i-2);
        for(int j = 0; j< l; j++){
            double diff = std::abs((res.get(i, j) - expected.get(i, j))/expected.get(i, j));
            if(diff > max_err) max_err = diff;
            if(diff > zero_thresh)err = true;
        }
    }
    if(err) res.print();
    std::cout<<"Max error: "<<max_err<<std::endl;
    assert(!err);
    return err;
}

std::vector<double> solve(const std::vector<double> rhs){

  double thresh = 1e-14;

  int sz = len; // Avoid rename for now

#ifdef DEBUG
  assert(sz == rhs.size());
#endif

  std::vector<double> rho, psi;
  rho.resize(len);
  psi.resize(len);

  // Forward sweep
  rho[0] = rhs[0];
  rho[1] = rhs[1] - lu_values[z][1] * rho[0];
  for(int i =2; i<sz; i++){
    rho[i] = rhs[i] - lu_values[z][i] * rho[i-1] - lu_values[e_st][i]*rho[i-2]/lu_values[x][i-2];
  }

  // Reverse
  psi[sz-1] = rho[sz-1]/lu_values[x][sz-1];
  psi[sz-2] = (rho[sz-2] - lu_values[y][sz-2]*psi[sz-1])/lu_values[x][sz-2];
  for(int i = sz-3; i > -1; i--){
    psi[i] = (rho[i] - lu_values[y][i] * psi[i+1] - lu_values[c_st][i] * psi[i+2])/lu_values[x][i];
  }

  return psi;
}

};
#endif

#include <algorithm>
#include <iostream>

#include "banded.hpp"

const double zero_thresh = 1e-14;

class reference_penta_thomas{

  std::vector<std::vector<double> > lu_values;
  int len;
  const int x=0, y=1, z=2, c_st=3, e_st=4;

  public:

  // Update stored values for a new array (including size changes)
  // Will completely wipe out any previous values
  void update_array(banded_general_penta new_arr){

    len = new_arr.len;
    const int a=2, b=3, c=4, d=1, e=0;

    lu_values.clear();
    lu_values.resize(5); // x, y, z
    for(int i = 0; i<5; i++) lu_values[i].resize(len);

    // Initial values
    lu_values[x][0] = new_arr.values[a][0];
    lu_values[y][0] = new_arr.values[b][0];
    lu_values[z][1] = new_arr.values[d][0] / lu_values[x][0];
    lu_values[y][1] = new_arr.values[b][1] - lu_values[z][1] * new_arr.values[c][0];
    lu_values[x][1] = new_arr.values[a][1] - lu_values[y][0] * lu_values[z][1];

    for(int i = 2; i < len; i++){
        lu_values[z][i] = ( new_arr.values[d][i-1] - (new_arr.values[e][i-2] * lu_values[y][i-2]/lu_values[x][i-2]))/lu_values[x][i-1];
        lu_values[y][i] = new_arr.values[b][i] -  lu_values[z][i]* new_arr.values[c][i-1];
        lu_values[x][i] = new_arr.values[a][i] - (lu_values[y][i-1] * lu_values[z][i]) - (new_arr.values[e][i-2] * new_arr.values[c][i-2]/lu_values[x][i-2]);
    }

    // Renumber for simpler maths
    for(int i = 0; i < len-2; i++){
        lu_values[c_st][i] = new_arr.values[c][i];
        lu_values[e_st][i+2] = new_arr.values[e][i];
    }

#ifdef DEBUG
    assert(std::all_of(lu_values[x].begin(), lu_values[x].end(), [](int i){return std::abs(i)<zero_thresh}));
#endif

  }

bool verify_stored_array(banded_general_penta expected){
    banded_general_penta L(len), U(len), res(len);

    // Construct them
    for(int i = 0; i< len; i++){
        L.values[2][i] = 1.0;
        if(i < len-1){
          L.values[1][i] = lu_values[z][i+1];
        }
        if(i < len-2){
          L.values[0][i] = lu_values[e_st][i+2]/ lu_values[x][i];
        }
    }

    for(int i = 0; i< len; i++){
        U.values[2][i] = lu_values[x][i];
        U.values[3][i] = lu_values[y][i];
        if(i < len-2){
          U.values[4][i] = lu_values[c_st][i];
        }
    }

    const int a=2, b=3, c=4, d=1, e=0;

    // Multiply them
    res.values[a][0] = L.values[a][0]*U.values[a][0];
    res.values[a][1] = L.values[d][0]*U.values[b][0] + L.values[a][1]*U.values[a][1];

    res.values[d][0] = L.values[d][0]*U.values[a][0];
    res.values[b][0] = L.values[a][0]*U.values[b][0];

    for(int i = 0; i < len; i++){
        if(i > 1){
          res.values[a][i] = L.values[e][i-2]*U.values[c][i-2] + L.values[d][i-1]*U.values[b][i-1] + L.values[a][i]*U.values[a][i];
        }
        if(i > 0 && i < len-1){
          res.values[b][i] = L.values[d][i-1]*U.values[c][i-1] + L.values[a][i]*U.values[b][i];
          res.values[d][i] = L.values[e][i-1]*U.values[b][i-1] + L.values[d][i]*U.values[a][i];
        }
        if(i < len-2){
          res.values[c][i] = L.values[a][i]*U.values[c][i];
          res.values[e][i] = L.values[e][i]*U.values[a][i];
        }
    }

    bool err = false;
    for(int i = 0; i < 5; i++){
        for(int j = 0; j< res.values[i].size(); j++){
            if(std::abs(res.values[i][j] - expected.values[i][j]) > zero_thresh) err = true;
        }
    }
    assert(!err);
    if(err) res.print();

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
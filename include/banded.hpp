// Compact storage of low-bandwidth diagonal matrices
// Customised for dynamics matrices: some number of boundary rows
// and a central block

#ifndef __bandedHPP__
#define __bandedHPP__

#include <vector>
#include <cmath>
#include <iostream>

namespace elementary_matrices{

    // Values. Bandwidth is (sz*2 -1); Simplify by not nesting...
    const std::vector<double> Spring{-1, 2, -1,   2, -4, 2,   -1, 2, -1};
    const std::vector<double> SpringBanded{-1, 2, -1, 2, -1,  2, -4, 2,  -1};
};

template <int bandw> 
class banded_general{

  std::vector< std::vector<double> > values;

  public:
  const int len;

  banded_general(const int len_in):len(len_in){
    // Checks
#ifdef DEBUG
    assert(bandw % 2 == 1); // Bandwidth has to be odd (main diag plus 2*off-diag)
    assert(len > bandw); // Doesn't make any sense otherwise...
#endif 

    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
    values.resize(bandw);
    for(int i = 0; i < bandw; ++i){
        int l = len - std::abs(i-centr);
        values[i].resize(l);
    }
  }
  void clear(){
    for(int i = 0; i < bandw; i++){
        for(double & value : values[i]){
            value = 0.0;
        }
    }
  }

  void fill_from_otherstuff(){ //TODO
  };


  void fill_from_elemBanded(const std::vector<double> elem){
    // ADDS the elementary banded array contributions
    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
    const int elem_sz = std::sqrt(elem.size()); // Size of elementary
    const int elem_bw = elem_sz * 2 - 1;

    for(int j = 0; j < len - elem_sz+1; j++){
      int offset = 0;
      for(int e = 0; e < elem_sz; e++){
        for(int i = e; i < bandw-e; i++){
          values[i][j+e] += elem[offset + i];
        }
        offset += bandw - e*2 - 1; 
      }
    }
  }
  void add_identity_factor(double factor){
    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
    for(int i = 0; i < len; i++){
        values[centr][i] += factor;
    }
  }

  void print(){
    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
    for(int i =0; i < bandw; i++){
        int l = len - std::abs(i-centr);
        for(int j = 0; j < l; j++ ){
            std::cout<< values[i][j]<<"\t";
        }
        std::cout<<'\n';
    }
  }
  void pretty_print(){print();}

  double get(int i, int j) const{
    return values[i][j];
  }
  void set(int i, int j, double val){
    values[i][j] = val;
  }

};

template <int bandw>
class banded_repeating{
// A banded matrix where the middle rows are identical - saves storing them all
// The number of start and end rows is fixed, but the length is allowed to vary by
// changing the number of repeats

  private:
  std::vector< std::vector<double> > values;

  public:
  int len;
  int reps;
  const int hdr; // Number of initial rows
  const int ftr; // Number of final rows

  banded_repeating(const int len_in, const int hdr_in, const int ftr_in):len(len_in), hdr(hdr_in), ftr(ftr_in){
    // Checks
#ifdef DEBUG
    assert(bandw % 2 == 1); // Bandwidth has to be odd (main diag plus 2*off-diag)
    assert(len > bandw); // Doesn't make any sense otherwise...
    assert(len > hdr+ftr+1); // Must have at least one repeating row
#endif 

    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
    values.resize(bandw);
    reps = len - hdr - ftr;
    int stored_len = hdr + ftr + 1; // Length required to store elements
    for(int i = 0; i < bandw; ++i){
        int l = stored_len - std::abs(i-centr);
        values[i].resize(l);
    }
  }
  void clear(){
    for(int i = 0; i < bandw; i++){
        for(double & value : values[i]){
            value = 0.0;
        }
    }
  }

  void fill_from_otherstuff(){ //TODO
  };


  void fill_from_elemBanded(const std::vector<double> elem){
    // ADDS the elementary banded array contributions
    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing

    const int elem_sz = std::sqrt(elem.size()); // Size of elementary
    const int elem_bw = elem_sz * 2 - 1;

    const int stored_len = values[centr].size();

    for(int j = 0; j < stored_len - elem_sz+1; j++){
      int offset = 0;
      for(int e = 0; e < elem_sz; e++){
        for(int i = e; i < bandw-e; i++){
          values[i][j+e] += elem[offset + i];
        }
        offset += bandw - e*2 - 1; 
      }
    }
  }
  void add_identity_factor(double factor){
    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
    const int stored_len = values[centr].size();

    for(int i = 0; i < stored_len; i++){
        values[centr][i] += factor;
    }
  }

  void print(){
    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
    const int stored_len = values[centr].size();

    for(int i =0; i < bandw; i++){
        int l = stored_len - std::abs(i-centr);
        for(int j = 0; j < l; j++ ){
            std::cout<< values[i][j]<<"\t";
            if(j == hdr) std::cout<<":\t";
        }
        std::cout<<'\n';
    }
  }

  double get(int i, int j) const{
    if(j < hdr){
        return values[i][j];
    }else if(j > len - ftr - 1){
        return values[i][j - reps+1]; // End bit
    }else{
        return values[i][hdr];
    }
  }

  void pretty_print(){
    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
    for(int i = 0; i < bandw; i++){
      int l = len - std::abs(i-centr);
        for(int j = 0; j < l; j++){
            std::cout<<get(i, j)<<"\t";
        }
        std::cout<<'\n';
    }
  }


  banded_general<bandw> to_full(int len){
    auto tmp = banded_general<bandw>(len);
    for(int i = 0; i < bandw; i++){
        for(int j = 0; j < len; j++){
            if(j < hdr){
                tmp.values[i][j] = values[i][j];
            }else if(j >= len - ftr){
                tmp.values[i][j] = values[i][j - reps]; // End bit
            }else{
                tmp.values[i][j] = values[i][hdr];
            }
        }
    }
    
  }

};


using banded_general_penta = banded_general<5>;
using penta_repeating = banded_repeating<5>;

#endif
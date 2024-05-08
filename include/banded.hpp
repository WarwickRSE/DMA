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

  public:
  const int len;
  std::vector< std::vector<double> > values;

  banded_general(const int len_in):len(len_in){
    // Checks
    //len > bandw
    // bandw is odd
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

};

class banded_symm_repeating{

  const int len;
  const int bandw;
  const int hdr_id;
  const int hdr;

};


using banded_general_penta = banded_general<5>;

#endif
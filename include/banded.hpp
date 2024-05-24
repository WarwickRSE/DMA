// Compact storage of low-bandwidth diagonal matrices
// Customised for dynamics matrices: some number of boundary rows
// and a central block

#ifndef __bandedHPP__
#define __bandedHPP__

#include <vector>
#include <cmath>
#include <iostream>
#include <cassert>

namespace elementary_matrices{

    // Values. Bandwidth is (sz*2 -1); Simplify by not nesting...
    const std::vector<double> Spring{-1, 2, -1,   2, -4, 2,   -1, 2, -1};
    const std::vector<double> SpringBanded{-1, 2, -1, 2, -1,  2, -4, 2,  -1};
};

  // Tiny helper struct for banded indexing
  struct index_helper{
    int i, j;
    bool valid;
  };

class general_mat{
  // Convenience class for use in multiplications
  // Is neither performant nor generally useful!
  // Implements just enough for this use
  std::vector< std::vector<double> > values;
  public:
   const int len;
   general_mat(const int len_in):len(len_in){
    values.resize(len);
    for(int i = 0; i < len; ++i){
        values[i].resize(len);
    }
   }
    void clear(){
      for(int i = 0; i < len; i++){
        for(int j = 0; j < len; j++){
          values[i][j] = 0.0;
        }
      }
    }
    double get(int i, int j) const{
#ifdef DEBUG
      if(i >= len || j >= len) throw std::out_of_range("Out of range in get");
#endif
      return values[i][j];
    }
    void set(int i, int j, double val){
#ifdef DEBUG
      if(i >= len || j >= len) throw std::out_of_range("Out of range in set");
#endif
      values[i][j] = val;
    }
    friend int diag_bandw(const general_mat & mat){
      // Bandwith as in, for what N could this be stored as a N-diagonal matrix
      // Might be a simpler way but this is best guess for now
      int bw = 0;
      int centr = mat.len/2;
      // i is index of upper off diag
      for(int i=0; i<mat.len; i++){
        // j runs along the diagonal checking for non-zero in this diag or it's opposite
        for(int j = 0; j<mat.len-i; j++){
          if(mat.values[i][j] != 0.0 || mat.values[j][i] != 0.0){
            bw = std::max(bw, 2*(i-j)+1);
          }
        }
      }
      return bw; // By definition is either 0 or an odd number
    }

    void print()const{
      for(int i = 0; i < len; i++){
          for(int j = 0; j < len; j++){
              std::cout<<values[i][j]<<"\t";
          }
          std::cout<<'\n';
      }
    }
    const bool operator==(const general_mat & other) const{
      if(len != other.len) return false;
      for(int i = 0; i < len; i++){
        for(int j = 0; j < len; j++){
          if(values[i][j] != other.values[i][j]) return false;
        }
      }
      return true;
    }
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

  void multiply_all(double factor){
    for(int i = 0; i < bandw; i++){
        for(double & value : values[i]){
            value *= factor;
        }
    }
  }

  bool set_identity_row(int row){
    // Allow for any row
    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
    for(int i = 0; i < bandw; i++){
        int off = centr - i;
        if(off == 0){
          values[i][row] = 1.0;
        }
        else if(off > 0 && row-off >= 0){
            values[i][row-off] = 0.0;
        }else if(off < 0){
            values[i][row] = 0.0;
        }
    }
    return true;
  }

  void print()const{
    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
    for(int i =0; i < bandw; i++){
        int l = len - std::abs(i-centr);
        for(int j = 0; j < l; j++ ){
            std::cout<< values[i][j]<<"\t";
        }
        std::cout<<'\n';
    }
  }
  void pretty_print()const{print();}

  void expanded_print(bool full = false)const{
    // Print as full matrix
    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
    for(int i = 0; i < len; i++){
        for(int j = 0; j < len; j++){
            int col = j - i + centr;
            int off = centr - col;
            if( col >= centr) off = 0;
            if(col > -1 && col < bandw && i >= off){
               std::cout<<get(col, i-off)<<"\t";
            }else{
                std::cout<<"0\t";
            }
       }
        std::cout<<'\n';
    }
  }

  double get(int i, int j) const{
    // Get value in diagonal i, position j
#ifdef DEBUG
    if(i >= bandw || j >= len) throw std::out_of_range("Out of range banded get");
#endif
    return values[i][j];
  }
  void set(int i, int j, double val){
    // Set value in diagonal i, position j
#ifdef DEBUG
    if(i >= bandw || j >= len) throw std::out_of_range("Out of range banded set");
 #endif
    values[i][j] = val;
  }

  inline index_helper get_from_real(int i, int j)const{
    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
    index_helper inds;
    inds.i = j - i + centr;
    int off = i-j;
    inds.j = (off < 0 ? i : i-off);
    inds.valid = (inds.i >= 0 && inds.i < bandw && inds.j >= 0);
    return inds;
  }

  double get_real(int i, int j)const{
    // Get value at _real_ indices i, j
    index_helper inds = get_from_real(i, j);
    if(inds.valid){
      return get(inds.i, inds.j);
    }else{
      return 0.0;
    }
  }
  void set_real(int i, int j, double val){
    // Set value at _real_ indices i, j
    index_helper inds = get_from_real(i, j);
    if(inds.valid){
      set(inds.i, inds.j, val);
    }else{
#ifdef DEBUG
      assert(false, "Attempt to set non-existent element");
#endif
    }
  }

  general_mat to_general()const{
    // Expand into a normal matrix form
    auto tmp = general_mat(len);
    for(int i = 0; i < len; i++){
        for(int j = 0; j < len; j++){
          tmp.set(i, j, get_real(i, j));
        }
    }
    return tmp;
  }

  bool operator==(const banded_general & other) const{
    if(len != other.len) return false;
    for(int i = 0; i < bandw; i++){
      for(int j = 0; j < values[i].size(); j++){
        if(values[i][j] != other.values[i][j]) return false;
      }
    }
    return true;
  }
  friend bool operator==(const general_mat & gen, const banded_general & banded){
    if(gen.len != banded.len) return false;
    for(int i = 0; i < gen.len; i++){
      for(int j = 0; j < gen.len; j++){
        // Also checks for all out-of-band elements being zero because of getter implementation
        if(gen.get(i, j) != banded.get_real(i, j)) return false;
      }
    }
    return true;
  }
  friend bool operator==(const banded_general & banded, const general_mat & gen){
    return gen == banded;
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
  const int hdr; // Number of initial rows
  const int ftr; // Number of final rows

  inline int reps()const{return len - hdr - ftr;}

  banded_repeating(const int len_in, const int hdr_in, const int ftr_in):len(len_in), hdr(hdr_in), ftr(ftr_in){
    // Checks
#ifdef DEBUG
    assert(bandw % 2 == 1); // Bandwidth has to be odd (main diag plus 2*off-diag)
    assert(len > bandw); // Doesn't make any sense otherwise...
    assert(len > hdr+ftr+1); // Must have at least one repeating row
#endif 

    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
    values.resize(bandw);
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
  void change_len(int new_len){

#ifdef DEBUG
    assert(new_len > hdr+ftr+1); // Must have at least one repeating row
    assert(new_len > bandw); // Doesn't make any sense otherwise...
#endif

    len = new_len;
  }

  void fill_from_otherstuff(){ //TODO
  };


  void fill_from_elemBanded(const std::vector<double> elem){
    // ADDS the elementary banded array contributions
    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing

    const int elem_sz = std::sqrt(elem.size()); // Size of elementary
    const int elem_bw = elem_sz * 2 - 1;
    assert(elem_sz-1 <= hdr);

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
  void multiply_all(double factor){
    for(int i = 0; i < bandw; i++){
        for(double & value : values[i]){
            value *= factor;
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

  bool set_row(int row, std::vector<double> vals){
    if(row < 0 || row >= hdr) return false;
    values[row] = vals;
    return true;
  }
  bool set_identity_row(int row){
    // Only allowed to manipulate initial rows
    if(row < 0 || row >= hdr) return false;
    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
    for(int i = 0; i < bandw; i++){
        int off = centr - i;
        if(off == 0){
          values[i][row] = 1.0;
        }
        else if(off > 0 && row-off >= 0){
            values[i][row-off] = 0.0;
        }else if(off < 0){
            values[i][row] = 0.0;
        }
    }
    return true;
  }

  void print()const{
    for(int i =0; i < bandw; i++){
        for(int j = 0; j < values[i].size(); j++ ){
            std::cout<< values[i][j]<<"\t";
            if(j == hdr) std::cout<<":\t";
        }
        std::cout<<'\n';
    }
  }

  inline int get_second(int j)const{
    // Get the index for stored values corresponding to value j in actual diagonal
    return (j < hdr ? j :(j > len-ftr-1 ? (j-reps()+1) : hdr));
    /* Inline replication of:
     if(j < hdr){
        j;
    }else if(j > len - ftr - 1){
        j - reps+1; // End bit
    }else{
        hdr; // Repeating bit - note hdr is length of hdr section, so is correct index
    */
  }

  double get(int i, int j) const{
    // Return the value at position j on diagonal i. j can run from 0 to len - l where l is the number of the diagonal
#ifdef DEBUG
    if(i >= bandw || j >= len) throw std::out_of_range("Out of range");
#endif
    return values[i][get_second(j)];
  }

  inline index_helper index_from_real(int i, int j)const{
    // Indices suitable for passing to get - have to compress j if indexing directly!
    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
    index_helper inds;
    inds.i = j - i + centr;
    int off = i-j;
    inds.j = (off < 0 ? i : i-off);
    inds.valid = (inds.i >= 0 && inds.i < bandw && inds.j >= 0);
    return inds;
  }

  double get_real(int i, int j)const{
    // Get value at _real_ indices i, j
    index_helper inds = index_from_real(i, j);
    if(inds.valid){
      return get(inds.i, inds.j);
    }else{
      return 0.0;
    }
  }

  void pretty_print()const{
    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
    for(int i = 0; i < bandw; i++){
      int l = len - std::abs(i-centr);
        for(int j = 0; j < l; j++){
            std::cout<<get(i, j)<<"\t";
        }
        std::cout<<'\n';
    }
  }
  void expanded_print(bool full = false)const{
    // Print as full matrix
    // full flag means print all the repeating rows too
    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
    for(int i = 0; i < len; i++){
      if(full || i < hdr + 1 || i > len - ftr - 3){
        for(int j = 0; j < len; j++){
            std::cout<<get_real(i, j)<<"\t";
        }
        std::cout<<'\n';
      }else if(i == hdr + 1){
            std::cout<<"........\n";
      }
    }
  }

  banded_general<bandw> to_full()const{
    auto tmp = banded_general<bandw>(len);
    for(int i = 0; i < bandw; i++){
        for(int j = 0; j < len; j++){
          tmp.set(i, j, get(i, j));
        }
    }
    return tmp;
  }

  general_mat to_general()const{
    // Expand into a normal matrix form
    auto tmp = general_mat(len);
    for(int i = 0; i < len; i++){
        for(int j = 0; j < len; j++){
          tmp.set(i, j, get_real(i, j));
        }
    }
    return tmp;
  }

  friend bool operator==(const banded_repeating & lhs,  const banded_general<bandw> & rhs){
    // Check for equality without repeats
    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
    if(lhs.len != rhs.len) return false;
    for(int i = 0; i < bandw; i++){
      int l = lhs.len - std::abs(i-centr);
      for(int j = 0; j < l; j++ ){
        // Have to check every value in other across repeating section
        if(lhs.get(i, j) != rhs.get(i, j)) return false;
      }
    }
    return true;
  }
  friend bool operator==(const banded_general<bandw> & lhs, const banded_repeating & rhs){
    return rhs == lhs;
  }
  friend bool operator==(const general_mat & gen, const banded_repeating & banded){
    if(gen.len != banded.len) return false;
    for(int i = 0; i < gen.len; i++){
      for(int j = 0; j < gen.len; j++){
        if(gen.get(i, j) != banded.get_real(i, j)) return false;
      }
    }
    return true;
  }
  friend bool operator==(const banded_repeating & banded, const general_mat & gen){
    return gen == banded;
  }

  bool operator==(const banded_repeating & other) const{
    if(len != other.len) return false;
    for(int i = 0; i < bandw; i++){
      for(int j = 0; j < values[i].size(); j++){
        if(values[i][j] != other.values[i][j]) return false;
      }
    }
    return true;
  }
};

general_mat matmul(const general_mat & L, const general_mat & U){
  assert(L.len == U.len);
  general_mat result(L.len);
  for(int i = 0; i < L.len; i++){
    for(int j = 0; j < L.len; j++){
      double tmp = 0.0;
      for(int k = 0; k < L.len; k++){
        tmp += L.get(i, k)*U.get(k, j);
      }
      result.set(i, j, tmp);
    }
  }
  return result;
}

template <int bandw>
general_mat matmul(const banded_general<bandw> & L, const banded_general<bandw> & U){
  // Multiply two banded matrices. Note that in general the result is NOT banded, so only a general matrix is returned
  assert(L.len == U.len);
  general_mat result(L.len);
  for(int i = 0; i < L.len; i++){
    for(int j = 0; j < L.len; j++){
      double tmp = 0.0;
      for(int k = 0; k < L.len; k++){
        tmp += L.get_real(i, k)*U.get_real(k, j);
      }
      result.set(i, j, tmp);
    }
  }
  return result;
}

template <int bandw>
general_mat matmul(const banded_repeating<bandw> & L, const banded_repeating<bandw> & U){
  // Multiply two repeating matrices. Note that in general the result is NOT banded, so only a general matrix is returned
  assert(L.len == U.len);
  general_mat result(L.len);
  for(int i = 0; i < L.len; i++){
    for(int j = 0; j < L.len; j++){
      double tmp = 0.0;
      for(int k = 0; k < L.len; k++){
        tmp += L.get_real(i, k)*U.get_real(k, j);
      }
      result.set(i, j, tmp);
    }
  }
  return result;
}

template <int bandw>
banded_general<bandw> from_general(const general_mat & gen){
  // Form a banded matrix from a general one, IFF bandw is sufficient
  int bandw_in = diag_bandw(gen);
  assert(bandw_in <= bandw); // Check it's storable!
  banded_general<bandw> res(gen.len);
  for(int i = 0; i < gen.len; i++){
    for(int j = 0; j < gen.len; j++){
      res.set_real(i, j, gen.get(i, j));
    }
  }
  return res;
}

using banded_general_penta = banded_general<5>;
using penta_repeating = banded_repeating<5>;

#endif
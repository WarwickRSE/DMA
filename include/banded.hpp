// Compact storage of low-bandwidth diagonal matrices
// Customised for dynamics matrices: some number of boundary rows
// and a central block

#ifndef __bandedHPP__
#define __bandedHPP__

#include <vector>
#include <cmath>
#include <iostream>
#include <cassert>
#include <type_traits>

namespace elementary_matrices{

    // Values. Bandwidth is (sz*2 -1); Simplify by not nesting...
    const std::vector<double> Spring{-1, 2, -1,   2, -4, 2,   -1, 2, -1};
    const std::vector<double> SpringBanded{-1, 2, -1, 2, -1,  2, -4, 2,  -1};
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

template <int bandw, bool repeating=false>
class banded_general{

  const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing

  typedef typename std::conditional<repeating, int, const int>::type int_c;
  typedef typename std::conditional<repeating, const int, int>::type int_c2;

  std::vector< std::vector<double> > values;

  inline int reps()const{return len - hdr - ftr;}

  public:
  int_c len; // Variable for repeating, const for not
  int_c2 hdr, ftr; // Const for repeating, redundant for not

  banded_general(const int len_in, const int hdr_in=0, const int ftr_in=0):len(len_in), hdr(hdr_in), ftr(ftr_in){
  // Checks
#ifdef DEBUG
    static_assert(bandw % 2 == 1); // Bandwidth has to be odd (main diag plus 2*off-diag)
    assert(len > bandw); // Doesn't make any sense otherwise...
    assert(!repeating || len > hdr+ftr+1); // Must have at least one repeating row if repeating
#endif 

    values.resize(bandw);
    for(int i = 0; i < bandw; ++i){
        int l;
        if constexpr(repeating){
          l = hdr + ftr + 1;
        }else{
          l = len;
        }
        l -= std::abs(i-centr);
        values[i].resize(l);
    }
  }
  banded_general(const banded_general & other):len(other.len), hdr(other.hdr), ftr(other.ftr){
    values = other.values;
  }
  banded_general & operator=(const banded_general & other){
    // Types must be conformable - for non-repeating that means same len
    // For repeating it means hdr and ftr same
    if constexpr(repeating){
      assert(hdr == other.hdr && ftr == other.ftr);
    }else{
      assert(len == other.len);
    }
    if(this != &other){
      values = other.values;
    }
    return *this;
  }

  void clear(){
    for(int i = 0; i < bandw; i++){
        for(double & value : values[i]){
            value = 0.0;
        }
    }
  }
  void change_len(int new_len){
    // Checked way to change the length. For unchecked, just do it
    if constexpr(!repeating){
      static_assert(repeating, "Cannot change length of non-repeating matrix");
    }
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
    const int elem_sz = std::sqrt(elem.size()); // Size of elementary
    const int elem_bw = elem_sz * 2 - 1;
    assert(!repeating || elem_sz-1 <= hdr);// Else the elementary matrix will not repeat correctly

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
    for(int i = 0; i < values[centr].size(); i++){
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

  void set_identity_row(int row){
    // For repeating, only the top rows; for non repeating any row
#ifdef DEBUG
    assert(row >= 0 || row < len);
    assert(!repeating || row < hdr);
#endif
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
  }

  void print()const{
    for(int i =0; i < bandw; i++){
        for(int j = 0; j < values[i].size(); j++ ){
            std::cout<< values[i][j]<<"\t";
            if(repeating && j == hdr) std::cout<<":\t";
        }
        std::cout<<'\n';
    }
  }

  void pretty_print()const{
    if constexpr(!repeating){
      print();
    }else{
      for(int i = 0; i < bandw; i++){
        int l = len - std::abs(i-centr); // Real diagonal length, not stored
          for(int j = 0; j < l; j++){
              std::cout<<get(i, j)<<"\t";
          }
          std::cout<<'\n';
      }
    }
  }

  void expanded_print(bool full = false)const{
    // Print as full matrix
    for(int i = 0; i < len; i++){
      bool print_row = !repeating || (full || i < hdr + 1 || i > len - ftr - 3);
      if(print_row){
        for(int j = 0; j < len; j++){
          std::cout<<get_real(i, j)<<"\t";
        }
        std::cout<<'\n';
      }else if(repeating && i == hdr + 1){
        std::cout<<"........\n";
      }
    }
  }

  double get(int i, int j) const{
    // Get value in diagonal i, position j. j can run from 0 to len - l where l is the number of the diagonal
#ifdef DEBUG
    const int l = len - std::abs(i-centr);
    if(i >= bandw || j >= l || i < 0 || j < 0) throw std::out_of_range("Out of range banded get");
#endif
    if constexpr(repeating){
      return values[i][(j < hdr ? j :(j > len-ftr-1 ? (j-reps()+1) : hdr))];
    }else{
      return values[i][j];
    }
  }
  /* NOTE: that horror show is the inline replication of:
     if(j < hdr){
        j;
    }else if(j > len - ftr - 1){
        j - reps+1; // End bit
    }else{
        hdr; // Repeating bit - note hdr is length of hdr section, so is correct index
    */

  // Tiny helper struct for banded indexing
  struct index_helper{
    int i, j;
    bool valid;
  };
  // Convert from real indices to diagonal ones
  inline index_helper index_from_real(int i, int j)const{
    index_helper inds;
    inds.i = j - i + centr;
    inds.j = (i < j ? i : j);
    inds.valid = (inds.i >= 0 && inds.i < bandw && inds.j >= 0 && inds.j < len- std::abs(j-i));
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

// Set values only for non-repeating matrices
  void set(int i, int j, double val){
    if constexpr(repeating){
      static_assert(!repeating, "Cannot set values in repeating matrix");
    }
    // Set value in diagonal i, position j
#ifdef DEBUG
    const int l = len - std::abs(i-centr);
    if(i >= bandw || j >= l || i < 0 || j < 0) throw std::out_of_range("Out of range banded set");
 #endif
    values[i][j] = val;
  }

  void set_real(int i, int j, double val){
    if constexpr(repeating){
      static_assert(!repeating, "Cannot set values in repeating matrix");
    }
    // Set value at _real_ indices i, j
    index_helper inds = index_from_real<bandw>(len, i, j);
    if(inds.valid){
      set(inds.i, inds.j, val);
    }else{
#ifdef DEBUG
      throw std::out_of_range("Attempt to set non-existent element");
#endif
    }
  }

   banded_general<bandw, false> to_full()const{
    if constexpr(!repeating){
      // No-op
      return *this;
    }
    auto tmp = banded_general<bandw, false>(len);
    for(int i = 0; i < bandw; i++){
        const int l = len - std::abs(i-centr);
        for(int j = 0; j < l; j++){
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

  friend bool operator==(const banded_general & lhs, const banded_general & rhs){
    if(lhs.len != rhs.len) return false;
    for(int i = 0; i < bandw; i++){
      // Only need to check the stored values
      for(int j = 0; j < lhs.values[i].size(); j++){
        if(lhs.values[i][j] != rhs.values[i][j]) return false;
      }
    }
    return true;
  }
  //Comparisons with other value for repeating...
  friend bool operator==(const banded_general<bandw, repeating> & lhs,  const banded_general<bandw, !repeating> & rhs){
    // Check for equality without repeats
    const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
    if(lhs.len != rhs.len) return false;
    for(int i = 0; i < bandw; i++){
      int l = lhs.len - std::abs(i-centr);
      for(int j = 0; j < l; j++ ){
        // Have to check every value in repeating section
        if(lhs.get(i, j) != rhs.get(i, j)) return false;
      }
    }
    return true;
  }

  //Equality with full matrices
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

template <int bandw, bool repeatinga, bool repeatingb>
general_mat matmul(const banded_general<bandw, repeatinga> & L, const banded_general<bandw, repeatingb> & U){
  // Multiply two banded matrices. Note that in general the result is NOT banded, so only a general matrix is returned
  assert(L.len == U.len);
  const int centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
  general_mat result(L.len);
  for(int i = 0; i < L.len; i++){
    for(int j = 0; j < L.len; j++){
      double tmp = 0.0;
      for(int k = i-centr; k < i+centr+1; k++){
        tmp += L.get_real(i, k)*U.get_real(k, j);
      }
      result.set(i, j, tmp);
    }
  }
  return result;
}

template <int bandw>
banded_general<bandw, false> from_general(const general_mat & gen){
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

using banded_general_penta = banded_general<5, false>;
using penta_repeating = banded_general<5, true>;

#endif
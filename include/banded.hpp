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
#include <algorithm>

namespace elementary_matrices{

    // Values. Bandwidth is (sz*2 -1); Simplify by not nesting...
    const std::vector<double> TriSpring{-1, 2, -1};
    const std::vector<double> Spring{-1, 2, -1,   2, -4, 2,   -1, 2, -1};
    const std::vector<double> SpringBanded{-1, 2, -1, 2, -1,  2, -4, 2,  -1};
};

class general_mat{
  // Convenience class for use in multiplications
  // Is neither performant nor generally useful!
  // Implements just enough for this use
  std::vector< std::vector<double> > values;
  public:
   const size_t len;
   general_mat(const size_t len_in):len(len_in){
    values.resize(len);
    for(size_t i = 0; i < len; ++i){
        values[i].resize(len);
    }
   }
    void clear(){
      for(size_t i = 0; i < len; i++){
        for(size_t j = 0; j < len; j++){
          values[i][j] = 0.0;
        }
      }
    }
    double get(size_t i, size_t j) const{
#ifdef DEBUG
      if(i >= len || j >= len) throw std::out_of_range("Out of range in get");
#endif
      return values[i][j];
    }
    void set(size_t i, size_t j, double val){
#ifdef DEBUG
      if(i >= len || j >= len) throw std::out_of_range("Out of range in set");
#endif
      values[i][j] = val;
    }
    friend size_t diag_bandw(const general_mat & mat){
      // Bandwith as in, for what N could this be stored as a N-diagonal matrix
      // Might be a simpler way but this is best guess for now
      size_t bw = 0;
      // i is index of upper off diag
      for(size_t i=0; i<mat.len; i++){
        // j runs along the diagonal checking for non-zero in this diag or it's opposite
        for(size_t j = 0; j<mat.len-i; j++){
          if(mat.values[i][j] != 0.0 || mat.values[j][i] != 0.0){
            bw = std::max(bw, 2*(i-j)+1);
          }
        }
      }
      return bw; // By definition is either 0 or an odd number
    }

    void print()const{
      for(size_t i = 0; i < len; i++){
          for(size_t j = 0; j < len; j++){
              std::cout<<values[i][j]<<"\t";
          }
          std::cout<<'\n';
      }
    }
    bool operator==(const general_mat & other) const{
      if(len != other.len) return false;
      for(size_t i = 0; i < len; i++){
        for(size_t j = 0; j < len; j++){
          if(values[i][j] != other.values[i][j]) return false;
        }
      }
      return true;
    }
    bool operator!=(const general_mat & other) const{
      return !(*this == other);
    }
};

template <size_t bandw, bool repeating=false>
class banded_general{

  const size_t centr = (bandw + 1)/2 - 1; // -1 for 0-indexing

  typedef typename std::conditional<repeating, size_t, const size_t>::type int_c;
  typedef typename std::conditional<repeating, const size_t, size_t>::type int_c2;

  std::vector< std::vector<double> > values;

  inline size_t reps()const{return len - hdr - ftr;}

  public:
  int_c len; // Variable for repeating, const for not
  int_c2 hdr, ftr; // Const for repeating, redundant for not

  inline size_t st_len(size_t i)const{return len - std::abs((long)(i-centr));}

  banded_general(size_t len_in, size_t hdr_in=0, size_t ftr_in=0):len(len_in), hdr(hdr_in), ftr(ftr_in){
  // Checks
#ifdef DEBUG
    static_assert(bandw % 2 == 1); // Bandwidth has to be odd (main diag plus 2*off-diag)
    assert(len >= bandw); // Doesn't make any sense otherwise...
    assert(!repeating || len >= hdr+ftr+1); // Must have at least one repeating row if repeating
#endif 

    values.resize(bandw);
    for(size_t i = 0; i < bandw; ++i){
        size_t l;
        if constexpr(repeating){
          l = hdr + ftr + 1;
        }else{
          l = len;
        }
        l -= std::abs((long)(i-centr));
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
    for(size_t i = 0; i < bandw; i++){
        for(double & value : values[i]){
            value = 0.0;
        }
    }
  }
  void change_len(size_t new_len){
    // Checked way to change the length. For unchecked, just do it
    if constexpr(!repeating){
      static_assert(repeating, "Cannot change length of non-repeating matrix");
    }
#ifdef DEBUG
    assert(new_len >= hdr+ftr+1); // Must have at least one repeating row
    assert(new_len >= bandw); // Doesn't make any sense otherwise...
#endif
    len = new_len;
  }

  void fill_from_otherstuff(){ //TODO
  };


  void fill_from_elemBanded(const std::vector<double> elem){
    // ADDS the elementary banded array contributions
    const size_t elem_sz = std::sqrt(elem.size()); // Size of elementary
    assert(!repeating || elem_sz-1 <= hdr);// Else the elementary matrix will not repeat correctly

    const size_t stored_len = values[centr].size();
    for(size_t j = 0; j < stored_len - elem_sz+1; j++){
      int offset = 0;
      for(size_t e = 0; e < elem_sz; e++){
        for(size_t i = e; i < bandw-e; i++){
          values[i][j+e] += elem[offset + i];
        }
        offset += bandw - e*2 - 1; 
      }
    }
  }
  void add_identity_factor(double factor){
    for(size_t i = 0; i < values[centr].size(); i++){
        values[centr][i] += factor;
    }
  }

  void multiply_all(double factor){
    for(size_t i = 0; i < bandw; i++){
        for(double & value : values[i]){
            value *= factor;
        }
    }
  }

  void set_identity_row(size_t row){
    // For repeating, only the top rows; for non repeating any row
#ifdef DEBUG
    assert(row < len);
    assert(!repeating || row < hdr);
#endif
    for(size_t i = 0; i < bandw; i++){
        int off = centr - i;
        if(off == 0){
          values[i][row] = 1.0;
        }
        else if(off > 0 && row >= off){
            values[i][row-off] = 0.0;
        }else if(off < 0){
            values[i][row] = 0.0;
        }
    }
  }

  void print()const{
    for(size_t i =0; i < bandw; i++){
        for(size_t j = 0; j < values[i].size(); j++ ){
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
      for(size_t i = 0; i < bandw; i++){
          for(size_t j = 0; j < st_len(i); j++){
              std::cout<<get(i, j)<<"\t";
          }
          std::cout<<'\n';
      }
    }
  }

  void expanded_print(bool full = false)const{
    // Print as full matrix
    for(size_t i = 0; i < len; i++){
      bool print_row = !repeating || (full || i < hdr + 1 || i > len - ftr - 3);
      if(print_row){
        for(size_t j = 0; j < len; j++){
          std::cout<<get_real(i, j)<<"\t";
        }
        std::cout<<'\n';
      }else if(repeating && i == hdr + 1){
        std::cout<<"........\n";
      }
    }
  }

  double get(size_t i, size_t j) const{
    // Get value in diagonal i, position j. j can run from 0 to len - l where l is the number of the diagonal
#ifdef DEBUG
    if(i >= bandw || j >= st_len(i)) throw std::out_of_range("Out of range banded get");
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
    size_t i, j;
    bool valid;
  };
  // Convert from real indices to diagonal ones
  inline index_helper index_from_real(size_t i, size_t j)const{
    index_helper inds;
    inds.i = j - i + centr;
    inds.j = (i < j ? i : j);
    inds.valid = (j+centr >= i && inds.i < bandw && inds.j < st_len(inds.i));
    return inds;
  }

  double get_real(size_t i, size_t j)const{
    // Get value at _real_ indices i, j
    index_helper inds = index_from_real(i, j);
    if(inds.valid){
      return get(inds.i, inds.j);
    }else{
      return 0.0;
    }
  }

// Set values only for non-repeating matrices
  void set(size_t i, size_t j, double val){
    if constexpr(repeating){
      static_assert(!repeating, "Cannot set values in repeating matrix");
    }
    // Set value in diagonal i, position j
#ifdef DEBUG
    if(i >= bandw || j >= st_len(i)) throw std::out_of_range("Out of range banded set");
 #endif
    values[i][j] = val;
  }

  void set_real(size_t i, size_t j, double val){
    if constexpr(repeating){
      static_assert(!repeating, "Cannot set values in repeating matrix");
    }
    // Set value at _real_ indices i, j
    index_helper inds = index_from_real(len, i, j);
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
    for(size_t i = 0; i < bandw; i++){
        for(size_t j = 0; j < st_len(i); j++){
          tmp.set(i, j, get(i, j));
        }
    }
    return tmp;
  }

  general_mat to_general()const{
    // Expand into a normal matrix form
    auto tmp = general_mat(len);
    for(size_t i = 0; i < len; i++){
        for(size_t j = 0; j < len; j++){
          tmp.set(i, j, get_real(i, j));
        }
    }
    return tmp;
  }

  friend bool operator==(const banded_general & lhs, const banded_general & rhs){
    if(lhs.len != rhs.len) return false;
    for(size_t i = 0; i < bandw; i++){
      // Only need to check the stored values
      for(size_t j = 0; j < lhs.values[i].size(); j++){
        if(lhs.values[i][j] != rhs.values[i][j]) return false;
      }
    }
    return true;
  }
  friend bool operator!=(const banded_general & lhs, const banded_general & rhs){
    return !(lhs == rhs);
  }
  //Comparisons with other value for repeating...
  friend bool operator==(const banded_general<bandw, repeating> & lhs,  const banded_general<bandw, !repeating> & rhs){
    // Check for equality without repeats
    if(lhs.len != rhs.len) return false;
    for(size_t i = 0; i < bandw; i++){
      size_t l = lhs.st_len(i);
      for(size_t j = 0; j < l; j++ ){
        // Have to check every value in repeating section
        if(lhs.get(i, j) != rhs.get(i, j)) return false;
      }
    }
    return true;
  }
  friend bool operator!=(const banded_general<bandw, repeating> & lhs,  const banded_general<bandw, !repeating> & rhs){
    return !(lhs == rhs);
  }

  //Equality with full matrices
  friend bool operator==(const general_mat & gen, const banded_general & banded){
    if(gen.len != banded.len) return false;
    for(size_t i = 0; i < gen.len; i++){
      for(size_t j = 0; j < gen.len; j++){
        // Also checks for all out-of-band elements being zero because of getter implementation
        if(gen.get(i, j) != banded.get_real(i, j)) return false;
      }
    }
    return true;
  }
  friend bool operator==(const banded_general & banded, const general_mat & gen){
    return gen == banded;
  }
  friend bool operator!=(const general_mat & gen, const banded_general & banded){
    return !(gen == banded);
  }
  friend bool operator!=(const banded_general & banded, const general_mat & gen){
    return !(gen == banded);
  }

};

inline general_mat matmul(const general_mat & L, const general_mat & U){
  assert(L.len == U.len);
  general_mat result(L.len);
  for(size_t i = 0; i < L.len; i++){
    for(size_t j = 0; j < L.len; j++){
      double tmp = 0.0;
      for(size_t k = 0; k < L.len; k++){
        tmp += L.get(i, k)*U.get(k, j);
      }
      result.set(i, j, tmp);
    }
  }
  return result;
}

inline std::vector<double> matvecmult(const general_mat & A, std::vector<double> & x){
  // Calculate A x for matrix A and vector x
  assert(A.len == x.size());
  std::vector<double> b(A.len);
  for(size_t i = 0; i < A.len; i++){
    double tmp = 0.0;
    for(size_t j = 0; j < A.len; j++){
      tmp += A.get(i, j)*x[j];
    }
    b[i] = tmp;
  }
  return b;
}

template <size_t bandw, bool repeatinga, bool repeatingb>
general_mat matmul(const banded_general<bandw, repeatinga> & L, const banded_general<bandw, repeatingb> & U){
  // Multiply two banded matrices. Note that in general the result is NOT banded, so only a general matrix is returned
  assert(L.len == U.len);
  const size_t centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
  general_mat result(L.len);
  for(size_t i = 0; i < L.len; i++){
    for(size_t j = 0; j < L.len; j++){
      double tmp = 0.0;
      for(size_t k = std::max(0, (int)(i-centr)); k < std::min(i+centr+1, L.len); k++){
        tmp += L.get_real(i, k)*U.get_real(k, j);
      }
      result.set(i, j, tmp);
    }
  }
  return result;
}
template <size_t bandw, bool repeating>
inline std::vector<double> matvecmult(const banded_general<bandw, repeating> & A, std::vector<double> & x){
  // Calculate A x for matrix A and vector x
  assert(A.len == x.size());
  const size_t centr = (bandw + 1)/2 - 1; // -1 for 0-indexing
  std::vector<double> b(A.len);
  for(size_t i = 0; i < A.len; i++){
    double tmp = 0.0;
    for(size_t j = std::max(0, (int)(i-centr)); j < std::min(i+centr+1, A.len); j++){
      tmp += A.get_real(i, j)*x[j];
    }
    b[i] = tmp;
  }
  return b;
}

template <size_t bandw>
banded_general<bandw, false> from_general(const general_mat & gen){
  // Form a banded matrix from a general one, IFF bandw is sufficient
  size_t bandw_in = diag_bandw(gen);
  assert(bandw_in <= bandw); // Check it's storable!
  banded_general<bandw> res(gen.len);
  for(size_t i = 0; i < gen.len; i++){
    for(size_t j = 0; j < gen.len; j++){
      res.set_real(i, j, gen.get(i, j));
    }
  }
  return res;
}

using banded_general_penta = banded_general<5, false>;
using penta_repeating = banded_general<5, true>;
using banded_general_tri = banded_general<3, false>;
using tri_repeating = banded_general<3, true>;

#endif
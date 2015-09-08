#ifndef MATRIX_CONV_H
#define MATRIX_CONV_H

#include <vector>

#include <assert.h>

#include <Eigen/Dense>
#include <zjucad/matrix/matrix.h>

namespace zsw{
  template<typename VAL_TYPE>
    void vvec3ToZjumat(const std::vector<Eigen::Matrix<VAL_TYPE, 3, 1> > &vvec3,
                       zjucad::matrix::matrix<VAL_TYPE> &zju_mat)
    {
      int cols = vvec3.size();
      assert(cols > 0);
      zju_mat.resize(3,cols);

      for (int i=0; i<cols; ++i) {
        zju_mat(0,i) = vvec3[i](0);
        zju_mat(1,i) = vvec3[i](1);
        zju_mat(2,i) = vvec3[i](2);
      }
    }

  template<typename VAL_TYPE>
    void zjumat3ToVvec(const zjucad::matrix::matrix<VAL_TYPE> &zju_mat,
                       std::vector<Eigen::Matrix<VAL_TYPE, 3, 1> > &vvec3)
    {
      assert(zju_mat.size(1) == 3);
      int cols = zju_mat.size(2);

      vvec3.reserve(cols); vvec3.clear();
      for (int i=0; i<cols; ++i) {
        Eigen::Matrix<VAL_TYPE, 3, 1> tmp_v(zju_mat(0,i), zju_mat(1,i),zju_mat(2,i));
        vvec3.push_back(tmp_v);
      }
    }

  template<typename VAL_TYPE>
    void vvec4ToZjumat(const std::vector<Eigen::Matrix<VAL_TYPE, 4, 1> > &vvec,
                       zjucad::matrix::matrix<VAL_TYPE> &zju_mat)
    {
      int cols = vvec.size();
      assert(cols > 0);
      zju_mat.resize(4, cols);
      for (int i=0; i<cols; ++i) {
        zju_mat(0,i) = vvec[i](0);
        zju_mat(1,i) = vvec[i](1);
        zju_mat(2,i) = vvec[i](2);
        zju_mat(3,i) = vvec[i](3);
      }
    }

  template<typename VAL_TYPE>
    void zjumat4ToVvec(const zjucad::matrix::matrix<VAL_TYPE> &zju_mat,
                       std::vector<Eigen::Matrix<VAL_TYPE, 4, 1> > &vvec4)
    {
      assert(zju_mat.size(1) == 4);
      int cols = zju_mat.size(2);
      vvec4.reserve(cols); vvec4.clear();

      for (int i=0; i<cols; ++i) {
        Eigen::Matrix<VAL_TYPE, 4, 1> tmp_v(zju_mat(0, i), zju_mat(1,i), zju_mat(2,i), zju_mat(3,i));
        vvec4.push_back(tmp_v);
      }
    }

  template<typename VAL_TYPE, typename INT_TYPE>
    void eigenSparse2hjSparse(const Eigen::SparseMatrix<VAL_TYPE> &A, hj::sparse::csc<VAL_TYPE,INT_TYPE> &B)
  {
    A.makeCompressed();
    B.resize(A.rows(), A.cols(), A.nonZeros());

    std::copy(A.innerIndexPtr(), A.innerIndexPtr()+A.nonZeros(), B.idx().begin());
    std::copy(A.outerIndexPtr(), A.outIndexPtr()+A.cols()+1, B.ptr().begin());
    std::copy(A.valuePtr(), A.valuePtr()+A.nonZeros(), B.val().begin());
  }

  template<typename VAL_TYPE, typename INT_TYPE>
    void hjSparse2eigenSparse(const hj::sparse::csc<VAL_TYPE, INT_TYPE> &A, Eigen::SparseMatrix<VAL_TYPE> &B)
  {
    std::vector<Eigen::Triplet<VAL_TYPE> > vec(nnz(A));
    int idx = 0;
    for (size_t j=0; j<A.size(2); ++j) {
      for (size_t cnt = A.ptr()[j]; cnt<A.ptr()[j+1]; ++cnt) {
        vec[idx++] = Eigen::Triplet<VAL_TYPE>(A.idx()[cnt], j, A.val()[cnt]);
      }
    }
    B.resize(A.size(1), A.size(2));
    B.reserve(nnz(A));
    B.setFromTriplets(vec.begin(), vec.end());
    B.makeCompressed();
  }

  template<typename VAL_TYPE>
  void vec2Zjumat3n(const std::vector<VAL_TYPE> &vec, zjucad::matrix::matrix<VAL_TYPE> &zju_mat)
  {
    assert(vec.size()%3==0);
    zju_mat.resize(3, vec.size()/3);
    std::copy(vec.begin(), vec.end(), zju_mat.begin());
  }

  template<typename VAL_TYPE>
  void zjumat3n2Vec(const zjucad::matrix::matrix<VAL_TYPE> &zju_mat, std::vector<VAL_TYPE> &vec)
  {
      vec.resize(zju_mat.size(1)*zju_mat.size(2));
      std::copy(zju_mat.begin(), zju_mat.end(), vec.begin());
  }
}

#endif /* MATRIX_CONV_H */

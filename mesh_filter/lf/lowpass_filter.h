#ifndef LOWPASS_FILTER_H
#define LOWPASS_FILTER_H

#include <armadillo>
#include "../data.h"

namespace zsw
{
  namespace mesh {
    class ZSW_API LowPassFilter {
    public:
      LowPassFilter(std::shared_ptr<zsw::data::FilterDataInterFace> df);
      void filter();
      void run();
    private:
      void calcOpMat(std::vector<zsw::FakeSet<size_t>> &ring, Eigen::Matrix<zsw::Scalar, Eigen::Dynamic, Eigen::Dynamic> &mat);
      void calcOpMat(std::vector<zsw::FakeSet<size_t>> &ring, Eigen::SparseMatrix<zsw::Scalar> &mat);
      Eigen::Matrix<zsw::Scalar, Eigen::Dynamic, Eigen::Dynamic> mat_l_;
      Eigen::SparseMatrix<zsw::Scalar> smat_l_;
      Eigen::Matrix<zsw::Scalar, Eigen::Dynamic, 3> x_h_;
      Eigen::Matrix<zsw::Scalar, Eigen::Dynamic, 3> p_;
      Eigen::Matrix<zsw::Scalar, Eigen::Dynamic, 3> res_p_;
      std::shared_ptr<zsw::data::FilterDataInterFace> df_;
    };
  }
  template<typename Scalar>
    void eigen2arma(const Eigen::SparseMatrix<Scalar> &emat, arma::SpMat<Scalar> &amat)
    {
      size_t row = emat.rows();
      size_t col = emat.cols();
      // amat.resize(row,col);
      for(size_t i=0; i<row; ++i) {
        for(size_t j=0; j<col; ++j) {
          if(emat.coeff(i,j)!=0) {
            amat(i,j) = emat.coeff(i,j);
          }
        }
      }
    }
  template<typename Scalar>
    void arma2eigen(const arma::vec avec, Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& evec)
    {
      evec.resize(avec.size());
      std::copy(avec.begin(), avec.end(), evec.begin());
    }
}

#endif /* LOWPASS_FILTER_H */

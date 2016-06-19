#ifndef LH_FILTER_H
#define LH_FILTER_H

#include "../data.h"

namespace zsw {
  namespace mesh {
    class ZSW_API LHFilter
    {
    public:
      LHFilter(std::shared_ptr<zsw::data::LHFDataInterFace> df);
      int run();
    private:
      void calcOpMat(std::vector<zsw::FakeSet<size_t>> &ring, Eigen::SparseMatrix<zsw::Scalar> &mat);
      void buildEquation();
      int solveEquation();
      Eigen::SparseMatrix<zsw::Scalar> mat_l_;
      Eigen::SparseMatrix<zsw::Scalar> mat_h_;
      Eigen::SparseMatrix<zsw::Scalar> A_;
      Eigen::Matrix<zsw::Scalar, Eigen::Dynamic, 1> b_;
      Eigen::Matrix<zsw::Scalar, Eigen::Dynamic, 1> p_;
      Eigen::Matrix<zsw::Scalar, Eigen::Dynamic, 1> res_p_;
      std::shared_ptr<zsw::data::LHFDataInterFace> df_;
    };
  }
}


#endif /* LH_FILTER_H */

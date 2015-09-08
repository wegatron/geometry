#ifndef MATRIX_TYPE_H
#define MATRIX_TYPE_H

#include <Eigen/Dense>
#include <zswlib/config.h>

namespace zsw {
  typedef Eigen::Matrix<zsw::Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;
  typedef Eigen::Matrix<zsw::Scalar, Eigen::Dynamic, 1> VectorXs;
  typedef Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXi;
  typedef Eigen::Matrix<size_t, Eigen::Dynamic, 1> VectorXi; 
};


#endif /* MATRIX_TYPE_H */

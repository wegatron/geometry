#ifndef SMOOTHER_H
#define SMOOTHER_H

#include <Eigen/Dense>
#include <zswlib/config.h>

namespace zsw{
  void laplaceSmooth(Eigen::Matrix<zsw::Scalar,3,1> &pt, const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &one_ring,
                     const Eigen::Matrix<zsw::Scalar,3,1> ori_pt, const zsw::Scalar dis_limit);

}
#endif /* SMOOTHER_H */

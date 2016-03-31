#include "smoother.h"
#include <iostream>
#include <zswlib/const_val.h>

void zsw::laplaceSmooth(Eigen::Matrix<zsw::Scalar,3,1> &pt, const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &one_ring,
                        const Eigen::Matrix<zsw::Scalar,3,1> &ori_pt, const size_t n, const zsw::Scalar dis_limit)
{
  Eigen::Matrix<zsw::Scalar,3,1> mean_pt=Eigen::Matrix<zsw::Scalar,3,1>::Zero();
  for(const Eigen::Matrix<zsw::Scalar,3,1> &pt : one_ring) {    mean_pt+=pt;  }
  Eigen::Matrix<zsw::Scalar,3,1> vn = 1.0/n * (1.0/one_ring.size() * mean_pt - pt);
  // add weight
  Eigen::Matrix<zsw::Scalar,3,1> m = pt-ori_pt;
  zsw::Scalar a=vn.squaredNorm();
  zsw::Scalar b=2*vn.dot(m);
  zsw::Scalar c=m.squaredNorm()-dis_limit*dis_limit;
  zsw::Scalar t=(-b+sqrt(b*b-4*a*c))/(2*a);
  assert(fabs((vn*t+pt-ori_pt).norm()-dis_limit) < zsw::const_val::eps);
  t=std::min(1.0, t);
  pt=vn*t+pt;
}

#include "sampling.h"

#include <iostream>
#include <Eigen/Geometry>
#include <zswlib/const_val.h>

void zsw::sampleTriangle(const Eigen::Matrix<zsw::Scalar, 3, 3> &tri_points, const zsw::Scalar r,
                         std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &samples)
{
  Eigen::Matrix<zsw::Scalar,3,1> v0=tri_points.block<3,1>(0,0);
  Eigen::Matrix<zsw::Scalar,3,1> v1=tri_points.block<3,1>(0,1);
  Eigen::Matrix<zsw::Scalar,3,1> v2=tri_points.block<3,1>(0,2);
  Eigen::Matrix<zsw::Scalar,3,1> n=v1-v0;
  Eigen::Matrix<zsw::Scalar,3,1> m=v2-v0;
  zsw::Scalar step_n=r/n.norm();
  zsw::Scalar step_m=r/m.norm();
  bool sn_flg=true;
  for(zsw::Scalar sn=0; sn_flg; sn+=step_n) {
    if(sn>1) { sn=1; sn_flg=false; }
    const Eigen::Matrix<zsw::Scalar,3,1> tmp_n=sn*n;
    const zsw::Scalar max_sm=1-sn-zsw::const_val::eps;
    bool sm_flg=true;
    for(zsw::Scalar sm=0; sm_flg; sm+=step_m) {
      if(sm>max_sm) { sm=max_sm; sm_flg=false; }
      samples.push_back(v0+tmp_n+sm*m);
    }
  }
}
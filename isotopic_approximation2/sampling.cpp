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
    if(sn>1-zsw::const_val::eps) { sn=1; sn_flg=false; }
    const Eigen::Matrix<zsw::Scalar,3,1> tmp_n=sn*n;
    const zsw::Scalar max_sm=1-sn-zsw::const_val::eps;
    bool sm_flg=true;
    for(zsw::Scalar sm=0; sm_flg; sm+=step_m) {
      if(sm>max_sm) { sm=max_sm; sm_flg=false; }
      samples.push_back(v0+tmp_n+sm*m);
    }
  }
}

void zsw::sampleTriangleRefTriangle(const Eigen::Matrix<zsw::Scalar, 3, 3> &sample_tri_pts,
                                    const Eigen::Matrix<zsw::Scalar, 3, 3> &ref_tri_pts, const zsw::Scalar r,
                                    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &samples)
{
  Eigen::Matrix<zsw::Scalar,3,1> rv0=ref_tri_pts.block<3,1>(0,0);
  Eigen::Matrix<zsw::Scalar,3,1> rv1=ref_tri_pts.block<3,1>(0,1);
  Eigen::Matrix<zsw::Scalar,3,1> rv2=ref_tri_pts.block<3,1>(0,2);
  Eigen::Matrix<zsw::Scalar,3,1> rn=rv1-rv0;
  Eigen::Matrix<zsw::Scalar,3,1> rm=rv2-rv0;

  zsw::Scalar step_n=r/rn.norm();
  zsw::Scalar step_m=r/rm.norm();

  Eigen::Matrix<zsw::Scalar,3,1> v0=sample_tri_pts.block<3,1>(0,0);
  Eigen::Matrix<zsw::Scalar,3,1> v1=sample_tri_pts.block<3,1>(0,1);
  Eigen::Matrix<zsw::Scalar,3,1> v2=sample_tri_pts.block<3,1>(0,2);
  Eigen::Matrix<zsw::Scalar,3,1> n=v1-v0;
  Eigen::Matrix<zsw::Scalar,3,1> m=v2-v0;

  bool sn_flg=true;
  for(zsw::Scalar sn=0; sn_flg; sn+=step_n) {
    if(sn>1-zsw::const_val::eps) { sn=1; sn_flg=false; }
    const Eigen::Matrix<zsw::Scalar,3,1> tmp_n=sn*n;
    const zsw::Scalar max_sm=1-sn-zsw::const_val::eps;
    bool sm_flg=true;
    for(zsw::Scalar sm=0; sm_flg; sm+=step_m) {
      if(sm>max_sm) { sm=max_sm; sm_flg=false; }
      samples.push_back(v0+tmp_n+sm*m);
    }
  }
}

void zsw::sampleTriangleAndDeformedTriangle(const Eigen::Matrix<zsw::Scalar,3,3> &deformed_tri_pts,
                                            const Eigen::Matrix<zsw::Scalar,3,3> &ori_tri_pts, const zsw::Scalar ori_r,
                                            std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_samples,
                                            std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &ori_samples)
{
  Eigen::Matrix<zsw::Scalar,3,1> ori_v0 = ori_tri_pts.block<3,1>(0,0);
  Eigen::Matrix<zsw::Scalar,3,1> ori_v1 = ori_tri_pts.block<3,1>(0,1);
  Eigen::Matrix<zsw::Scalar,3,1> ori_v2 = ori_tri_pts.block<3,1>(0,2);
  Eigen::Matrix<zsw::Scalar,3,1> ori_n = ori_v1-ori_v0;
  Eigen::Matrix<zsw::Scalar,3,1> ori_m = ori_v2-ori_v0;

  zsw::Scalar step_n=ori_r/ori_n.norm();
  zsw::Scalar step_m=ori_r/ori_m.norm();

  Eigen::Matrix<zsw::Scalar,3,1> deformed_v0 = deformed_tri_pts.block<3,1>(0,0);
  Eigen::Matrix<zsw::Scalar,3,1> deformed_v1 = deformed_tri_pts.block<3,1>(0,1);
  Eigen::Matrix<zsw::Scalar,3,1> deformed_v2 = deformed_tri_pts.block<3,1>(0,2);
  Eigen::Matrix<zsw::Scalar,3,1> deformed_n = deformed_v1 - deformed_v0;
  Eigen::Matrix<zsw::Scalar,3,1> deformed_m = deformed_v2 - deformed_v0;

  bool sn_flg=true;
  for(zsw::Scalar sn=0; sn_flg; sn+=step_n) {
    if(sn>1-zsw::const_val::eps) { sn=1; sn_flg=false; }
    const Eigen::Matrix<zsw::Scalar,3,1> ori_tmp_n=sn*ori_n;
    const Eigen::Matrix<zsw::Scalar,3,1> deformed_tmp_n=sn*deformed_n;
    const zsw::Scalar max_sm=1-sn-zsw::const_val::eps;
    bool sm_flg=true;
    for(zsw::Scalar sm=0; sm_flg; sm+=step_m) {
      if(sm>max_sm) { sm=max_sm; sm_flg=false; }
      ori_samples.push_back(ori_v0+ori_tmp_n+sm*ori_m);
      deformed_samples.push_back(deformed_v0+deformed_tmp_n+sm*deformed_m);
    }
  }
}

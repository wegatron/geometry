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

  void zsw::sampleTriangleAndRefTriangle(const Eigen::Matrix<zsw::Scalar,3,3> &sample_tri_pts,
                                    const Eigen::Matrix<zsw::Scalar,3,3> &ref_tri_pts, const zsw::Scalar ref_r,
                                    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &samples,
                                    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &ref_samples)
  {
    Eigen::Matrix<zsw::Scalar,3,1> rv0=ref_tri_pts.block<3,1>(0,0);
    Eigen::Matrix<zsw::Scalar,3,1> rv1=ref_tri_pts.block<3,1>(0,1);
    Eigen::Matrix<zsw::Scalar,3,1> rv2=ref_tri_pts.block<3,1>(0,2);
    Eigen::Matrix<zsw::Scalar,3,1> rn=rv1-rv0;
    Eigen::Matrix<zsw::Scalar,3,1> rm=rv2-rv0;

    zsw::Scalar step_n=ref_r/rn.norm();
    zsw::Scalar step_m=ref_r/rm.norm();

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
        ref_samples.push_back(rv0+tmp_n+sm*rm);
      }
    }
  }

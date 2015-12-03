#include "constraint.h"

#include <iostream>
#include <zswlib/const_val.h>

void zsw::KernelRegionJudger::addConstraint(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> &v1,
                       const Eigen::Matrix<zsw::Scalar,3,1> &v2, const Eigen::Matrix<zsw::Scalar,3,1> &vr)
{
    vec_v0.push_back(v0);
    Eigen::Matrix<zsw::Scalar,3,1> va=v1-v0;
    Eigen::Matrix<zsw::Scalar,3,1> vb=v2-v0;
    Eigen::Matrix<zsw::Scalar,3,1> vn=va.cross(vb);

    if(vn.norm()<zsw::const_val::eps) {
      std::cerr << "nv norm too small:" << vn.norm();
    }
    //assert(vn.norm()>zsw::const_val::eps)
    vn.normalize();
    if(vn.dot(vr-v0) < 0) { vn=-vn; }
    vec_vn.push_back(vn);
}

bool zsw::KernelRegionJudger::judge(const Eigen::Matrix<zsw::Scalar,3,1> &pt)
{
#if 1
  for(size_t i=0; i<vec_v0.size(); ++i) {
    if(vec_vn[i].dot(pt-vec_v0[i]) < 10*zsw::const_val::eps) {
      return false;
    }
  }
#else
  std::cerr << "only judge for condition " << condition_i_ << std::endl;
  std::cerr << "vec_vn:" <<  vec_vn[condition_i_].transpose() << std::endl;
  std::cerr << "vec_v0:" << vec_v0[condition_i_].transpose() << std::endl;
  std::cerr << "pt:" << pt.transpose() << std::endl;
  if(condition_i_!=-1 && vec_vn[condition_i_].dot(pt-vec_v0[condition_i_]) < zsw::const_val::eps) {
    return false;
  }
#endif
  return true;
}

void zsw::NormalConditionJudger::addConstraint(const size_t vid,
                                               const Eigen::Matrix<zsw::Scalar,3,1> &bev0,
                                               const Eigen::Matrix<zsw::Scalar,3,1> &bev1,
                                               const Eigen::Matrix<zsw::Scalar,3,1> &normal)
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
}

bool zsw::NormalConditionJudger::judge(const Eigen::Matrix<zsw::Scalar,3,1> &pt)
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  return false;
}

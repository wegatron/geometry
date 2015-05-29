#include "integrate.h"

using namespace zsw;

Eigen::Vector3d zsw::AdVectorIntegrator::operator()(const double* pos) const
{
  if(vfs_.size() < 3) {
    Eigen::Vector3d vel;
    vfs_[vfs_.size()-1]->val(pos, vel.data());
    return vel*h_;
  }

  Eigen::Vector3d ori_pos, tmp_pos;
  std::copy(pos, pos+3, ori_pos.data());
  Eigen::Vector3d k1;
  vfs_[vfs_.size()-3]->val(pos, k1.data());

  Eigen::Vector3d k2;
  tmp_pos = ori_pos + h_*k1;
  vfs_[vfs_.size()-2]->val(tmp_pos.data(), k2.data());

  Eigen::Vector3d k3;
  tmp_pos = ori_pos +2*h_*(-k1+2*k2);
  vfs_[vfs_.size()-1]->val(tmp_pos.data(), k3.data());
  return (k1+4*k2+k3)/6*h_;
}

void zsw::AdVectorIntegrator::pushVectorField(std::shared_ptr<VectorField> vf)
{
  #if 0
  if(vfs_.size()>=3) {
    vfs_.erase(vfs_.begin()+vfs_.size()-1);
  }
  # endif
  vfs_.push_back(vf);
}

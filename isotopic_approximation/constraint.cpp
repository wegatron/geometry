#include "constraint.h"

#include <iostream>
#include <fstream>
#include <zswlib/const_val.h>
#include <zswlib/mesh/vtk.h>

void zsw::KernelRegionJudger::addConstraint(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> &v1,
                                            const Eigen::Matrix<zsw::Scalar,3,1> &v2, const Eigen::Matrix<zsw::Scalar,3,1> &vr)
{
  if(!isgood_) { return; }
  vec_v0.push_back(v0);
  Eigen::Matrix<zsw::Scalar,3,1> va=v1-v0;
  Eigen::Matrix<zsw::Scalar,3,1> vb=v2-v0;
  Eigen::Matrix<zsw::Scalar,3,1> vn=va.cross(vb);

  if(vn.norm()<zsw::const_val::eps) {
    std::cerr << "nv norm too small:" << vn.norm() << std::endl;;
    isgood_=false;
    return;
  }
  //assert(vn.norm()>zsw::const_val::eps)
  vn.normalize();
  if(vn.dot(vr-v0) < 0) { vn=-vn; }
  vec_vn.push_back(vn);

  #if 1
  if(debug_) {
    std::cerr << "add krj:" << v0.transpose() << "#" << v1.transpose() << "#" << v2.transpose() << "#" << vr.transpose() << std::endl;
    std::cerr << "normal:" << vn.transpose() << std::endl;
  }
  #endif
}

bool zsw::KernelRegionJudger::judge(const Eigen::Matrix<zsw::Scalar,3,1> &pt)
{
#if 1
  if(!isgood_) { return false; }
  for(size_t i=0; i<vec_v0.size(); ++i) {
    if(vec_vn[i].dot(pt-vec_v0[i]) < precision_) {
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

void zsw::NormalConditionJudger::addConstraint(const Eigen::Matrix<zsw::Scalar,3,1> &bev0,
                                               const Eigen::Matrix<zsw::Scalar,3,1> &bev1,
                                               const Eigen::Matrix<zsw::Scalar,3,1> &normal)
{
  bev_[0].push_back(bev0);
  bev_[1].push_back(bev1);
  normals_.push_back(normal);
}

bool zsw::NormalConditionJudger::judge(const Eigen::Matrix<zsw::Scalar,3,1> &pt)
{
  for(size_t i=0; i<bev_[0].size(); ++i) {
    Eigen::Matrix<zsw::Scalar,3,1> cur_normal = (bev_[0][i]-pt).cross(bev_[1][i]-pt);
    cur_normal.normalize();
    if(cur_normal.dot(normals_[i]) < tol_) { return false; }
  }
  return true;
}

bool zsw::normalCondition(
                          const std::vector<zsw::Vertex> &vertices,
                          const std::vector<Eigen::Matrix<size_t,3,1>> &bound_tris,
                          const Eigen::Matrix<zsw::Scalar,3,1> &pt,
                          const zsw::Scalar pt_val,
                          const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bi_jpts,
                          const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bo_jpts,
                          std::shared_ptr<zsw::Flann<zsw::Scalar>> jpts_ptr_bi,
                          std::shared_ptr<zsw::Flann<zsw::Scalar>> jpts_ptr_bo)
{
  // for each tet
  for(const Eigen::Matrix<size_t,3,1> &b_tr : bound_tris) {
    Eigen::Matrix<zsw::Scalar,3,1> nv;
    for(size_t i=0; i<3; ++i) {
      if(vertices[b_tr[i]].pt_type_==zsw::INNER_POINT) { nv[i]=-1.0-pt_val; }
      else if(vertices[b_tr[i]].pt_type_==zsw::ZERO_POINT) { nv[i]=0.0-pt_val; }
      else { nv[i]=1.0-pt_val; }
    }

    Eigen::Matrix<zsw::Scalar,3,Eigen::Dynamic> tet_pts(3,4);
    tet_pts.block<3,1>(0,0)=pt; tet_pts.block<3,1>(0,1)=vertices[b_tr[0]].pt_;
    tet_pts.block<3,1>(0,2)=vertices[b_tr[1]].pt_; tet_pts.block<3,1>(0,3)=vertices[b_tr[2]].pt_;
    Eigen::Matrix<zsw::Scalar,3,Eigen::Dynamic> ori_tet_pts=tet_pts;
    const Eigen::Matrix<zsw::Scalar,3,1> v0=pt;
    Eigen::Matrix<zsw::Scalar,3,3> A;
    A.block<3,1>(0,0)=tet_pts.block<3,1>(0,1)-v0;
    A.block<3,1>(0,1)=tet_pts.block<3,1>(0,2)-v0;
    A.block<3,1>(0,2)=tet_pts.block<3,1>(0,3)-v0;
    Eigen::PartialPivLU<Eigen::Matrix<zsw::Scalar,3,3>> pplu; pplu.compute(A);

    // scale 70%
    Eigen::Matrix<zsw::Scalar,3,1> center=tet_pts.block<3,1>(0,0);
    center+=tet_pts.block<3,1>(0,1); center+=tet_pts.block<3,1>(0,2);
    center+=tet_pts.block<3,1>(0,3); center*=0.25;
    tet_pts*=0.3;
    tet_pts+= 0.7*center*Eigen::Matrix<zsw::Scalar,1,4>::Ones();

    {
      std::vector<size_t> bi_indices;
      std::vector<zsw::Scalar> bi_dist;
      jpts_ptr_bi->queryNearest(tet_pts, bi_indices, bi_dist);
      for(size_t ind : bi_indices) {
        Eigen::Matrix<zsw::Scalar,3,1> ans = pplu.solve(bi_jpts[ind]-v0);
        if((A*ans-(bi_jpts[ind]-v0)).norm()>zsw::const_val::eps) { std::cout << __FILE__ << __LINE__ << std::endl; abort(); }
        if(pt_val+ans.dot(nv)>0.5) {
          return false;
        }
      }
    }

    {
      std::vector<size_t> bo_indices;
      std::vector<zsw::Scalar> bo_dist;
      jpts_ptr_bi->queryNearest(tet_pts, bo_indices, bo_dist);
      for(size_t ind : bo_indices) {
        Eigen::Matrix<zsw::Scalar,3,1> ans = pplu.solve(bo_jpts[ind]-v0);
        if((A*ans-(bo_jpts[ind]-v0)).norm()>zsw::const_val::eps) { std::cout << __FILE__ << __LINE__ << std::endl; abort(); }
        if(pt_val+ans.dot(nv)<-0.5) { return false; }
      }
    }
  }
  return true;
}

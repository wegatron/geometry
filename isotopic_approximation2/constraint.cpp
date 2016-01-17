#include "constraint.h"

#include <iostream>
#include <fstream>
#include <zswlib/const_val.h>
#include <zswlib/mesh/vtk.h>

#include "basic_data_structure.h"

void zsw::KernelRegionJudger::addConstraint(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> &v1,
                                            const Eigen::Matrix<zsw::Scalar,3,1> &v2, const Eigen::Matrix<zsw::Scalar,3,1> &vr)
{
  if(!isgood_) { return; }
  vec_v0_.push_back(v0);
  Eigen::Matrix<zsw::Scalar,3,1> va=v1-v0;
  Eigen::Matrix<zsw::Scalar,3,1> vb=v2-v0;
  Eigen::Matrix<zsw::Scalar,3,1> vn=va.cross(vb);

  if(vn.norm()<10*zsw::const_val::eps) {
    isgood_=false;
    return;
  }
  //assert(vn.norm()>zsw::const_val::eps)
  vn.normalize();
  if(vn.dot(vr-v0) < 0) { vn=-vn; }
  vec_vn_.push_back(vn);
}

bool zsw::KernelRegionJudger::judge(const Eigen::Matrix<zsw::Scalar,3,1> &pt)
{
  if(!isgood_) { return false; }
  for(size_t i=0; i<vec_v0_.size(); ++i) {
    if(vec_vn_[i].dot(pt-vec_v0_[i]) < precision_) {      return false;    }
  }
  return true;
}

bool zsw::normalCondition(
                          const Eigen::Matrix<zsw::Scalar,4,1> &val,
                          const Eigen::Matrix<zsw::Scalar,3,4> &scaled_tri_pts,
                          const Eigen::Matrix<zsw::Scalar,3,4> &tri_pts,
                          const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &inner_jpts,
                          const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &outer_jpts,
                          std::shared_ptr<zsw::Flann<zsw::Scalar>> inner_kdtree_ptr,
                          std::shared_ptr<zsw::Flann<zsw::Scalar>> outer_kdtree_ptr,
                          bool debug_flag,
                          const std::string *filepath_ptr)
{
  // find the four inner_jpts and four outer_jpts
  std::vector<size_t> in_indices;
  std::vector<zsw::Scalar> in_dist;
  inner_kdtree_ptr->queryNearest(scaled_tri_pts, in_indices, in_dist);

  assert(fabs((inner_jpts[in_indices[0]]-scaled_tri_pts.block<3,1>(0,0)).squaredNorm()-in_dist[0])<zsw::const_val::eps);
  assert(fabs((inner_jpts[in_indices[1]]-scaled_tri_pts.block<3,1>(0,1)).squaredNorm()-in_dist[1])<zsw::const_val::eps);
  assert(fabs((inner_jpts[in_indices[2]]-scaled_tri_pts.block<3,1>(0,2)).squaredNorm()-in_dist[2])<zsw::const_val::eps);
  assert(fabs((inner_jpts[in_indices[3]]-scaled_tri_pts.block<3,1>(0,3)).squaredNorm()-in_dist[3])<zsw::const_val::eps);

  std::vector<size_t> out_indices;
  std::vector<zsw::Scalar> out_dist;
  outer_kdtree_ptr->queryNearest(scaled_tri_pts, out_indices, out_dist);

  assert(fabs((outer_jpts[out_indices[0]]-scaled_tri_pts.block<3,1>(0,0)).squaredNorm()-out_dist[0])<zsw::const_val::eps);
  assert(fabs((outer_jpts[out_indices[1]]-scaled_tri_pts.block<3,1>(0,1)).squaredNorm()-out_dist[1])<zsw::const_val::eps);
  assert(fabs((outer_jpts[out_indices[2]]-scaled_tri_pts.block<3,1>(0,2)).squaredNorm()-out_dist[2])<zsw::const_val::eps);
  assert(fabs((outer_jpts[out_indices[3]]-scaled_tri_pts.block<3,1>(0,3)).squaredNorm()-out_dist[3])<zsw::const_val::eps);

  if(debug_flag){
    assert(filepath_ptr!=nullptr);
    std::vector<Eigen::Matrix<zsw::Scalar,3,4>> cell_out;
    cell_out.push_back(tri_pts);
    cell_out.push_back(scaled_tri_pts);
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> pts;
    pts.push_back(inner_jpts[in_indices[0]]);
    pts.push_back(inner_jpts[in_indices[1]]);
    pts.push_back(inner_jpts[in_indices[2]]);
    pts.push_back(inner_jpts[in_indices[3]]);
    pts.push_back(outer_jpts[out_indices[0]]);
    pts.push_back(outer_jpts[out_indices[1]]);
    pts.push_back(outer_jpts[out_indices[2]]);
    pts.push_back(outer_jpts[out_indices[3]]);
    writeCellsAndPoints(*filepath_ptr, cell_out, pts);
  }

  // judge if can classify the eight point
  Eigen::Matrix<zsw::Scalar,4,4> A;
  A.block<3,4>(0,0)=tri_pts;
  A(3,0)=A(3,1)=A(3,2)=A(3,3)=1;
  Eigen::PartialPivLU<Eigen::Matrix<zsw::Scalar,4,4>> pplu; pplu.compute(A);
  for(size_t in=0;in<4; ++in) {
    Eigen::Matrix<zsw::Scalar,4,1> b;
    b.block<3,1>(0,0)=inner_jpts[in_indices[in]]; b(3,0)=1;
    Eigen::Matrix<zsw::Scalar,4,1> x=pplu.solve(b);
    if(val.dot(x)>0) {
        if(debug_flag) {
            std::cerr << "dot:" << val.dot(x) << std::endl;
            std::cerr << "inner jpt:" << inner_jpts[in_indices[in]].transpose() << std::endl;
            std::cerr << "val:" << val.transpose() << std::endl;
            std::cerr << "x:" << x.transpose() << std::endl;
            std::cerr << "A" << A << std::endl;
          }
      return false;
    }
  }
  for(size_t out=0; out<4; ++out) {
    Eigen::Matrix<zsw::Scalar,4,1> b;
    b.block<3,1>(0,0)=outer_jpts[out_indices[out]]; b(3,0)=1;
    Eigen::Matrix<zsw::Scalar,4,1> x=pplu.solve(b);
    if(val.dot(x)<0) {
        if(debug_flag) {
            std::cerr << "dot:" << val.dot(x) << std::endl;
            std::cerr << "outer jpt:" << outer_jpts[out_indices[out]].transpose() << std::endl;
            std::cerr << "val:" << val.transpose() << std::endl;
            std::cerr << "x:" << x.transpose() << std::endl;
            std::cerr << "A" << A << std::endl;
          }
      return false;
    }
  }
  return true;
}

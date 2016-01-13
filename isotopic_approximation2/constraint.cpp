#include "constraint.h"

#include <iostream>
#include <fstream>
#include <zswlib/const_val.h>
#include <zswlib/mesh/vtk.h>

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
// bool zsw::normalCondition(
//                           const std::vector<zsw::Vertex> &vertices,
//                           const std::vector<Eigen::Matrix<size_t,3,1>> &bound_tris,
//                           const Eigen::Matrix<zsw::Scalar,3,1> &pt,
//                           const zsw::Scalar pt_val,
//                           const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bi_jpts,
//                           const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bo_jpts,
//                           std::shared_ptr<zsw::Flann<zsw::Scalar>> jpts_ptr_bi,
//                           std::shared_ptr<zsw::Flann<zsw::Scalar>> jpts_ptr_bo)
// {
//   // for each tet
//   for(const Eigen::Matrix<size_t,3,1> &b_tr : bound_tris) {
//     if(pt_val<0 && vertices[b_tr[0]].pt_type_!=OUTER_POINT && vertices[b_tr[1]].pt_type_!=OUTER_POINT
//        && vertices[b_tr[2]].pt_type_!=OUTER_POINT) { continue; }
//     if(pt_val>0 && vertices[b_tr[0]].pt_type_!=INNER_POINT && vertices[b_tr[1]].pt_type_!=INNER_POINT
//        && vertices[b_tr[2]].pt_type_!=INNER_POINT) { continue; }
//     Eigen::Matrix<zsw::Scalar,3,1> nv;
//     for(size_t i=0; i<3; ++i) {
//       if(vertices[b_tr[i]].pt_type_==zsw::INNER_POINT) { nv[i]=-1.0-pt_val; }
//       else if(vertices[b_tr[i]].pt_type_==zsw::ZERO_POINT) { nv[i]=0.0-pt_val; }
//       else { nv[i]=1.0-pt_val; }
//     }

//     Eigen::Matrix<zsw::Scalar,3,Eigen::Dynamic> tet_pts(3,4);
//     tet_pts.block<3,1>(0,0)=pt; tet_pts.block<3,1>(0,1)=vertices[b_tr[0]].pt_;
//     tet_pts.block<3,1>(0,2)=vertices[b_tr[1]].pt_; tet_pts.block<3,1>(0,3)=vertices[b_tr[2]].pt_;
// #if 0
//     const Eigen::Matrix<zsw::Scalar,3,Eigen::Dynamic> ori_tet_pts=tet_pts;
// #endif

//     const Eigen::Matrix<zsw::Scalar,3,1> v0=pt;
//     Eigen::Matrix<zsw::Scalar,3,3> A;
//     A.block<3,1>(0,0)=tet_pts.block<3,1>(0,1)-v0;
//     A.block<3,1>(0,1)=tet_pts.block<3,1>(0,2)-v0;
//     A.block<3,1>(0,2)=tet_pts.block<3,1>(0,3)-v0;
//     Eigen::PartialPivLU<Eigen::Matrix<zsw::Scalar,3,3>> pplu; pplu.compute(A);

//     // scale 70%
//     Eigen::Matrix<zsw::Scalar,3,1> center=tet_pts.block<3,1>(0,0);
//     center+=tet_pts.block<3,1>(0,1); center+=tet_pts.block<3,1>(0,2);
//     center+=tet_pts.block<3,1>(0,3); center*=0.25;
//     tet_pts*=0.3;
//     tet_pts+= 0.7*center*Eigen::Matrix<zsw::Scalar,1,4>::Ones();

//     {
//       std::vector<size_t> bi_indices;
//       std::vector<zsw::Scalar> bi_dist;
//       jpts_ptr_bi->queryNearest(tet_pts, bi_indices, bi_dist);
//       for(size_t ind : bi_indices) {
//         Eigen::Matrix<zsw::Scalar,3,1> ans = pplu.solve(bi_jpts[ind]-v0);
//         assert((A*ans-(bi_jpts[ind]-v0)).norm()<zsw::const_val::eps);
//         if(pt_val+ans.dot(nv)>0) {
// #if 0
//           {
//             static size_t npo=0;
//             ++npo;
//             {
//               size_t tet_id[4]={0,1,2,3};
//               std::ofstream ofs0("/home/wegatron/tmp/debug_no_scale_tet_"+std::to_string(npo)+".vtk");
//               tet2vtk(ofs0, ori_tet_pts.data(), 4, tet_id, 1);
//               std::ofstream ofs1("/home/wegatron/tmp/debug_scaled_tet_"+std::to_string(npo)+".vtk");
//               tet2vtk(ofs1, tet_pts.data(), 4, tet_id, 1);
//             }
//             {
//               Eigen::Matrix<zsw::Scalar,3,4> nti;
//               nti.block<3,1>(0,0)=bi_jpts[bi_indices[0]];
//               nti.block<3,1>(0,1)=bi_jpts[bi_indices[1]];
//               nti.block<3,1>(0,2)=bi_jpts[bi_indices[2]];
//               nti.block<3,1>(0,3)=bi_jpts[bi_indices[3]];
//               size_t tet_id[4]={0,1,2,3};
//               std::ofstream ofs("/home/wegatron/tmp/debug_nti_"+std::to_string(npo)+".vtk");
//               tet2vtk(ofs, nti.data(), 4, tet_id, 1);
//             }
//             std::cerr << "pt:" << bi_jpts[ind].transpose() << std::endl;
//             std::cerr << "ans:" << ans.transpose() << std::endl;
//             std::cerr << "nv:" << nv.transpose() << std::endl;
//             std::cerr << "A:" << A << std::endl;
//             std::cerr << "nti ret!!!" << std::endl;
//           }
// #endif
//           return false;
//         }
//       }
//     }

//     {
//       std::vector<size_t> bo_indices;
//       std::vector<zsw::Scalar> bo_dist;
//       jpts_ptr_bo->queryNearest(tet_pts, bo_indices, bo_dist);
//       for(size_t ind : bo_indices) {
//         Eigen::Matrix<zsw::Scalar,3,1> ans = pplu.solve(bo_jpts[ind]-v0);
//         assert((A*ans-(bo_jpts[ind]-v0)).norm()<zsw::const_val::eps);
//         if(pt_val+ans.dot(nv)<0) {
// #if 0
//           {
//             static size_t npo=0;
//             ++npo;
//             {
//               size_t tet_id[4]={0,1,2,3};
//               std::ofstream ofs0("/home/wegatron/tmp/debug_no_scale_tet_"+std::to_string(npo)+".vtk");
//               tet2vtk(ofs0, ori_tet_pts.data(), 4, tet_id, 1);
//               std::ofstream ofs1("/home/wegatron/tmp/debug_scaled_tet_"+std::to_string(npo)+".vtk");
//               tet2vtk(ofs1, tet_pts.data(), 4, tet_id, 1);
//             }
//             {
//               Eigen::Matrix<zsw::Scalar,3,4> nto;
//               nto.block<3,1>(0,0)=bo_jpts[bo_indices[0]];
//               nto.block<3,1>(0,1)=bo_jpts[bo_indices[1]];
//               nto.block<3,1>(0,2)=bo_jpts[bo_indices[2]];
//               nto.block<3,1>(0,3)=bo_jpts[bo_indices[3]];
//               size_t tet_id[4]={0,1,2,3};
//               std::ofstream ofs("/home/wegatron/tmp/debug_nto_"+std::to_string(npo)+".vtk");
//               tet2vtk(ofs, nto.data(), 4, tet_id, 1);
//             }
//             std::cerr << "nto ret!!!" << std::endl;
//           }
// #endif
//           return false;
//         }
//       }
//     }
//   }
//   return true;
// }
  return true;
}

bool zsw::normalCondition(
                          const Eigen::Matrix<zsw::Scalar,4,1> &val,
                          const Eigen::Matrix<zsw::Scalar,3,4> &tri_pts,
                          const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &inner_jpts,
                          const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &outer_jpts,
                          std::shared_ptr<zsw::Flann<zsw::Scalar>> inner_kdtree_ptr,
                          std::shared_ptr<zsw::Flann<zsw::Scalar>> outer_kdtree_ptr)
{
  // find the four inner_jpts and four outer_jpts
  std::vector<size_t> in_indices;
  std::vector<zsw::Scalar> in_dist;
  inner_kdtree_ptr->queryNearest(tri_pts, in_indices, in_dist);

  assert(fabs((inner_jpts[in_indices[0]]-tri_pts.block<3,1>(0,0)).squaredNorm()-in_dist[0])<zsw::const_val::eps);
  assert(fabs((inner_jpts[in_indices[1]]-tri_pts.block<3,1>(0,1)).squaredNorm()-in_dist[1])<zsw::const_val::eps);
  assert(fabs((inner_jpts[in_indices[2]]-tri_pts.block<3,1>(0,2)).squaredNorm()-in_dist[2])<zsw::const_val::eps);
  assert(fabs((inner_jpts[in_indices[3]]-tri_pts.block<3,1>(0,3)).squaredNorm()-in_dist[3])<zsw::const_val::eps);

  std::vector<size_t> out_indices;
  std::vector<zsw::Scalar> out_dist;
  outer_kdtree_ptr->queryNearest(tri_pts, out_indices, out_dist);

  assert(fabs((outer_jpts[out_indices[0]]-tri_pts.block<3,1>(0,0)).squaredNorm()-out_dist[0])<zsw::const_val::eps);
  assert(fabs((outer_jpts[out_indices[1]]-tri_pts.block<3,1>(0,1)).squaredNorm()-out_dist[1])<zsw::const_val::eps);
  assert(fabs((outer_jpts[out_indices[2]]-tri_pts.block<3,1>(0,2)).squaredNorm()-out_dist[2])<zsw::const_val::eps);
  assert(fabs((outer_jpts[out_indices[3]]-tri_pts.block<3,1>(0,3)).squaredNorm()-out_dist[3])<zsw::const_val::eps);

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
      // std::cerr << "inner jpts:" << in_indices[in] << std::endl;
      // std::cerr << "pt:" << inner_jpts[in_indices[in]].transpose() << std::endl;
      return false;
    }
  }
  for(size_t out=0; out<4; ++out) {
    Eigen::Matrix<zsw::Scalar,4,1> b;
    b.block<3,1>(0,0)=outer_jpts[out_indices[out]]; b(3,0)=1;
    Eigen::Matrix<zsw::Scalar,4,1> x=pplu.solve(b);
    if(val.dot(x)<0) {
      // std::cerr << "outer jpts:" << out_indices[out] << std::endl;
      // std::cerr << "pt:" << inner_jpts[out_indices[out]].transpose() << std::endl;
      return false;
    }
  }
  return true;
}

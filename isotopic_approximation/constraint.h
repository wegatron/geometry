#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <zswlib/config.h>
#include <zswlib/const_val.h>
#include <zswlib/mesh/zsw_flann2.h>
#include "basic_data_structure.h"

namespace zsw
{

  /// \brief kernel region judger.
  ///
  ///  judge wether the point is in kernel region.
  ///
  class KernelRegionJudger final
  {
  public:
    KernelRegionJudger(const zsw::Scalar precision=10*zsw::const_val::eps) { precision_=precision; isgood_=true; }
    void addConstraint(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> &v1,
                       const Eigen::Matrix<zsw::Scalar,3,1> &v2, const Eigen::Matrix<zsw::Scalar,3,1> &vr);
    bool judge(const Eigen::Matrix<zsw::Scalar,3,1> &pt);
  private:
    zsw::Scalar precision_;
    bool isgood_;
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> vec_v0;
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> vec_vn;
  };


  /// \brief Normal Condition Judger
  ///
  /// judge the boundry or zero surface's new faces' normals is in a difference bound.
  ///
  class NormalConditionJudger final
  {
  public:
    NormalConditionJudger(const zsw::Scalar tol) {
      tol_=tol;
    }

    /// \brief add normal constraint
    ///
    /// A detailed description, it should be 2 lines at least.
    ///
    ///
    /// \param be, boundary edge
    /// \param normal, expected normal merge point with that boundary edge
    /// \param vid, 0 or 1, the collapse edge's vertex id
    void addConstraint(const Eigen::Matrix<zsw::Scalar,3,1> &bev0,
                       const Eigen::Matrix<zsw::Scalar,3,1> &bev1,
                       const Eigen::Matrix<zsw::Scalar,3,1> &normal);

    bool judge(const Eigen::Matrix<zsw::Scalar,3,1> &pt);
  private:
    zsw::Scalar tol_;
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bev_[2];
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> normals_;
  };

  bool normalCondition(
                       const std::vector<zsw::Vertex> &vertices,
                       const std::vector<Eigen::Matrix<size_t,3,1>> &bound_tris,
                       const Eigen::Matrix<zsw::Scalar,3,1> &pt,
                       const zsw::Scalar pt_vals,
                       const std::vector<zsw::JudgePoint> &bi_jpts,
                       const std::vector<zsw::JudgePoint> &bo_jpts,
                       std::shared_ptr<zsw::Flann2<zsw::Scalar,2>> jpts_ptr_bi,
                       std::shared_ptr<zsw::Flann2<zsw::Scalar,2>> jpts_ptr_bo);

  /* bool normalCond(const Eigen::Matrix<zsw::Scalar,3,4> &tet_pts, const Eigen::Matrix<zsw::Scalar,4,1> &pt_vals, */
  /*                 const std::vector<zsw::JudgePoint> &bi_jpts, */
  /*                 const std::vector<zsw::JudgePoint> &bo_jpts, */
  /*                 std::shared_ptr<zsw::Flann2<zsw::Scalar,2>> jpts_ptr_bi, */
  /*                 std::shared_ptr<zsw::Flann2<zsw::Scalar,2>> jpts_ptr_bo) */
  /* { */
  /*   const Eigen::Matrix<zsw::Scalar,3,1> v0=tet_pts.block<3,1>(0,0); */
  /*   Eigen::Matrix<zsw::Scalar,3,3> A; */
  /*   A.block<3,1>(0,0)=tet_pts.block<3,1>(0,1)-v0; */
  /*   A.block<3,1>(0,1)=tet_pts.block<3,1>(0,2)-v0; */
  /*   A.block<3,1>(0,2)=tet_pts.block<3,1>(0,3)-v0; */
  /*   Eigen::PartialPivLU<Eigen::Matrix<zsw::Scalar,3,3>> pplu; pplu.compute(A); */

  /*   Eigen::Matrix<zsw::Scalar,3,1> nv; */
  /*   for(size_t i=0; i<3; ++i) {      nv[i]=pt_vals[i+1]-pt_vals[0];    } */

  /*   Eigen::Matrix<zsw::Scalar,3,1> center=v0; center+=tet_pts.block<3,1>(0,1); */
  /*   center+=tet_pts.block<3,1>(0,2); center+=tet_pts.block<3,1>(0,3); center*=0.25; */
  /*   Eigen::Matrix<zsw::Scalar,3,4> scaled_tet_pts=tet_pts*0.7+0.3*center*Eigen::Matrix<zsw::Scalar,1,4>::Ones(); */

  /*   // bi points */
  /*   { */
  /*     std::vector<size_t> bi_indices; */
  /*     std::vector<zsw::Scalar> bi_dist; */
  /*     jpts_ptr_bi->queryNearest(scaled_tet_pts, bi_indices, bi_dist); */
  /*     for(size_t ind : bi_indices) { */
  /*       Eigen::Matrix<zsw::Scalar,3,1> ans = pplu.solve(bi_jpts[ind].pt_-v0); */
  /*       if((A*ans-(bi_jpts[ind].pt_-v0)).norm()>zsw::const_val::eps) { std::cout << __FILE__ << __LINE__ << std::endl; abort(); } */
  /*       if(pt_vals[0]+ans.dot(nv)>0) { return false; } */
  /*     } */
  /*   } */

  /*   // bo points */
  /*   { */
  /*   } */
  /* } */
}


#endif /* CONSTRAINT_H */

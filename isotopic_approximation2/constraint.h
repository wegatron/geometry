#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <zswlib/config.h>
#include <zswlib/const_val.h>
#include "basic_data_structure.h"
#include "kdtree_warp.h"

namespace zsw
{

  /// \brief kernel region judger.
  ///
  ///  judge wether the point is in kernel region.
  ///
  class KernelRegionJudger final
  {
  public:
    KernelRegionJudger(const zsw::Scalar precision=10*zsw::const_val::eps) { precision_=precision; isgood_=true;  }
    void addConstraint(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> &v1,
                       const Eigen::Matrix<zsw::Scalar,3,1> &v2, const Eigen::Matrix<zsw::Scalar,3,1> &vr);
    bool judge(const Eigen::Matrix<zsw::Scalar,3,1> &pt);
  private:
    bool isgood_;
    zsw::Scalar precision_;
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> vec_v0_;
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> vec_vn_;
  };

  bool normalCondition(
                       const Eigen::Matrix<zsw::Scalar,4,1> &val,
                       const Eigen::Matrix<zsw::Scalar,3,4> &scaled_tri_pts,
                       const Eigen::Matrix<zsw::Scalar,3,4> &tri_pts,
                       const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &inner_jpts,
                       const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &outer_jpts,
                       const zsw::KdTreeWarper &inner_kdtree,
                       const zsw::KdTreeWarper &outer_kdtree,
                       bool debug_flag=false,
                       const std::string *filepath_ptr=nullptr);
}


#endif /* CONSTRAINT_H */

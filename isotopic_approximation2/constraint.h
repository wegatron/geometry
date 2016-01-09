#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <zswlib/config.h>
#include <zswlib/const_val.h>
#include <zswlib/mesh/zsw_flann.h>
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

  /* bool normalCondition( */
  /*                      const std::vector<zsw::Vertex> &vertices, */
  /*                      const std::vector<Eigen::Matrix<size_t,3,1>> &bound_tris, */
  /*                      const Eigen::Matrix<zsw::Scalar,3,1> &pt, */
  /*                      const zsw::Scalar pt_vals, */
  /*                      const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bi_jpts, */
  /*                      const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bo_jpts, */
  /*                      std::shared_ptr<zsw::Flann<zsw::Scalar>> jpts_ptr_bi, */
  /*                      std::shared_ptr<zsw::Flann<zsw::Scalar>> jpts_ptr_bo); */
}


#endif /* CONSTRAINT_H */

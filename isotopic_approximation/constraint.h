#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <vector>
#include <Eigen/Dense>
#include <zswlib/config.h>
#include <zswlib/const_val.h>

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

}


#endif /* CONSTRAINT_H */

#ifndef MLS_H
#define MLS_H

#include <vector>
#include <Eigen/Dense>
#include <memory>

#include <zswlib/config.h>

namespace zsw
{
  class DisWeightFunc
  {
  public:
    virtual void calcWeight(
                            const Eigen::Matrix<zsw::Scalar,3,1> &vt,
                            const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &vs,
                            std::vector<zsw::Scalar> &weight) = 0;
  };

  class GaussDisWeightFunc : public DisWeightFunc
  {
  public:
    GaussDisWeightFunc() { c_ = -0.1; }
    virtual void calcWeight(
                            const Eigen::Matrix<zsw::Scalar,3,1> &vt,
                            const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &vs,
                            std::vector<zsw::Scalar> &weight);
  private:
    zsw::Scalar c_;
  };

  class DeformFunc
  {
  public:

  DeformFunc() : dis_weight_func_(new GaussDisWeightFunc()) {}

    /// \brief calc DeformedPos of a vertex according nearby vertices and their deformed positions.
    ///
    /// using moving least square
    ///
    /// \param vs known vertices
    /// \param vs's deformed position
    /// \param vt the vertex we want to calculate
    /// \return dvt the deformed vertex position
    ///
    void calcDeformedPos(
                         const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &vs,
                         const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &dvs,
                         const Eigen::Matrix<zsw::Scalar,3,1> &vt,
                         Eigen::Matrix<zsw::Scalar,3,1> &dvt);
  private:
    void updateAB(const Eigen::Matrix<zsw::Scalar,3,1> &vt,
                  const Eigen::Matrix<zsw::Scalar,3,1> &v,
                  const Eigen::Matrix<zsw::Scalar,3,1> &dv,
                  const zsw::Scalar weitght);

    std::shared_ptr<zsw::DisWeightFunc> dis_weight_func_;
    Eigen::Matrix<zsw::Scalar,12,12> A_;
    Eigen::Matrix<zsw::Scalar,12,1> B_;
  };
}

#endif /* MLS_H */

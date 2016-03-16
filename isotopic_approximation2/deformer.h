#ifndef DEFORMER_H
#define DEFORMER_H

#include <vector>
#include <Eigen/Dense>
#include <memory>

#include <zswlib/config.h>
#include <zswlib/mesh/mesh_type.h>
#include <zswlib/matrix_type.h>
#include <zswlib/mesh/vtk.h>
#include "kdtree_warp.h"

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
    GaussDisWeightFunc() { c_ = -0.2; }
    virtual void calcWeight(
                            const Eigen::Matrix<zsw::Scalar,3,1> &vt,
                            const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &vs,
                            std::vector<zsw::Scalar> &weight);
  private:
    zsw::Scalar c_;
  };

  class LocalDeformFunc
  {
  public:
  LocalDeformFunc() : dis_weight_func_(new GaussDisWeightFunc()) {}
    /// \brief calc DeformedPos of a vertex according nearby vertices and their deformed positions.
    ///
    /// using moving least square
    ///
    /// \param vs known vertices
    /// \param vs's deformed position
    /// \param vt the vertex we want to calculate
    /// \return dvt the deformed vertex position
    ///
    virtual void calcDeformedPos(
                                 const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &vs,
                                 const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &dvs,
                                 const Eigen::Matrix<zsw::Scalar,3,1> &vt,
                                 Eigen::Matrix<zsw::Scalar,3,1> &dvt) = 0;

  protected:
    std::shared_ptr<DisWeightFunc> dis_weight_func_;
  };

  class LocalTranslateDeformFunc :  public LocalDeformFunc
  {
  public:
    LocalTranslateDeformFunc() {}
    void calcDeformedPos(
                         const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &vs,
                         const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &dvs,
                         const Eigen::Matrix<zsw::Scalar,3,1> &vt,
                         Eigen::Matrix<zsw::Scalar,3,1> &dvt);
  };

  class LocalVectorFieldDeformFunc
  {
  public:
    virtual void calcDeformedPos(const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &vs,
                                 const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &dvs,
                                 const Eigen::Matrix<zsw::Scalar,3,1> &vt,
                                 Eigen::Matrix<zsw::Scalar,3,1> &dvt);
  };

  class Deformer
  {
  public:
    Deformer(const zsw::mesh::TriMesh &ori_mesh, const zsw::mesh::TriMesh &deformed_mesh, const zsw::Scalar sample_r);
    const std::vector<zsw::Vector3s>& getRefVs() const { return ref_vs_; }
    const std::vector<zsw::Vector3s>& getRefDvs() const { return ref_dvs_; }
    virtual void deformTo(const std::vector<zsw::Vector3s> &vs, std::vector<zsw::Vector3s> &dvs) = 0;
    virtual void deformBack(const std::vector<zsw::Vector3s> &dvs, std::vector<zsw::Vector3s> &vs) = 0;
  protected:
    std::vector<zsw::Vector3s> ref_vs_;
    std::vector<zsw::Vector3s> ref_dvs_;
    zsw::KdTreeWarper vs_kdt_;
    zsw::KdTreeWarper dvs_kdt_;
  };
  
  class LocalDeformer : public Deformer
  {
  public:
    LocalDeformer(const zsw::mesh::TriMesh &ori_mesh, const zsw::mesh::TriMesh &deformed_mesh,
                  const zsw::Scalar sample_r, std::shared_ptr<LocalDeformFunc> df);
    virtual void deformTo(const std::vector<zsw::Vector3s> &vs, std::vector<zsw::Vector3s> &dvs);
    virtual void deformBack(const std::vector<zsw::Vector3s> &dvs, std::vector<zsw::Vector3s> &vs);
  private:
    std::shared_ptr<LocalDeformFunc> df_;
    zsw::Scalar ref_r_;
  };
}

#endif /* DEFORMER_H */

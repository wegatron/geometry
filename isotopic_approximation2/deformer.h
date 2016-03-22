#ifndef DEFORMER_H
#define DEFORMER_H

#include <vector>
#include <Eigen/Dense>
#include <memory>

#include <zswlib/config.h>
#include <zswlib/mesh/mesh_type.h>
#include <zswlib/matrix_type.h>
#include <zswlib/mesh/vtk.h>
#include <zswlib/mesh/zsw_flann.h>

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

  class Deformer
  {
  public:
    Deformer(const zsw::mesh::TriMesh &ori_mesh, const zsw::mesh::TriMesh &deformed_mesh,
             const zsw::Scalar sample_r, const size_t near_count,
             std::shared_ptr<zsw::DisWeightFunc> dis_weight_func);

    const std::vector<zsw::Vector3s>& getRefVs() const { return ref_vs_; }

    const std::vector<zsw::Vector3s>& getRefDvs() const { return ref_dvs_; }

    void deformTo(std::shared_ptr<std::vector<zsw::Vector3s>> sample_out,
                  std::shared_ptr<std::vector<zsw::Vector3s>> sample_in,
                  std::shared_ptr<std::vector<zsw::Vector3s>> sample_out_d,
                  std::shared_ptr<std::vector<zsw::Vector3s>> sample_in_d);

    void deformBack(const std::vector<zsw::Vector3s> &ptsd, std::vector<zsw::Vector3s> &pts_bk) const;

    virtual void deformTo(const std::vector<zsw::Vector3s> &vs, std::vector<zsw::Vector3s> &dvs) = 0;

  protected:
    size_t near_count_;
    std::shared_ptr<zsw::DisWeightFunc> dis_weight_func_;
    std::vector<zsw::Vector3s> ref_vs_normal_;
    std::vector<zsw::Vector3s> ref_dvs_normal_;
    std::vector<zsw::Vector3s> ref_vs_;
    std::vector<zsw::Vector3s> ref_dvs_;

    std::shared_ptr<std::vector<zsw::Vector3s>> sample_out_;
    std::shared_ptr<std::vector<zsw::Vector3s>> sample_in_;
    std::shared_ptr<std::vector<zsw::Vector3s>> sample_out_d_;
    std::shared_ptr<std::vector<zsw::Vector3s>> sample_in_d_;
    std::shared_ptr<zsw::Flann<zsw::Scalar>> vs_ann_;
    std::shared_ptr<zsw::Flann<zsw::Scalar>> out_d_ann_;
    std::shared_ptr<zsw::Flann<zsw::Scalar>> in_d_ann_;
  };

  class LocalVectorFieldDeformer : public Deformer
  {
  public:
    LocalVectorFieldDeformer(const zsw::mesh::TriMesh &ori_mesh, const zsw::mesh::TriMesh &deformed_mesh,
                             const zsw::Scalar sample_r, const size_t near_count,
                             std::shared_ptr<zsw::DisWeightFunc> dis_weight_func);

    void deformTo(std::shared_ptr<std::vector<zsw::Vector3s>> sample_out,
                  std::shared_ptr<std::vector<zsw::Vector3s>> sample_in,
                  std::shared_ptr<std::vector<zsw::Vector3s>> sample_out_d,
                  std::shared_ptr<std::vector<zsw::Vector3s>> sample_in_d)
    {
      Deformer::deformTo(sample_out, sample_in, sample_out_d,
                         sample_in_d);
    }

    /* void deformBack(const std::vector<zsw::Vector3s> &ptsd, const std::vector<std::vector<size_t>> &adjs, */
    /*                 std::vector<zsw::Vector3s> &pts_bk) const; */

    void deformTo(const std::vector<zsw::Vector3s> &vs, std::vector<zsw::Vector3s> &dvs);

    const std::vector<Eigen::Matrix<zsw::Scalar,3,3>>& getJac() const { return jac_; }

    void writeValidRefVs() const;
  private:
    void calcJac();
    void resolveInvalidJacobian(const std::vector<std::vector<size_t>> &indices);
    std::vector<Eigen::Matrix<zsw::Scalar,3,3>> jac_;
    std::vector<bool> valid_;
  };
}

#endif /* DEFORMER_H */

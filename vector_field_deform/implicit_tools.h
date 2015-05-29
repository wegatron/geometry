#ifndef IMPLICIT_TOOLS_H
#define IMPLICIT_TOOLS_H

#include "integrate.h"

namespace zsw{

  void writeVtk(const std::string& file_path, Eigen::Matrix<double, 3, Eigen::Dynamic> &verts,
                const Eigen::Matrix<size_t, 3, Eigen::Dynamic>& tris);

  class VfDeformer final
  {
  public:
    VfDeformer() {}
    void loadModel(const std::string& file_path);
    void saveModel(const std::string& file_path);
    void setVectorFieldIntegrator(std::shared_ptr<VectorFieldIntegrator> vf_integrator) { vf_integrator_ = vf_integrator; }
    std::shared_ptr<VectorFieldIntegrator> getVectorFieldIntegrator() { return vf_integrator_; }
    void pushVectorFieldAndDeform(std::shared_ptr<VectorField> vf);

    Eigen::Matrix<double, 3, Eigen::Dynamic>& getVerts() { return verts_; }
    Eigen::Matrix<size_t, 3, Eigen::Dynamic>& getTris() { return tris_; }
  private:
    std::shared_ptr<VectorFieldIntegrator> vf_integrator_;
    Eigen::Matrix<double, 3, Eigen::Dynamic> verts_; // default column major
    Eigen::Matrix<size_t, 3, Eigen::Dynamic> tris_; // index start from 0
  };

  class ImplicitTool
  {
  public:
    void setDeformer(std::shared_ptr<VfDeformer> deformer);
    void setTimeSlice(size_t time_slice);
  protected:
    virtual void updateVectorFieldAndDeform() = 0;
    std::shared_ptr<VfDeformer> deformer_;
    size_t time_slice_;
    double r_[2];
  };

  class SphereDeformTool final : public ImplicitTool
  {
  public:
    SphereDeformTool(const double *center, const double ri, const double ro) {
      center_[0] = center[0]; center_[1] = center[1]; center_[2] = center[2];
      r_[0] = ri; r_[1] = ro;
      time_slice_  = 100;
    }
    void translateAndDeform(const double *trans_vec);
  protected:
    void updateVectorFieldAndDeform();
  private:
    void calcU(const Eigen::Vector3d &u_dest, Eigen::Vector3d &u0, Eigen::Vector3d &u1);
    double center_[3];
    // translate vector
    const double* trans_vec_;
  };

  class BendDeformTool final : public ImplicitTool
  {
  public:
    BendDeformTool(const double *b, const double *a, const double *center, const double ri, const double ro);
    void rotateAndDeform(const double theta);
  protected:
    void updateVectorFieldAndDeform();
  private:
    double theta_;
    std::shared_ptr<zsw::VectorField> vf_;
    double a_[3], b_[3], center_[3];
    double ri_, ro_;
  };

  class TwistDeformTool final : public ImplicitTool
  {
  public:
    TwistDeformTool(const double *a, const double *center, const double ri, const double ro);
    void twistAndDeform(const double theta);
  protected:
    virtual void updateVectorFieldAndDeform();
  private:
    std::shared_ptr<VectorField> vf_;
    double angle_v_;
  };
}

#endif /* IMPLICIT_TOOLS_H */

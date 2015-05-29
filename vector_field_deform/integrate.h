#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <vector>
#include <Eigen/Dense>

#include "vector_field.h"

namespace zsw {
  class VectorFieldIntegrator
  {
  public:
    virtual Eigen::Vector3d operator()(const double* pos) const = 0;
    virtual void pushVectorField(std::shared_ptr<VectorField> vf) = 0;
    virtual void setStep(const double h) { h_ = h; }
  protected:
    double h_;
  };

  class AdVectorIntegrator final : public VectorFieldIntegrator
  {
  public:
    AdVectorIntegrator() {
      h_=0.01; // the value is from Fernando
    }
    Eigen::Vector3d operator()(const double* pos) const;
    void pushVectorField(std::shared_ptr<VectorField> vf);
  private:
    std::vector<std::shared_ptr<VectorField>> vfs_;
  };
}

#endif /* INTEGRATE_H */

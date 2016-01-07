#ifndef BOUND_SPHERE_H
#define BOUND_SPHERE_H

#include <vector>
#include <Eigen/Dense>
#include <zswlib/config.h>

namespace zsw{
  class BoundSphere
  {
  public:
    BoundSphere(const std::string &filepath, const zsw::Scalar scale, const Eigen::Matrix<zsw::Scalar,3,1> &transform);
    const std::vector<Eigen::Matrix<zsw::Scalar,3,1>>& getVertices() { return vertices_; }
  private:
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> vertices_;
  };
}

#endif /* BOUND_SPHERE_H */

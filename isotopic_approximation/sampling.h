#ifndef SAMPLING_H
#define SAMPLING_H

#include <vector>
#include <zswlib/mesh/mesh_type.h>

namespace zsw
{

  class Sampler  final
  {
  public:
    Sampler() {}
    void sampleSigmaDense(const zsw::mesh::TriMesh &tm, const zsw::Scalar sigma, std::vector<Eigen::Matrix<zsw::Scalar, 3, 1>> &samples);
#ifdef NDEBUG
  private:
#endif
    void sampleTriangle(Eigen::Matrix<zsw::Scalar, 3, 3> tri_points, const zsw::Scalar sigma, std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &samples);
    void calcLocalCoordinate(const Eigen::Matrix<zsw::Scalar, 3, 3> &tri_points, Eigen::Matrix<zsw::Scalar, 3, 1> &translate, Eigen::Matrix<zsw::Scalar, 3, 3> &rotate);

    /***
     * check if sample point is inner the triangle,
     * if not move onto the triangle, make sure the new circle cover the area of the triangle covered by the earlier circle.
     */
    void resolvePoint(const Eigen::Matrix<zsw::Scalar,3,3> &tri_points, Eigen::Matrix<zsw::Scalar,3,1> &sample_point);
  };
}

#endif /* SAMPLING_H */

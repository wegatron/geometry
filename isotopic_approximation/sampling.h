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
     * if the point should move to the triangle point, then discard this point, as each triangle point will be sample at last.
     */
    bool resolvePoint(const Eigen::Matrix<zsw::Scalar,3,3> &tri_points, Eigen::Matrix<zsw::Scalar,3,1> &sample_point);

    /***
     * check if vr and vt is in the same side of line v0, v1
     * all the point are in the same plane.
     */
    bool sameSide(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> &v1, const Eigen::Matrix<zsw::Scalar,3,1> &vr, const Eigen::Matrix<zsw::Scalar,3,1> &vt);

    void projectToLine(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> &v1, Eigen::Matrix<zsw::Scalar,3,1> &sample_point);
  };

}

#endif /* SAMPLING_H */

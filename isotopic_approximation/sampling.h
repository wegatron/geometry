#ifndef SAMPLING_H
#define SAMPLING_H

#include <zswlib/mesh/mesh_type.h>

namespace zsw
{
  class Sampler
  {
  public:
    Eigen::Matrix<zsw::Scalar, 3, Eigen::Dynamic> sampleSigmaDense(const zsw::mesh::TriMesh &tm);
  private:
    void sampleTriangle(const Eigen::Matrix<zsw::Scalar, 3, 3> vertices, std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &samples);
    void calcLocalCoordinate(const Eigen::Matrix<zsw::Scalar, 3, 3> vertices, Eigen::Matrix<zsw::Scalar, 3, 1> & translate, Eigen::Matrix<zsw::Scalar, 3, 3> &rotate);

    /***
     * calculate the point vt on line (vl0, vl1) which has the minest distance to point vc
     **/
    void calcMinPos(const Eigen::Matrix<zsw::Scalar, 2,1> &vl0, const Eigen::Matrix<zsw::Scalar,2,1> &vl1, const Eigen::Matrix<zsw::Scalar,2,1> &vc, Eigen::Matrix<zsw::Scalar,2,1> &vt);

    /***
     * judge if the point vr and vc is in the same side of line (vl0, vl1)
     */
    bool sameSideLine(const Eigen::Matrix<zsw::Scalar, 2,1> &vl0, const Eigen::Matrix<zsw::Scalar,2,1> &vl1, const Eigen::Matrix<zsw::Scalar,2,1> &vr, const Eigen::Matrix<zsw::Scalar,2,1> &vc);
  };
}

#endif /* SAMPLING_H */

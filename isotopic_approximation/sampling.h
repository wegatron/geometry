#ifndef SAMPLING_H
#define SAMPLING_H

#include <vector>
#include <zswlib/mesh/mesh_type.h>
#include "cgal_common.h"
#include "basic_data_structure.h"

namespace zsw
{
  void sampleTriangle(const Eigen::Matrix<zsw::Scalar, 3, 3> &tri_points, const zsw::Scalar r,
                      const zsw::Scalar val_exp, std::vector<zsw::JudgePoint> &samples);

  void sampleTet(const Eigen::Matrix<zsw::Scalar,3,4> &tet_points, const zsw::Scalar r,
                 std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &samples);

  class Sampler  final
  {
  public:
    Sampler() {}
    void sampleSigmaDense(const zsw::mesh::TriMesh &tm, const zsw::Scalar sigma, std::vector<Eigen::Matrix<zsw::Scalar, 3, 1>> &samples);
#ifdef ZSW_NDEBUG
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

    void projectToLine(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> &v1, Eigen::Matrix<zsw::Scalar,3,1> &sample_point);
  };

  template<typename SCALAR, size_t DIMENSION, size_t N>
    void calcBbox(const Eigen::Matrix<SCALAR, DIMENSION, N> &points, Eigen::Matrix<SCALAR, DIMENSION, 2> &bbox)
  {
    assert(points.rows()==DIMENSION && points.cols()==N);
    assert(bbox.rows()==DIMENSION);
    assert(N>0);
    bbox.template block<DIMENSION, 1>(0,0);
    bbox.template block<DIMENSION,1>(0,0)=points.template block<DIMENSION,1>(0,0);
    bbox.template block<DIMENSION,1>(0,1) = points.template block<DIMENSION,1>(0,0);
    for(size_t i=1; i<N; ++i) {
      for(size_t di=0; di<DIMENSION; ++di) {
        if(points(di,i) < bbox(di, 0)) {
          bbox(di,0) = points(di,i);
        } else if(points(di,i) > bbox(di,1)) {
          bbox(di,1) = points(di,i);
        }
      }
    }
  }
  bool inTet(const Point &pt, const Eigen::Matrix<zsw::Scalar,3,4> &tet_points);

  /***
   * check if vr and vt is in the same side of line v0, v1
   * all the point are in the same plane.
   */
  bool sameSide2D(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> &v1, const Eigen::Matrix<zsw::Scalar,3,1> &vr, const Eigen::Matrix<zsw::Scalar,3,1> &vt);

  bool sameSide3D(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> v1, const Eigen::Matrix<zsw::Scalar,3,1> &v2, const Eigen::Matrix<zsw::Scalar,3,1> &vr, const Eigen::Matrix<zsw::Scalar,3,1> &vt);

}



#endif /* SAMPLING_H */

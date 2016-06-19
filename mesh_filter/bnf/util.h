#ifndef UTIL_H
#define UTIL_H

#include <Eigen/Dense>

namespace zsw
{
  namespace mesh {
    /***
     * \brief 计算一个三角形的面积。
     * \param tri_pts 三个顶点的坐标。
     * \return 三角形的面积
     **/
    template<typename Scalar>
      Scalar calcFaceArea(const Eigen::Matrix<Scalar, 3, 3> &tri_pts)
      {
        Eigen::Matrix<Scalar, 3, 1> v[2] = {
          tri_pts.block<3,1>(0,1) - tri_pts.block<3,1>(0,0),
          tri_pts.block<3,1>(0,2) - tri_pts.block<3,1>(0,0)
        };
        return 0.5*v[0].cross(v[1]).norm();
      }

      /***
       * \brief 计算一个mesh的所有三角形面积。
       * \param trimesh 输入mesh
       * \return face_area 每个三角形的面积。
       * */
    template<typename Scalar>
      void calcFaceArea(const zsw::mesh::TriMesh &trimesh, Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &face_area)
      {
        face_area.resize(trimesh.n_faces(),1);
#pragma omp parallel for
        for(size_t i=0; i<trimesh.n_faces(); ++i) {
          Eigen::Matrix<Scalar, 3, 3> pts;
          size_t j=0;
          for(zsw::mesh::TriMesh::CFVIter fv_it = trimesh.cfv_iter(zsw::mesh::TriMesh::FaceHandle(i)); fv_it.is_valid();
              ++fv_it) {
            pts.block<3,1>(0,j) = trimesh.point(fv_it); ++j;
          }
          face_area(i) = calcFaceArea(pts);
        }
      }

      /***
       * \brief 计算mesh每个三角形面片的重心坐标。
       * \param trimesh 输入mesh
       * \return face_center 每个三角形的重心坐标。
       * */
    template<typename Scalar>
      void calcFaceCenter(const zsw::mesh::TriMesh &trimesh, Eigen::Matrix<Scalar, 3, Eigen::Dynamic> &face_center)
      {
        face_center.resize(3, trimesh.n_faces());
#pragma omp parallel for
        for(size_t i=0; i<trimesh.n_faces(); ++i) {
          face_center.block<3,1>(0,i) = trimesh.calc_face_centroid(zsw::mesh::TriMesh::FaceHandle(i));
        }
      }
  }
}

#endif /* UTIL_H */

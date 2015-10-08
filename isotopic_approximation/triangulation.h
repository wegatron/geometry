#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#include <vector>
#include <Eigen/Dense>
#include <zswlib/mesh/mesh_type.h>

#include "cgal_common.h"

namespace zsw {

  class KernelRegionJudger
  {
  public:
    KernelRegionJudger();
    void addConstraint(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> &v1,
                       const Eigen::Matrix<zsw::Scalar,3,1> &v2, const Eigen::Matrix<zsw::Scalar,3,1> &vr);
    bool judge(const Point &pt);
  private:
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> vec_v0;
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> vec_v1;
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> vec_v2;
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> vec_vr;
  };
  class TetMesh
  {
  public:
    struct TetPoint
    {
      char point_type_; // outer_boundary 1, zero_set 0, inner_boundary -1
      Point pt_data_;
      std::vector<size_t> tet_ids_;
    };

    struct Tet
    {
      size_t pt_ids_[4]; // sorted
      std::vector<Point> sample_points_;
    };

    TetMesh(const std::vector<Point> bz_points, const std::vector<Point> &bo_points, const std::vector<Point> &bi_points,
            const zsw::Scalar sample_dense);

    /*** clean invalid points
     * update tets and edges(keep invalid tets and edges)
     */
    void cleanPoints();
    void simplify();
    void writeVtk(const std::string &filepath);
    void writeZeroSetSurface(const std::string &filepath);
  private:
    void updateKrj(size_t tet_id);
    bool collapseZEdge(std::pair<size_t,size_t> &edge);
    bool isValidTet(Tet &tet);
    bool isValidEdge(std::pair<size_t, size_t> &edge);
    zsw::Scalar sample_dense_;
    std::vector<Tet> tets_;
    std::vector<TetPoint> tet_points_;
    std::vector<std::pair<size_t, size_t>> edges_;
    std::vector<size_t> pf_;
  };
}


#endif /* TRIANGULATION_H */

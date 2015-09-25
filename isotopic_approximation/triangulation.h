#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#include <vector>
#include <Eigen/Dense>
#include <zswlib/mesh/mesh_type.h>

#include "cgal_common.h"

namespace zsw {

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

    TetMesh(const std::vector<Point> bz_points, const std::vector<Point> &bo_points, const std::vector<Point> &bi_points);

    /*** clean invalid points
     * update tets and edges(keep invalid tets and edges)
     */
    void cleanPoints();
    void simplify();
    void writeVtk(const std::string &filepath);
    void writeZeroSetSurface(const std::string &filepath);
  private:
    bool collapseZEdge(std::pair<size_t,size_t> &edge);
    bool isValidTet(Tet &tet);
    bool isValidEdge(std::pair<size_t, size_t> &edge);
    std::vector<Tet> tets_;
    std::vector<TetPoint> tet_points_;
    std::vector<std::pair<size_t, size_t>> edges_;
    std::vector<size_t> pf_;
  };
}


#endif /* TRIANGULATION_H */

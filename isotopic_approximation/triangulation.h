#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#include <vector>
#include <Eigen/Dense>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <zswlib/mesh/mesh_type.h>

namespace zsw {
  typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
  typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>    Vb;
  typedef CGAL::Triangulation_data_structure_3<Vb>                    Tds;
  //Use the Fast_location tag. Default or Compact_location works too.
  typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Delaunay;
  typedef Delaunay::Point                                             Point;

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
    bool isValidTet(Tet &tet);
    bool isValidEdge(std::pair<size_t, size_t> &edge);
    std::vector<Tet> tets_;
    std::vector<TetPoint> tet_points_;
    std::vector<std::pair<size_t, size_t>> edges_;
    std::vector<size_t> pf_;
  };
}


#endif /* TRIANGULATION_H */

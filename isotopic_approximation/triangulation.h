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
      bool iszero_; // is in zero set
      Point pt_data_;
      std::vector<size_t> tet_ids_;
    };

    struct Tet
    {
      size_t pt_ids_[4];
      std::vector<Point> sample_points_;
    };

    TetMesh(const zsw::mesh::TriMesh &zero_mesh, const zsw::mesh::TriMesh &bo_mesh, const zsw::mesh::TriMesh &bi_mesh);
    void simplify();
    void writeVtk(const std::string &filepath);
    void writeZeroSetSurface(const std::string &filepath);
  private:
    std::vector<TetPoint> tet_points_;
    Eigen::Matrix<size_t, 2, 1> edges_;
    std::vector<size_t> pf_;
  };
}


#endif /* TRIANGULATION_H */

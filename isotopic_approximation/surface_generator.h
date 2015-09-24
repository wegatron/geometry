#ifndef SURFACE_GENERATOR_H
#define SURFACE_GENERATOR_H

#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <zswlib/mesh/mesh_type.h>

namespace zsw
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
  typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>    Vb;
  typedef CGAL::Triangulation_data_structure_3<Vb>                    Tds;
  //Use the Fast_location tag. Default or Compact_location works too.
  typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Delaunay;
  typedef Delaunay::Point                                             Point;

  class SurfaceGenerator
  {
  public:
    SurfaceGenerator() {}
    void genSurface(const zsw::Scalar dis, zsw::mesh::TriMesh &tm, zsw::mesh::TriMesh &bo_mesh, zsw::mesh::TriMesh &bi_mesh);
    void genPoints(const zsw::Scalar dis, std::vector<Delaunay::Point> &bz_points, std::vector<Delaunay::Point> &bo_points,
                   std::vector<Delaunay::Point> &bi_points);
  };
}


#endif /* SURFACE_GENERATOR_H */

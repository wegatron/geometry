#ifndef SURFACE_GENERATOR_H
#define SURFACE_GENERATOR_H

#include <vector>
#include <zswlib/mesh/mesh_type.h>
#include "cgal_common.h"

namespace zsw
{
  class SurfaceGenerator
  {
  public:
    SurfaceGenerator() {}
    void genSurface(const zsw::Scalar dis, zsw::mesh::TriMesh &tm, zsw::mesh::TriMesh &bo_mesh, zsw::mesh::TriMesh &bi_mesh);
    void genPoints(const zsw::Scalar dis, zsw::mesh::TriMesh &tm, std::vector<Point> &bz_points,
                   std::vector<Delaunay::Point> &bo_points,
                   std::vector<Delaunay::Point> &bi_points);
  };
}


#endif /* SURFACE_GENERATOR_H */

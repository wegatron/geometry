#ifndef SURFACE_GENERATOR_H
#define SURFACE_GENERATOR_H

#include <vector>
#include <zswlib/mesh/mesh_type.h>
#include "cgal_common.h"

namespace zsw
{
    void genPoints(const zsw::Scalar dis, zsw::mesh::TriMesh &tm, std::vector<Point> &bz_points,
                   std::vector<Delaunay::Point> &bo_points,
                   std::vector<Delaunay::Point> &bi_points);

  /* class SurfaceGenerator */
  /* { */
  /* public: */
  /*   /\* void genSurface(const zsw::Scalar dis, zsw::mesh::TriMesh &tm, zsw::mesh::TriMesh &bo_mesh, zsw::mesh::TriMesh &bi_mesh); *\/ */
  /*   void genPoints(const zsw::Scalar dis, zsw::mesh::TriMesh &tm, std::vector<Point> &bz_points, */
  /*                  std::vector<Delaunay::Point> &bo_points, */
  /*                  std::vector<Delaunay::Point> &bi_points); */

  /*   /\* void genPoints2(const zsw::Scalar dis, std::vector<Point> &bz_points, *\/ */
  /*   /\*                 std::vector<Delaunay::Point> &bo_points, *\/ */
  /*   /\*                 std::vector<Delaunay::Point> &bi_points); *\/ */
  /*   void empty(); */
  /* private: */
  /*   char ch; */
  /* }; */
}


#endif /* SURFACE_GENERATOR_H */

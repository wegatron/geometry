#ifndef SURFACE_GENERATOR_H
#define SURFACE_GENERATOR_H

#include <zswlib/mesh/mesh_type.h>

namespace zsw
{
  class SurfaceGenerator
  {
  public:
    SurfaceGenerator() {}
    void genSurface(const zsw::Scalar dis, zsw::mesh::TriMesh &tm, zsw::mesh::TriMesh &bo_mesh, zsw::mesh::TriMesh &bi_mesh);
  };
}


#endif /* SURFACE_GENERATOR_H */

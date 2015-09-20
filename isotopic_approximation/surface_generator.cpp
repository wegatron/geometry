#include "surface_generator.h"

void zsw::SurfaceGenerator::genSurface(const zsw::Scalar dis, zsw::mesh::TriMesh &tm,
                                        std::pair<zsw::mesh::TriMesh, zsw::mesh::TriMesh> &surf_mesh)
{
  if(!tm.has_vertex_normals()) {
    tm.request_vertex_normals();
    tm.update_vertex_normals();
  }
}

#include "surface_generator.h"

void zsw::SurfaceGenerator::genSurface(const zsw::Scalar dis, zsw::mesh::TriMesh &tm,
                                       zsw::mesh::TriMesh &bo_mesh, zsw::mesh::TriMesh &bi_mesh)
{
  if(!tm.has_vertex_normals()) {
    tm.request_vertex_normals();
    tm.update_vertex_normals();
  }
  bo_mesh  = tm;  bi_mesh  = tm;
  for(int i=0; i<tm.n_vertices(); ++i) {
    zsw::mesh::TriMesh::VertexHandle pt_h(i);
    bo_mesh.set_point(pt_h, tm.point(pt_h)+dis*tm.normal(pt_h));
    bi_mesh.set_point(pt_h, tm.point(pt_h)-dis*tm.normal(pt_h));
  }
}

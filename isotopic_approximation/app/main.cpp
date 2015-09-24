#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include "../triangulation.h"
#include "../surface_generator.h"

int main(int argc, char *argv[])
{
  zsw::mesh::TriMesh input_mesh;
  if(!OpenMesh::IO::read_mesh(input_mesh, "/home/wegatron/workspace/geometry/data/dragon.obj")) {
    std::cerr << "[ERROR] can't read mesh!" << std::endl;
  }
  zsw::SurfaceGenerator sfg;
  zsw::mesh::TriMesh bo_mesh, bi_mesh;
  sfg.genSurface(0.05, input_mesh, bo_mesh, bi_mesh);
  zsw::TetMesh tm(input_mesh, bo_mesh, bi_mesh);
  tm.simplify();
  tm.writeVtk("/home/wegatron/tmp.vtk");
  return 0;
}

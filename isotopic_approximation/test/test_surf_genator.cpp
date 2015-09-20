#include <iostream>
#include <zswlib/const_val.h>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "../surface_generator.h"

using namespace std;

void test_genSurface()
{
  zsw::mesh::TriMesh tm;
  if(!OpenMesh::IO::read_mesh(tm, "/home/wegatron/workspace/geometry/data/dragon.obj")) {
    std::cerr << "[ERROR] can't read mesh: /home/wegatron/workspace/geometry/data/dragon.obj"<< std::endl;
    return;
  }
  zsw::SurfaceGenerator surf_gen;
  zsw::mesh::TriMesh bo_mesh, bi_mesh;
  surf_gen.genSurface((zsw::Scalar)0.2, tm, bo_mesh, bi_mesh);
  OpenMesh::IO::write_mesh(bo_mesh, "/home/wegatron/tmp_bo.obj");
  // OpenMesh::IO::write_mesh(bi_mesh, "/home/wegatron/tmp_bi.obj");
}

int main(int argc, char *argv[])
{
  test_genSurface();
  return 0;
}

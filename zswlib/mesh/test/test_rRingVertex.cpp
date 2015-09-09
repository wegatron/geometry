#include <iostream>
#include <boost/foreach.hpp>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_op.h>

int main(int argc, char *argv[])
{
  zsw::mesh::TriMesh tm;
  if(!OpenMesh::IO::read_mesh(tm, "/home/wegatron/data/3dmodels/sphere_out.obj")) {
    std::cerr << "read mesh: " << "/home/wegatron/data/3dmodels/sphere_out.obj" << std::endl;
    return __LINE__;
  }
  std::vector<zsw::FakeSet<std::size_t>> ring;
  zsw::mesh::rRingVertex(tm,2, ring);
  for(size_t vid : ring[1346]) {
    std::cout << vid << std::endl;
  }
  return 0;
}

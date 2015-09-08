#include <iostream>
#include <boost/foreach.hpp>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_op.h>

int main(int argc, char *argv[])
{
  zsw::mesh::TriMesh tm;
  if(!OpenMesh::IO::read_mesh(tm, "E:/workspace/geometry/RBMLS/result/plane2.obj")) {
    std::cerr << "read mesh: " << "E:/workspace/geometry/RBMLS/result/plane2.obj" << std::endl;
    return __LINE__;
  }
  std::vector<zsw::FakeSet<size_t>> ring;
  zsw::mesh::rRingVertex(tm, ring, 5);
  BOOST_FOREACH(size_t vid, ring[4]) {
    std::cout << vid << " ";
  }
  return 0;
}

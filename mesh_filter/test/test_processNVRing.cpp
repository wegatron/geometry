#include <boost/foreach.hpp>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "../lh_filter.h"

#include "../basic_op.h"

int main(int argc, char *argv[])
{
  zsw::mesh::TriMesh trimesh;
  if(!OpenMesh::IO::read_mesh(trimesh, "E:/workspace/geometry/RBMLS/result/plane2.obj")) {
    std::cerr << "unable to read mesh!" << std::endl;
    return __LINE__;
  }
  std::vector<zsw::FakeSet<size_t>> ring;
  zsw::mesh::processNVRing(trimesh, 1, ring);

  BOOST_FOREACH(size_t i, ring[3]) {
    std::cout << i << " ";
  }
  return 0;
}

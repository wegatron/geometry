#include <iostream>
#include <boost/foreach.hpp>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "../bnf/bilateral_normal_filter.h"

using namespace std;

int main(int argc, char *argv[])
{
  zsw::mesh::TriMesh tm;
  if(!OpenMesh::IO::read_mesh(tm, "e:/tmp/input.stl")) {
    std::cerr << "read mesh: " << "e:/tmp/input.stl" << std::endl;
    return __LINE__;
  }
  std::vector<zsw::FakeSet<size_t>> ring;
  zsw::rRingFacesByVertex(tm, 1, ring);
  BOOST_FOREACH(size_t fid, ring[13]) {
    std::cerr << fid << " ";
  }
  return 0;
}

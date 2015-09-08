#include <iostream>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/EdgeFlipMean.h>

int main(int argc, char *argv[])
{
  zsw::mesh::TriMesh tm;
  OpenMesh::IO::read_mesh(tm, "E:/workspace/geometry/bilateral_normal_filtering/result/2n/input.obj");
  Sn3DGraphics::EdgeFlipMean ef;
  ef.run(tm);
  OpenMesh::IO::write_mesh(tm, "E:/tmp/test_edge_flip.obj");
  return 0;
}

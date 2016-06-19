#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include "../bnf/EdgeFlipMean.h"

int main(int argc, char *argv[])
{
  zsw::mesh::TriMesh tm;
  if(!OpenMesh::IO::read_mesh(tm, argv[1])) {
    std::cerr << "Failed to read mesh: " << std::string(argv[1]) << std::endl;
    return 0;
  }
  Sn3DGraphics::EdgeFlipMean ef;
  ef.run(tm);
  OpenMesh::IO::write_mesh(tm, argv[2]);
  return 0;
}

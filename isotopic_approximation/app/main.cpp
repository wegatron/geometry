#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include "../triangulation.h"
#include "../surface_generator.h"

using namespace std;

int main(int argc, char *argv[])
{
  zsw::mesh::TriMesh input_mesh;
  if(!OpenMesh::IO::read_mesh(input_mesh, "/home/wegatron/workspace/geometry/data/dragon.obj")) {
    std::cerr << "[ERROR] can't read mesh!" << std::endl;
  }
  zsw::SurfaceGenerator sfg;
  vector<zsw::Point> bz_points, bo_points, bi_points;
  sfg.genPoints(0.05, bz_points, bo_points, bi_points);
  zsw::TetMesh tm(bz_points, bo_points, bi_points);
  tm.simplify();
  tm.writeVtk("/home/wegatron/tmp.vtk");
  return 0;
}

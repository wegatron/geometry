#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include "../triangulation.h"
#include "../surface_generator.h"

using namespace std;

int main(int argc, char *argv[])
{
  zsw::mesh::TriMesh input_mesh;
  std::string file_path="/home/wegatron/workspace/geometry/data/sphere.stl";
  if(!OpenMesh::IO::read_mesh(input_mesh, file_path)) {
    std::cerr << "[ERROR] can't read mesh!" << std::endl;
  }
  zsw::SurfaceGenerator sfg;
  vector<zsw::Point> bz_points, bo_points, bi_points;
  sfg.genPoints(0.2, input_mesh, bz_points, bo_points, bi_points);
  zsw::TetMesh tm(bz_points, bo_points, bi_points, 0.02);
  tm.simplify();
  tm.writeVtk("/home/wegatron/tmp.vtk");
  return 0;
}

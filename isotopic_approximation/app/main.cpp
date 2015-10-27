#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include "../triangulation.h"
#include "../surface_generator.h"

using namespace std;

void test(const std::string &file_path, const zsw::Scalar dis)
{
  zsw::mesh::TriMesh input_mesh;
  if(!OpenMesh::IO::read_mesh(input_mesh, file_path)) {
    std::cerr << "[ERROR] can't read mesh!" << std::endl;
  }
  zsw::SurfaceGenerator sfg;
  vector<zsw::Point> bz_points, bo_points, bi_points;
  sfg.genPoints(dis, input_mesh, bz_points, bo_points, bi_points);
  zsw::TetMesh tm(bz_points, bo_points, bi_points, 0.02);
  tm.simplify();
  tm.writeVtk("/home/wegatron/tmp_beam.vtk");
}

int main(int argc, char *argv[])
{
  test("/home/wegatron/workspace/geometry/data/beam.stl", 0.1);
  return 0;
}

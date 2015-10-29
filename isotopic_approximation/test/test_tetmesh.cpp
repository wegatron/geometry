#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include "../triangulation.h"
#include "../surface_generator.h"

using namespace std;

void test_triangulation(const string &file_path, const string &output_prefix, const zsw::Scalar dis)
{
  zsw::mesh::TriMesh input_mesh;
  if(!OpenMesh::IO::read_mesh(input_mesh, file_path)) {
    std::cerr << "[ERROR] can't read mesh!" << std::endl;
  }

  vector<zsw::Point> bz_points, bo_points, bi_points;
  zsw::genPoints(dis, input_mesh, bz_points, bo_points, bi_points);

  zsw::TetMesh tm(bz_points, bo_points, bi_points);
  tm.writeVtk(output_prefix+"out.vtk");
  tm.writeVtk(output_prefix+"in.vtk", 1);
  std::cerr << "!!!!!!!" << std::endl;
}

int main(int argc, char *argv[])
{
  test_triangulation("/home/wegatron/workspace/geometry/data/teaport.obj",
                     "/home/wegatron/tmp/test_triangulation/teaport_", 0.01);

  test_triangulation("/home/wegatron/workspace/geometry/data/beam.stl",
                     "/home/wegatron/tmp/test_triangulation/beam_", 0.1);

  test_triangulation("/home/wegatron/workspace/geometry/data/monkey.stl",
                     "/home/wegatron/tmp/test_triangulation/monkey_", 0.03);
  return 0;
}

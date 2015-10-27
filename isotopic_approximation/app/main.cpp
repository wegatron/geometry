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
  tm.writeSurface("/home/wegatron/tmp/ori_bo_sphere.obj", 1);
  tm.writeSurface("/home/wegatron/tmp/ori_zero_sphere.obj", 0);
  tm.writeSurface("/home/wegatron/tmp/ori_bi_sphere.obj", -1);

  tm.simplify();
  tm.writeSurface("/home/wegatron/tmp/res_sphere.obj", 0);
}

int main(int argc, char *argv[])
{
  //test("/home/wegatron/workspace/geometry/data/beam.stl", 0.1);
  test("/home/wegatron/workspace/geometry/data/sphere.stl", 0.07);
  return 0;
}

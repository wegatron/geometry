#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include "../triangulation.h"
#include "../surface_generator.h"

using namespace std;

void test(const std::string &file_path, const string &output_prefix, const zsw::Scalar dis)
{
  zsw::mesh::TriMesh input_mesh;
  if(!OpenMesh::IO::read_mesh(input_mesh, file_path)) {
    std::cerr << "[ERROR] can't read mesh!" << std::endl;
  }
  zsw::SurfaceGenerator sfg;
  vector<zsw::Point> bz_points, bo_points, bi_points;
  sfg.genPoints(dis, input_mesh, bz_points, bo_points, bi_points);
  zsw::TetMesh tm(bz_points, bo_points, bi_points, 0.02);
  tm.writeSurface(output_prefix+"bi.obj", 1);
  tm.writeSurface(output_prefix+"zero.obj", 0);
  tm.writeSurface(output_prefix+"bi.obj", -1);

  tm.simplify(dis/100);
  tm.writeSurface(output_prefix+"res.obj", 0);
}

int main(int argc, char *argv[])
{
  //test("/home/wegatron/workspace/geometry/data/beam.stl", 0.1);
  //test("/home/wegatron/workspace/geometry/data/sphere.stl", 0.07);
  test("/home/wegatron/workspace/geometry/data/teaport.obj", "/home/wegatron/tmp/teaport/", 0.01);
  return 0;
}

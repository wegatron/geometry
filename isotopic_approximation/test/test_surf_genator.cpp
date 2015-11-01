#include <iostream>
#include <zswlib/const_val.h>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "../surface_generator.h"

using namespace std;

void test_genSurface()
{
  zsw::mesh::TriMesh tm;
  if(!OpenMesh::IO::read_mesh(tm, "/home/wegatron/workspace/geometry/data/sphere.obj")) {
    std::cerr << "[ERROR] can't read mesh: /home/wegatron/workspace/geometry/data/sphere.obj"<< std::endl;
    return;
  }

  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bo_points;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bi_points;
  zsw::genPoints(0.01, tm, bo_points, bi_points);
}

int main(int argc, char *argv[])
{
  test_genSurface();
  return 0;
}

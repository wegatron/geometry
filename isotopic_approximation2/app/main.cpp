#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include <zswlib/error_ctrl.h>
#include <zswlib/zsw_clock_c11.h>
#include "../isotopic_approximation.h"
#include "../surface_generator.h"

using namespace std;

int main(int argc, char *argv[])
{
  std::ofstream ofs;
  zsw::mesh::TriMesh input_mesh;
  if(!OpenMesh::IO::read_mesh(input_mesh, "/home/wegatron/workspace/geometry/data/sphere.obj")) {
    std::cerr << "[ERROR] can't read mesh!" << std::endl;
    abort();
  }

  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bi_points;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bo_points;
  zsw::genPoints(1.5, input_mesh, bi_points, bo_points);
  zsw::Approximation appro;
  appro.init(.05, 0.5, bi_points, bo_points);

  appro.writeJudgePoints("/home/wegatron/tmp/judge_point.vtk");
  appro.writeTetMesh("/home/wegatron/tmp/tmp_tol.vtk", {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
  appro.simpTolerance();
  return 0;
}

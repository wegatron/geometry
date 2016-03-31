#include "basic_data_structure.h"

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
  appro.testCollapse();
  appro.writeTetMesh("/home/wegatron/tmp/after_collapse.vtk", {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
  return 0;
}

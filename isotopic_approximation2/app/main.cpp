#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include <zswlib/error_ctrl.h>
#include <zswlib/zsw_clock_c11.h>

#include "../shell.h"
#include "../isotopic_approximation.h"

using namespace std;

void test0(const std::string &file_path, const zsw::Scalar err_epsilon,
           const zsw::Scalar tri_sample_r, const zsw::Scalar tet_sample_r)
{
  std::ofstream ofs;
  zsw::mesh::TriMesh input_mesh;
  if(!OpenMesh::IO::read_mesh(input_mesh, file_path)) {
    std::cerr << "can't open file " << file_path << std::endl;
    abort();
  }
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> inner_jpts;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> outer_jpts;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bs_jpts;
  zsw::genAndSampleShell(input_mesh, err_epsilon, tri_sample_r, inner_jpts, outer_jpts, bs_jpts);
  zsw::Approximation appro;
  appro.init(err_epsilon, tri_sample_r, tet_sample_r, inner_jpts, outer_jpts, bs_jpts);
  appro.writeTetMesh("/home/wegatron/tmp/refine_res.vtk", {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
}

int main(int argc, char *argv[])
{
  test0("/home/wegatron/workspace/geometry/data/sphere.obj",0.1, 0.03, 0.03);
  //test0("/home/wegatron/workspace/geometry/data/bunny.obj",1.5, 0.5, 0.5);
  return 0;
}

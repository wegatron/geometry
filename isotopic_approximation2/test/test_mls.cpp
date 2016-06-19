#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include <zswlib/error_ctrl.h>
#include <zswlib/zsw_clock_c11.h>

#include "../shell.h"
#include "../isotopic_approximation.h"

using namespace std;

void genDeformShell()
{
  // if(!OpenMesh::IO::read_mesh(ori_mesh, ori_file_path)) {
  //   std::cerr << "can't open file " << ori_file_path << std::endl;
  //   abort();
  // }
  // std::vector<Eigen::Matrix<zsw::Scalar,3,1>> inner_jpts;
  // std::vector<Eigen::Matrix<zsw::Scalar,3,1>> outer_jpts;
  // std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bs_jpts;
  // zsw::genAndSampleShell(ori_mesh, err_epsilon, tri_sample_r, inner_jpts, outer_jpts, bs_jpts);
  // zsw::mesh::TriMesh deformed_mesh;
  // if(!OpenMesh::IO::read_mesh(deformed_mesh, deformed_file_path)) {
  //   std::cerr << "can't open file " << deformed_file_path << std::endl;
  //   abort();
  // }
  // zsw::sampleMesh(ori_mesh, deformed_mesh, ori_sample_pts, deformed_sample_pts);
  // KdTreeWarper kd_tree; kd_tree.buildTree(ori_sample_pts[0].data(), ori_sample_pts.size());
  // std::vector<Eigen::Matrix<zsw::Scalar,3,1>> deformed_outer_jpts;
  // for(size_t i=0; i<outer_jpts.size(); ++i) {
  //   kd_tree.query()
  // }
}

int main(int argc, char *argv[])
{

  return 0;
}

#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include <zswlib/error_ctrl.h>
#include <zswlib/zsw_clock_c11.h>

#include "../shell.h"
#include "../isotopic_approximation.h"
using namespace std;

void test0(const std::string &file_path,
           const std::string &output_dir,
           const zsw::Scalar err_epsilon,
           const zsw::Scalar tri_sample_r,
           const zsw::Scalar tet_sample_r)
{
  zsw::mesh::TriMesh input_mesh;
  if(!OpenMesh::IO::read_mesh(input_mesh, file_path)) {
    std::cerr << "can't open file " << file_path << std::endl;
    abort();
  }
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> inner_jpts;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> outer_jpts;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bs_jpts;
  zsw::genAndSampleShell(input_mesh, err_epsilon, tri_sample_r, inner_jpts, outer_jpts, bs_jpts);
  zsw::writePoints(output_dir+"inner_jpts.vtk", inner_jpts);
  zsw::writePoints(output_dir+"outer_jpts.vtk", outer_jpts);
  zsw::writePoints(output_dir+"bs_jpts.vtk", bs_jpts);
  zsw::Approximation appro;
  appro.setTmpOutDir(output_dir);
  appro.init(err_epsilon, tri_sample_r, tet_sample_r, inner_jpts, outer_jpts, bs_jpts);
  appro.writeTetMesh(output_dir+"refine_res.vtk", {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
  appro.simpTolerance();
  appro.writeTetMesh("/home/wegatron/tmp/after_simp_tol.vtk", {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
  // appro.mutuallTessellation();
  // appro.writeTetMesh(output_dir+"mutuall_tessellation.vtk", {zsw::ignore_bbox, zsw::ignore_out});
  // appro.simpZeroSurface();
  // appro.writeTetMesh(output_dir+"simped_zero_surf.vtk", {zsw::ignore_bbox, zsw::ignore_out});
  // appro.simp(output_dir);
  appro.writeTetMesh(output_dir+"simped_final.vtk", {zsw::ignore_bbox, zsw::ignore_out});
}

int main(int argc, char *argv[])
{
  test0(std::string(argv[1]), std::string(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]));
  //test0("/home/wegatron/workspace/geometry/data/sphere.obj", "/home/wegatron/tmp/", 0.1, 0.03, 0.03);
  // test0("/home/wegatron/workspace/geometry/data/cylinder_smoothed.obj", "/home/wegatron/tmp/", 0.3, 0.1, 0.1);
  //test0("/home/wegatron/workspace/geometry/data/bunny.obj",1.5, 0.5, 0.5);
  return 0;
}

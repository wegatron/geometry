#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include <zswlib/error_ctrl.h>
#include <zswlib/zsw_clock_c11.h>

#include "../shell.h"
#include "../isotopic_approximation.h"

using namespace std;

void appro_nowv(
                const size_t version,
                const std::string &ori_file_path,
                const std::string &deformed_file_path,
                const std::string &output_dir,
                const zsw::Scalar err_epsilon,
                const zsw::Scalar tri_sample_r,
                const zsw::Scalar tet_sample_r)
{
  RESETLOG_FILE(output_dir+"log");
  zsw::mesh::TriMesh input_mesh;
  if(!OpenMesh::IO::read_mesh(input_mesh, ori_file_path)) {
    std::cerr << "can't open file " << ori_file_path << std::endl;
    abort();
  }
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> inner_jpts;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> outer_jpts;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bs_jpts;
  zsw::Scalar g_scale=1.0;
  zsw::mesh::TriMesh deformed_mesh;
  if(!OpenMesh::IO::read_mesh(deformed_mesh, deformed_file_path)) {
    std::cerr << "can't open file " << deformed_file_path << std::endl;
    abort();
  }
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> deformed_inner_jpts;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> deformed_outer_jpts;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> deformed_bs_jpts;
  zsw::genAndSampleShellD(input_mesh, deformed_mesh, err_epsilon, tri_sample_r, inner_jpts, outer_jpts, bs_jpts,
                          deformed_inner_jpts, deformed_outer_jpts, deformed_bs_jpts, g_scale);
  zsw::writePoints(output_dir+"deformed_inner_jpts.vtk", deformed_inner_jpts);
  zsw::writePoints(output_dir+"deformed_outer_jpts.vtk", deformed_outer_jpts);
  zsw::writePoints(output_dir+"deformed_bs_jpts.vtk", deformed_bs_jpts);
  zsw::writePoints(output_dir+"ori_inner_jpts.vtk", inner_jpts);
  zsw::writePoints(output_dir+"ori_outer_jpts.vtk", outer_jpts);
  zsw::writePoints(output_dir+"ori_bs_jpts.vtk", bs_jpts);
  zsw::Approximation appro;
  appro.setGscale(g_scale);
  appro.setTetByRef(version == 2);
  appro.setTmpOutDir(output_dir);
  appro.initD(err_epsilon, tri_sample_r, tet_sample_r, inner_jpts, outer_jpts, bs_jpts,
              deformed_inner_jpts, deformed_outer_jpts, deformed_bs_jpts);
  appro.writeTetMesh(output_dir+"refine_res.vtk", {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
  appro.simp(output_dir);
  appro.writeTetMesh(output_dir+"simped_final.vtk", {zsw::ignore_bbox, zsw::ignore_out});
}

void appro_oriv(const std::string &ori_file_path,
                const std::string &output_dir,
                const zsw::Scalar err_epsilon,
                const zsw::Scalar tri_sample_r,
                const zsw::Scalar tet_sample_r)
{
  RESETLOG_FILE(output_dir+"log");
  zsw::mesh::TriMesh input_mesh;
  if(!OpenMesh::IO::read_mesh(input_mesh, ori_file_path)) {
    std::cerr << "can't open file " << ori_file_path << std::endl;
    abort();
  }
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> inner_jpts;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> outer_jpts;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bs_jpts;
  zsw::genAndSampleShell(input_mesh, err_epsilon, tri_sample_r, inner_jpts, outer_jpts, bs_jpts);
  zsw::writePoints(output_dir+"ori_inner_jpts.vtk", inner_jpts);
  zsw::writePoints(output_dir+"ori_outer_jpts.vtk", outer_jpts);
  zsw::writePoints(output_dir+"ori_bs_jpts.vtk", bs_jpts);
  zsw::Approximation appro;
  appro.setTmpOutDir(output_dir);
  appro.init(err_epsilon, tri_sample_r, tet_sample_r, inner_jpts, outer_jpts, bs_jpts);
  appro.writeTetMesh(output_dir+"refine_res.vtk", {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
  appro.simp(output_dir);
  appro.writeTetMesh(output_dir+"simped_final.vtk", {zsw::ignore_bbox, zsw::ignore_out});
}

int main(int argc, char *argv[])
{
  size_t v = atoi(argv[1]);
  if(v == 0) {    appro_oriv(std::string(argv[2]), std::string(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]));  }
  else { appro_nowv(v, std::string(argv[2]), std::string(argv[3]), std::string(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7])); }
  return 0;
}

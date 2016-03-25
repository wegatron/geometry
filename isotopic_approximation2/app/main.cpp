#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include <zswlib/error_ctrl.h>
#include <zswlib/zsw_clock_c11.h>

#include "../shell.h"
#include "../isotopic_approximation.h"
#include "../deformer.h"
#include "../basic_op.h"
#include "../bound_sphere.h"

using namespace std;

void appro_v0(const std::string &ori_file_path,
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
  appro.writeZeroSurface(output_dir+"simped_final.obj");
}

void appro_v1(
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
  //appro.setTetByRef(version == 2);
  appro.setTmpOutDir(output_dir);
  appro.initD(err_epsilon, tri_sample_r, tet_sample_r, inner_jpts, outer_jpts, bs_jpts,
              deformed_inner_jpts, deformed_outer_jpts, deformed_bs_jpts);
  appro.writeTetMesh(output_dir+"refine_res.vtk", {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
  appro.simp(output_dir);
  appro.writeTetMesh(output_dir+"simp_final_tol.vtk", {zsw::ignore_bbox, zsw::ignore_self_out});
  appro.writeTetMesh(output_dir+"simp_final_deformed_tol.vtk", {zsw::ignore_bbox, zsw::ignore_self_out}, nullptr, false);
  appro.writeTetMesh(output_dir+"simped_final.vtk", {zsw::ignore_bbox, zsw::ignore_out});
  //appro.writeTetMesh(output_dir+"simped_final_deformed.vtk", {zsw::ignore_bbox, zsw::ignore_out}, nullptr, false);
  appro.writeZeroSurface(output_dir+"simped_final.obj");
}

void appro_v3(const std::string &ori_file_path,
                  const std::string &deformed_file_path,
                  const std::string &output_dir,
                  const zsw::Scalar err_epsilon,
                  const zsw::Scalar tri_sample_r,
                  const zsw::Scalar tet_sample_r,
                  const size_t near_count)
{
  RESETLOG_FILE(output_dir+"log");
  zsw::mesh::TriMesh mesh_ori;
  if(!OpenMesh::IO::read_mesh(mesh_ori, ori_file_path)) {
    std::cerr << "can't open file " << ori_file_path << std::endl;
    abort();
  }
  zsw::mesh::TriMesh mesh_d;
  if(!OpenMesh::IO::read_mesh(mesh_d, deformed_file_path)) {
    std::cerr << "can't open file " << deformed_file_path << std::endl;
    abort();
  }
  std::shared_ptr<std::vector<Eigen::Matrix<zsw::Scalar,3,1>>> samples_out(new std::vector<Eigen::Matrix<zsw::Scalar,3,1>>());
  std::shared_ptr<std::vector<Eigen::Matrix<zsw::Scalar,3,1>>> samples_in(new std::vector<Eigen::Matrix<zsw::Scalar,3,1>>());;
  std::shared_ptr<std::vector<Eigen::Matrix<zsw::Scalar,3,1>>> samples_outd(new std::vector<Eigen::Matrix<zsw::Scalar,3,1>>());;
  std::shared_ptr<std::vector<Eigen::Matrix<zsw::Scalar,3,1>>> samples_ind(new std::vector<Eigen::Matrix<zsw::Scalar,3,1>>());;

  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bs_jpts;
  zsw::genAndSampleShell(mesh_ori, err_epsilon, tri_sample_r, *samples_in, *samples_out);
  zsw::writePoints(output_dir+"ori_inner_jpts.vtk", *samples_in);
  zsw::writePoints(output_dir+"ori_outer_jpts.vtk", *samples_out);

  zsw::LocalVectorFieldDeformer deformer(mesh_ori, mesh_d, tri_sample_r, near_count, nullptr);
  deformer.deformTo(samples_out, samples_in, samples_outd, samples_ind);

  zsw::writePoints(output_dir+"deformed_inner_jpts.vtk", *samples_ind);
  zsw::writePoints(output_dir+"deformed_outer_jpts.vtk", *samples_outd);

  // bs_jpts
  Eigen::Matrix<zsw::Scalar,3,2> bbox;
  zsw::calcBBOX(*samples_outd, bbox);
  Eigen::Matrix<zsw::Scalar,3,1> center = (bbox.block<3,1>(0,0) + bbox.block<3,1>(0,1)) * 0.5;
  zsw::Scalar radius = (bbox.block<3,1>(0,0) - bbox.block<3,1>(0,1)).norm() * 0.5;
  zsw::BoundSphere bs("bound_sphere.obj", radius, center);
  bs_jpts = bs.getVertices();

  zsw::calcBBOX(*samples_out, bbox);
  zsw::Scalar radius_ori = (bbox.block<3,1>(0,0) - bbox.block<3,1>(0,1)).norm() * 0.5;
  zsw::Scalar gscale = radius / radius_ori;

  zsw::writePoints(output_dir+"bs_jpts.vtk", bs_jpts);
  zsw::Approximation appro;
  appro.setTmpOutDir(output_dir);
  appro.init(err_epsilon, gscale*tri_sample_r, gscale * tet_sample_r, *samples_ind, *samples_outd, bs_jpts);
  appro.writeTetMesh(output_dir+"refine_res.vtk", {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
  appro.simp(output_dir);
  appro.writeTetMesh(output_dir+"simped_final_d.vtk", {zsw::ignore_bbox, zsw::ignore_out});
  appro.writeZeroSurface(output_dir+"simped_final_d.obj");
  // ----------------------------------------------------------------------------------------------------------------------------------
  zsw::mesh::TriMesh mesh_fd;
  if(!OpenMesh::IO::read_mesh(mesh_fd, output_dir+"simped_final_d.obj")) {
    std::cout << __FILE__ << __LINE__ << std::endl; abort();
  }
  if(!mesh_fd.has_vertex_normals()) {
    mesh_fd.request_face_normals();
    mesh_fd.request_vertex_normals();
    mesh_fd.update_normals();
  }
  std::vector<zsw::Vector3s> pts, pts_bk, pt_vns;
  for(auto vit=mesh_fd.vertices_begin(); vit!=mesh_fd.vertices_end(); ++vit) {
    pts.push_back(mesh_fd.point(vit));
    pt_vns.push_back(mesh_fd.normal(vit));
  }
  pts_bk.resize(pts.size());
  deformer.deformBack(pts, pt_vns, pts_bk);
  size_t i=0;
  for(auto vit=mesh_fd.vertices_begin(); vit!=mesh_fd.vertices_end(); ++vit) { mesh_fd.set_point(vit, pts_bk[i++]); }
  OpenMesh::IO::write_mesh(mesh_fd, output_dir+"simped_final_o.obj");
}

int main(int argc, char *argv[])
{
  size_t v = atoi(argv[1]);
  if(v == 0) {    appro_v0(std::string(argv[2]), std::string(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]));  }
  else if(v == 1) { appro_v1(v, std::string(argv[2]), std::string(argv[3]), std::string(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7])); }
  else if(v == 3) {
    appro_v3(std::string(argv[2]), std::string(argv[3]), std::string(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]), atoi(argv[8]));
  } else {
    std::cout << __FILE__ << __LINE__ << std::endl;
  }
  PRINT_FUNCTION_COST();
  return 0;
}

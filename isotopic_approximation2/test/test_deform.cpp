#define BOOST_AUTO_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "../deformer.h"
#include <zswlib/mesh/mesh_type.h>
#include <zswlib/matrix_type.h>
#include <zswlib/mesh/vtk.h>

#include "../sampling.h"
#include "../shell.h"
#include "../isotopic_approximation.h"

using namespace std;

BOOST_AUTO_TEST_SUITE(vf_deform)

void deformBack(const zsw::Scalar err_epsilon, const zsw::Scalar tri_sample_r, const size_t near_count,
                const std::string &file_o, const std::string &file_d, const std::string &file_b,
                const std::string &file_prefix)
{
  zsw::mesh::TriMesh mesh_ori;
  if(!OpenMesh::IO::read_mesh(mesh_ori, file_o)) {
    std::cout << __FILE__ << __LINE__ << std::endl;
    abort();
  }
  zsw::mesh::TriMesh mesh_d;
  if(!OpenMesh::IO::read_mesh(mesh_d, file_d)) {
    std::cout << __FILE__ << __LINE__ << std::endl;
    abort();
  }
  std::cout << __FILE__ << __LINE__ << std::endl;
  std::shared_ptr<std::vector<Eigen::Matrix<zsw::Scalar,3,1>>> samples_out(new std::vector<Eigen::Matrix<zsw::Scalar,3,1>>());
  std::shared_ptr<std::vector<Eigen::Matrix<zsw::Scalar,3,1>>> samples_in(new std::vector<Eigen::Matrix<zsw::Scalar,3,1>>());;
  std::shared_ptr<std::vector<Eigen::Matrix<zsw::Scalar,3,1>>> samples_outd(new std::vector<Eigen::Matrix<zsw::Scalar,3,1>>());;
  std::shared_ptr<std::vector<Eigen::Matrix<zsw::Scalar,3,1>>> samples_ind(new std::vector<Eigen::Matrix<zsw::Scalar,3,1>>());;

  zsw::genAndSampleShell(mesh_ori, err_epsilon, tri_sample_r, *samples_in, *samples_out);
  zsw::writePoints(file_prefix+"_ori_inner_jpts.vtk", *samples_in);
  zsw::writePoints(file_prefix+"_ori_outer_jpts.vtk", *samples_out);

  std::cout << __FILE__ << __LINE__ << std::endl;

  zsw::LocalVectorFieldDeformer deformer(mesh_ori, mesh_d, tri_sample_r, near_count, nullptr);
  deformer.deformTo(samples_out, samples_in, samples_outd, samples_ind);

  std::cout << __FILE__ << __LINE__ << std::endl;

  zsw::writePoints(file_prefix+"_deformed_inner_jpts.vtk", *samples_ind);
  zsw::writePoints(file_prefix+"_deformed_outer_jpts.vtk", *samples_outd);
  // ----------------------------------------------------------------------------------------------------------------------------------
  zsw::mesh::TriMesh mesh_b;
  //if(!OpenMesh::IO::read_mesh(mesh_b, "/home/wegatron/tmp/simped_ellipsoid.ply")) {
  std::cout << "file_b:" << file_b << std::endl;
  if(!OpenMesh::IO::read_mesh(mesh_b, file_b)) {
    std::cout << __FILE__ << __LINE__ << std::endl;
    abort();
  }

  std::cout << __FILE__ << __LINE__ << std::endl;
  if(!mesh_b.has_vertex_normals()) {
    mesh_b.request_face_normals();
    mesh_b.request_vertex_normals();
    mesh_b.update_normals();
  }
  std::vector<zsw::Vector3s> pts, pts_bk, pt_vns;
  for(auto vit=mesh_b.vertices_begin(); vit!=mesh_b.vertices_end(); ++vit) {
    pts.push_back(mesh_b.point(vit));
    pt_vns.push_back(mesh_b.normal(vit));
  }
  pts_bk.resize(pts.size());
  deformer.deformBack(pts, pt_vns, pts_bk);
  std::cout << __FILE__ << __LINE__ << std::endl;
  size_t i=0;
  for(auto vit=mesh_b.vertices_begin(); vit!=mesh_b.vertices_end(); ++vit) { mesh_b.set_point(vit, pts_bk[i++]); }
  OpenMesh::IO::write_mesh(mesh_b, file_prefix+"_bk.obj");
}

BOOST_AUTO_TEST_CASE(ellipsoid_deform_back0)
{
  zsw::Scalar err_epsilon = 0.03;
  zsw::Scalar tri_sample_r = 0.015;
  zsw::Scalar tet_sample_r = 0.06;
  size_t near_count = 100;
  zsw::mesh::TriMesh mesh_ori;
  if(!OpenMesh::IO::read_mesh(mesh_ori, "/home/wegatron/workspace/geometry/data/ellipsoid.obj")) {
    std::cout << __FILE__ << __LINE__ << std::endl;
    abort();
  }
  zsw::mesh::TriMesh mesh_d;
  if(!OpenMesh::IO::read_mesh(mesh_d, "/home/wegatron/workspace/geometry/data/ellipsoid_deformed.1.obj")) {
    std::cout << __FILE__ << __LINE__ << std::endl;
    abort();
  }
  std::shared_ptr<std::vector<Eigen::Matrix<zsw::Scalar,3,1>>> samples_out(new std::vector<Eigen::Matrix<zsw::Scalar,3,1>>());
  std::shared_ptr<std::vector<Eigen::Matrix<zsw::Scalar,3,1>>> samples_in(new std::vector<Eigen::Matrix<zsw::Scalar,3,1>>());;
  std::shared_ptr<std::vector<Eigen::Matrix<zsw::Scalar,3,1>>> samples_outd(new std::vector<Eigen::Matrix<zsw::Scalar,3,1>>());;
  std::shared_ptr<std::vector<Eigen::Matrix<zsw::Scalar,3,1>>> samples_ind(new std::vector<Eigen::Matrix<zsw::Scalar,3,1>>());;

  zsw::genAndSampleShell(mesh_ori, err_epsilon, tri_sample_r, *samples_in, *samples_out);
  zsw::writePoints("/home/wegatron/tmp/deform_test/ori_inner_jpts.vtk", *samples_in);
  zsw::writePoints("/home/wegatron/tmp/deform_test/ori_outer_jpts.vtk", *samples_out);

  zsw::LocalVectorFieldDeformer deformer(mesh_ori, mesh_d, tri_sample_r, near_count, nullptr);
  deformer.deformTo(samples_out, samples_in, samples_outd, samples_ind);

  zsw::writePoints("/home/wegatron/tmp/deform_test/deformed_inner_jpts.vtk", *samples_ind);
  zsw::writePoints("/home/wegatron/tmp/deform_test/deformed_outer_jpts.vtk", *samples_outd);
  // ----------------------------------------------------------------------------------------------------------------------------------
  zsw::mesh::TriMesh mesh_b;
  //if(!OpenMesh::IO::read_mesh(mesh_b, "/home/wegatron/tmp/simped_ellipsoid.ply")) {
  if(!OpenMesh::IO::read_mesh(mesh_b, "/home/wegatron/workspace/geometry/data/ellipsoid_deformed.1.obj")) {
    std::cout << __FILE__ << __LINE__ << std::endl;
    abort();
  }

  if(!mesh_b.has_vertex_normals()) {
    mesh_b.request_face_normals();
    mesh_b.request_vertex_normals();
    mesh_b.update_normals();
  }
  std::vector<zsw::Vector3s> pts, pts_bk, pt_vns;
  for(auto vit=mesh_b.vertices_begin(); vit!=mesh_b.vertices_end(); ++vit) {
    pts.push_back(mesh_b.point(vit));
    pt_vns.push_back(mesh_b.normal(vit));
  }
  pts_bk.resize(pts.size());
  deformer.deformBack(pts, pt_vns, pts_bk);
  size_t i=0;
  for(auto vit=mesh_b.vertices_begin(); vit!=mesh_b.vertices_end(); ++vit) { mesh_b.set_point(vit, pts_bk[i++]); }
  OpenMesh::IO::write_mesh(mesh_b, "/home/wegatron/tmp/bk.obj");
}

BOOST_AUTO_TEST_CASE(ellipsoid_deform_back1)
{
  zsw::Scalar err_epsilon = 0.1;
  zsw::Scalar tri_sample_r = 0.05;
  size_t near_count = 100;
  std::string file_o = "/home/wegatron/workspace/geometry/data/ellipsoid.obj";
  std::string file_d = "/home/wegatron/workspace/geometry/data/ellipsoid_deformed.1.obj";
  std::string file_b = "/home/wegatron/tmp/approximate/ellipsoid_v3/simped_final_d.obj";
  std::string file_prefix = "/home/wegatron/tmp/deform_test/df1";
  deformBack(err_epsilon, tri_sample_r, near_count, file_o, file_d, file_b, file_prefix);
}

BOOST_AUTO_TEST_SUITE_END()

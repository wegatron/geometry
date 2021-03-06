#define BOOST_AUTO_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "../deformer.h"
#include <zswlib/mesh/mesh_type.h>
#include <zswlib/matrix_type.h>
#include <zswlib/mesh/vtk.h>

#include "../sampling.h"
#include "../isotopic_approximation.h"
#include "../kdtree_warp.h"

using namespace std;

void sampleMesh(const std::string &ori_file, const std::string &deformed_file,
                std::vector<zsw::Vector3s> &ori_samples, std::vector<zsw::Vector3s> &deformed_samples,
                std::vector<zsw::Vector3s> &ori_shell_vs, const zsw::Scalar err_epsilon,
                const zsw::Scalar sample_r)
{
  zsw::mesh::TriMesh ori_mesh;
  if(!OpenMesh::IO::read_mesh(ori_mesh, ori_file)) {
    std::cerr << "[ERROR] can't read ori mesh" << std::endl;
    abort();
  }
  if(!ori_mesh.has_vertex_normals()) {
    ori_mesh.request_face_normals();
    ori_mesh.request_vertex_normals();
    ori_mesh.update_normals();
  }
  zsw::mesh::TriMesh deformed_mesh;
  if(!OpenMesh::IO::read_mesh(deformed_mesh, deformed_file)) {
    std::cerr << "[ERROR] can't read deform mesh" << std::endl;
    abort();
  }
  Eigen::Matrix<zsw::Scalar,3,3> ori_tri, deformed_tri, shell_tri;
  const size_t nf=ori_mesh.n_faces();
  for(size_t fi=0; fi<nf; ++fi) {
    zsw::mesh::TriMesh::FaceHandle fh=zsw::mesh::TriMesh::FaceHandle(int(fi));
    size_t vi=0;
    for(zsw::mesh::TriMesh::CFVIter fv_it=ori_mesh.cfv_iter(fh); fv_it.is_valid(); ++fv_it) {
      Eigen::Matrix<zsw::Scalar,3,1> offset = ori_mesh.normal(*fv_it)*err_epsilon;
      ori_tri.block<3,1>(0,vi) = ori_mesh.point(*fv_it);
      shell_tri.block<3,1>(0, vi++) = ori_mesh.point(*fv_it) + offset;
    }
    vi=0;
    for(zsw::mesh::TriMesh::CFVIter fv_it=deformed_mesh.cfv_iter(fh); fv_it.is_valid(); ++fv_it) {
      deformed_tri.block<3,1>(0,vi++) = deformed_mesh.point(*fv_it);
    }
    zsw::sampleTriangleAndDeformedTriangle(deformed_tri, ori_tri, sample_r, deformed_samples, ori_samples);
    zsw::sampleTriangle(shell_tri, sample_r, ori_shell_vs);
  }
}

void sampleOutShell(zsw::mesh::TriMesh &ori_mesh, const zsw::Scalar err_epsilon,
                    std::vector<zsw::Vector3s> &ori_shell_vs, zsw::Scalar sample_r)
{
  if(!ori_mesh.has_vertex_normals()) {
    ori_mesh.request_face_normals();
    ori_mesh.request_vertex_normals();
    ori_mesh.update_normals();
  }
  ori_shell_vs.clear();
  const size_t nf=ori_mesh.n_faces();
  Eigen::Matrix<zsw::Scalar,3,3> shell_tri;
  for(size_t fi=0; fi<nf; ++fi) {
    zsw::mesh::TriMesh::FaceHandle fh=zsw::mesh::TriMesh::FaceHandle(int(fi));
    size_t vi=0;
    for(zsw::mesh::TriMesh::CFVIter fv_it=ori_mesh.cfv_iter(fh); fv_it.is_valid(); ++fv_it) {
      Eigen::Matrix<zsw::Scalar,3,1> offset = ori_mesh.normal(*fv_it)*err_epsilon;
      shell_tri.block<3,1>(0, vi++) = ori_mesh.point(*fv_it) + offset;
    }
    zsw::sampleTriangle(shell_tri, sample_r, ori_shell_vs);
  }
}

void deform(
            zsw::Deformer &deformer,
            zsw::mesh::TriMesh &ori_mesh,
            const std::string &prefix,
            const zsw::Scalar err_epsilon,
            const zsw::Scalar sample_r)
{
  std::vector<zsw::Vector3s> ori_shell_vs, deformed_shell_vs, deformed_back_vs;
  sampleOutShell(ori_mesh, err_epsilon, ori_shell_vs, sample_r);
  zsw::writePoints(prefix+"_shell_vs_ori.vtk", ori_shell_vs);
  std::cout << "ref_vs size = " << deformer.getRefVs().size() << std::endl;
  zsw::writePoints(prefix+"_ref_vs.vtk", deformer.getRefVs());
  zsw::writePoints(prefix+"_ref_dvs.vtk", deformer.getRefDvs());
  deformer.deformTo(ori_shell_vs, deformed_shell_vs);
  // // writeout sample points
  zsw::writePoints(prefix+"_res_dvs.vtk", deformed_shell_vs);
}

// BOOST_AUTO_TEST_SUITE(lt_deform)

// BOOST_AUTO_TEST_CASE(ellipsoid_xscale)
// {
//   const zsw::Scalar err_epsilon = 0.05;
//   const zsw::Scalar sample_r = 0.025;
//   zsw::mesh::TriMesh ori_mesh;
//   if(!OpenMesh::IO::read_mesh(ori_mesh, "/home/wegatron/workspace/geometry/data/ellipsoid.obj")) {
//     std::cerr << "[ERROR] can't read ori mesh" << std::endl;
//     abort();
//   }

//   zsw::mesh::TriMesh deformed_mesh;
//   if(!OpenMesh::IO::read_mesh(deformed_mesh, "/home/wegatron/workspace/geometry/data/ellipsoid_opt.obj")) {
//     std::cerr << "[ERROR] can't read deform mesh" << std::endl;
//     abort();
//   }

//   std::shared_ptr<zsw::DisWeightFunc> dwf(new zsw::GaussDisWeightFunc());
//   zsw::LocalTranslateDeformer deformer(ori_mesh, deformed_mesh, sample_r, dwf);
//   deform(deformer, ori_mesh, "/home/wegatron/tmp/deform_test/ltd_x_scale", err_epsilon, sample_r);
// }

// BOOST_AUTO_TEST_CASE(ellipsoid_allscale)
// {
//   const zsw::Scalar err_epsilon = 0.05;
//   const zsw::Scalar sample_r = 0.025;
//   zsw::mesh::TriMesh ori_mesh;
//   if(!OpenMesh::IO::read_mesh(ori_mesh, "/home/wegatron/workspace/geometry/data/ellipsoid.obj")) {
//     std::cerr << "[ERROR] can't read ori mesh" << std::endl;
//     abort();
//   }

//   zsw::mesh::TriMesh deformed_mesh;
//   if(!OpenMesh::IO::read_mesh(deformed_mesh, "/home/wegatron/workspace/geometry/data/ellipsoid_scale.obj")) {
//     std::cerr << "[ERROR] can't read deform mesh" << std::endl;
//     abort();
//   }

//   std::shared_ptr<zsw::DisWeightFunc> dwf(new zsw::GaussDisWeightFunc());
//   zsw::LocalTranslateDeformer deformer(ori_mesh, deformed_mesh, sample_r, dwf);
//   deform(deformer, ori_mesh, "/home/wegatron/tmp/deform_test/ltd_all_scale", err_epsilon, sample_r);
// }


// BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(vf_deform)

BOOST_AUTO_TEST_CASE(vf_jac_valid)
{
  const zsw::Scalar err_epsilon = 0.07;
  const zsw::Scalar sample_r = 0.07;
  zsw::mesh::TriMesh ori_mesh;
  if(!OpenMesh::IO::read_mesh(ori_mesh, "/home/wegatron/workspace/geometry/data/ball.obj")) {
    std::cerr << "[ERROR] can't read ori mesh" << std::endl;
    abort();
  }

  zsw::mesh::TriMesh deformed_mesh;
  if(!OpenMesh::IO::read_mesh(deformed_mesh, "/home/wegatron/workspace/geometry/data/ball_scale.obj")) {
    std::cerr << "[ERROR] can't read deform mesh" << std::endl;
    abort();
  }

  std::shared_ptr<zsw::DisWeightFunc> dwf(new zsw::GaussDisWeightFunc());
  zsw::LocalVectorFieldDeformer deformer(ori_mesh, deformed_mesh, sample_r, 100, dwf);
  const std::vector<Eigen::Matrix<zsw::Scalar,3,3>> &jac = deformer.getJac();
  Eigen::Matrix<zsw::Scalar,3,3> exp_jac = 2 * Eigen::Matrix<zsw::Scalar,3,3>::Identity();
  deform(deformer, ori_mesh, "/home/wegatron/tmp/deform_test/lvfd_ball_scale", err_epsilon, sample_r);
  for(size_t i=0; i<jac.size(); ++i) {
    BOOST_CHECK(jac[i].allFinite());
    BOOST_CHECK((jac[i]- exp_jac).norm() < 1e-3);
    if((jac[i]- exp_jac).norm() > 1e-3) { std::cout << "jac:\n" << jac[i] << std::endl; }
  }
}

BOOST_AUTO_TEST_CASE(vf_jac_basic_valid_jac)
{
  const zsw::Scalar err_epsilon = 0.05;
  const zsw::Scalar sample_r = 0.05;
  zsw::mesh::TriMesh ori_mesh;
  if(!OpenMesh::IO::read_mesh(ori_mesh, "/home/wegatron/workspace/geometry/data/ball2.obj")) {
    std::cerr << "[ERROR] can't read ori mesh" << std::endl;
    abort();
  }

  zsw::mesh::TriMesh deformed_mesh;
  if(!OpenMesh::IO::read_mesh(deformed_mesh, "/home/wegatron/workspace/geometry/data/ball_scale2.obj")) {
    std::cerr << "[ERROR] can't read deform mesh" << std::endl;
    abort();
  }

  std::shared_ptr<zsw::DisWeightFunc> dwf(new zsw::GaussDisWeightFunc());
  zsw::LocalVectorFieldDeformer deformer(ori_mesh, deformed_mesh, sample_r, 100, dwf);
  const std::vector<Eigen::Matrix<zsw::Scalar,3,3>> &jac = deformer.getJac();
  Eigen::Matrix<zsw::Scalar,3,3> exp_jac = 2 * Eigen::Matrix<zsw::Scalar,3,3>::Identity();
  deform(deformer, ori_mesh, "/home/wegatron/tmp/deform_test/lvfd_ball_scale2", err_epsilon, sample_r);
  for(size_t i=0; i<jac.size(); ++i) {
    BOOST_CHECK(jac[i].allFinite());
    BOOST_CHECK((jac[i]- exp_jac).norm() < 1e-3);
    if((jac[i]- exp_jac).norm() > 1e-3) { std::cout << "jac:\n" << jac[i] << std::endl; }
  }
}

BOOST_AUTO_TEST_CASE(vf_jac_basic_resolve_invalid_jac)
{
  const zsw::Scalar err_epsilon = 0.03;
  const zsw::Scalar sample_r = 0.02;
  zsw::mesh::TriMesh ori_mesh;
  if(!OpenMesh::IO::read_mesh(ori_mesh, "/home/wegatron/workspace/geometry/data/cube.obj")) {
    std::cerr << "[ERROR] can't read ori mesh" << std::endl;
    abort();
  }

  zsw::mesh::TriMesh deformed_mesh;
  if(!OpenMesh::IO::read_mesh(deformed_mesh, "/home/wegatron/workspace/geometry/data/cube2.obj")) {
    std::cerr << "[ERROR] can't read deform mesh" << std::endl;
    abort();
  }

  std::shared_ptr<zsw::DisWeightFunc> dwf(new zsw::GaussDisWeightFunc());
  zsw::LocalVectorFieldDeformer deformer(ori_mesh, deformed_mesh, sample_r, 100, dwf);
  const std::vector<Eigen::Matrix<zsw::Scalar,3,3>> &jac = deformer.getJac();
  Eigen::Matrix<zsw::Scalar,3,3> exp_jac = 2*Eigen::Matrix<zsw::Scalar,3,3>::Identity();
  for(size_t i=0; i<jac.size(); ++i) {
    // BOOST_CHECK(jac[i].allFinite());
    if(!jac[i].allFinite()) { continue; }
    BOOST_CHECK((jac[i]- exp_jac).norm() < 1e-3);
    if((jac[i]- exp_jac).norm() > 1e-3) {
      std::cout << "jac:\n" << jac[i] << std::endl;
    }
  }
}

BOOST_AUTO_TEST_CASE(vf_cube_scale)
{
  const zsw::Scalar err_epsilon = 0.03;
  const zsw::Scalar sample_r = 0.02;
  zsw::mesh::TriMesh ori_mesh;
  if(!OpenMesh::IO::read_mesh(ori_mesh, "/home/wegatron/workspace/geometry/data/cube.obj")) {
    std::cerr << "[ERROR] can't read ori mesh" << std::endl;
    abort();
  }

  zsw::mesh::TriMesh deformed_mesh;
  if(!OpenMesh::IO::read_mesh(deformed_mesh, "/home/wegatron/workspace/geometry/data/cube2.obj")) {
    std::cerr << "[ERROR] can't read deform mesh" << std::endl;
    abort();
  }

  std::shared_ptr<zsw::DisWeightFunc> dwf(new zsw::GaussDisWeightFunc());
  zsw::LocalVectorFieldDeformer deformer(ori_mesh, deformed_mesh, sample_r, 130, dwf);
  const std::vector<Eigen::Matrix<zsw::Scalar,3,3>> &jac = deformer.getJac();
  for(size_t i=0; i<jac.size(); ++i) {
    BOOST_REQUIRE(jac[i].allFinite());
  }
  deform(deformer, ori_mesh, "/home/wegatron/tmp/deform_test/lvfd_cube_scale", err_epsilon, sample_r);
}

BOOST_AUTO_TEST_CASE(vf_banana_scale)
{
  const zsw::Scalar err_epsilon = 0.06;
  const zsw::Scalar sample_r = 0.04;
  zsw::mesh::TriMesh ori_mesh;
  if(!OpenMesh::IO::read_mesh(ori_mesh, "/home/wegatron/workspace/geometry/data/banana.obj")) {
    std::cerr << "[ERROR] can't read ori mesh" << std::endl;
    abort();
  }

  zsw::mesh::TriMesh deformed_mesh;
  if(!OpenMesh::IO::read_mesh(deformed_mesh, "/home/wegatron/workspace/geometry/data/banana2.obj")) {
    std::cerr << "[ERROR] can't read deform mesh" << std::endl;
    abort();
  }

  std::shared_ptr<zsw::DisWeightFunc> dwf(new zsw::GaussDisWeightFunc());
  zsw::LocalVectorFieldDeformer deformer(ori_mesh, deformed_mesh, sample_r, 130, dwf);
  const std::vector<Eigen::Matrix<zsw::Scalar,3,3>> &jac = deformer.getJac();
  for(size_t i=0; i<jac.size(); ++i) {
    BOOST_REQUIRE(jac[i].allFinite());
  }
  deform(deformer, ori_mesh, "/home/wegatron/tmp/deform_test/lvfd_banana_scale", err_epsilon, sample_r);
}

BOOST_AUTO_TEST_CASE(vf_banana_deform)
{
  const zsw::Scalar err_epsilon = 0.06;
  const zsw::Scalar sample_r = 0.04;
  zsw::mesh::TriMesh ori_mesh;
  if(!OpenMesh::IO::read_mesh(ori_mesh, "/home/wegatron/workspace/geometry/data/banana.obj")) {
    std::cerr << "[ERROR] can't read ori mesh" << std::endl;
    abort();
  }

  zsw::mesh::TriMesh deformed_mesh;
  if(!OpenMesh::IO::read_mesh(deformed_mesh, "/home/wegatron/workspace/geometry/data/banana_deformed.obj")) {
    std::cerr << "[ERROR] can't read deform mesh" << std::endl;
    abort();
  }

  std::shared_ptr<zsw::DisWeightFunc> dwf(new zsw::GaussDisWeightFunc());
  zsw::LocalVectorFieldDeformer deformer(ori_mesh, deformed_mesh, sample_r, 130, dwf);
  const std::vector<Eigen::Matrix<zsw::Scalar,3,3>> &jac = deformer.getJac();
  for(size_t i=0; i<jac.size(); ++i) {
    BOOST_REQUIRE(jac[i].allFinite());
  }
  deform(deformer, ori_mesh, "/home/wegatron/tmp/deform_test/lvfd_banana_deformed", err_epsilon, sample_r);
}

BOOST_AUTO_TEST_CASE(ellipsoid_xscale)
{
  const zsw::Scalar err_epsilon = 0.05;
  const zsw::Scalar sample_r = 0.025;
  zsw::mesh::TriMesh ori_mesh;
  if(!OpenMesh::IO::read_mesh(ori_mesh, "/home/wegatron/workspace/geometry/data/ellipsoid.obj")) {
    std::cerr << "[ERROR] can't read ori mesh" << std::endl;
    abort();
  }

  zsw::mesh::TriMesh deformed_mesh;
  if(!OpenMesh::IO::read_mesh(deformed_mesh, "/home/wegatron/workspace/geometry/data/ellipsoid_opt.obj")) {
    std::cerr << "[ERROR] can't read deform mesh" << std::endl;
    abort();
  }

  std::shared_ptr<zsw::DisWeightFunc> dwf(new zsw::GaussDisWeightFunc());
  zsw::LocalVectorFieldDeformer deformer(ori_mesh, deformed_mesh, sample_r, 100, dwf);
  deform(deformer, ori_mesh, "/home/wegatron/tmp/deform_test/lvfd_x_scale", err_epsilon, sample_r);
}

BOOST_AUTO_TEST_CASE(ellipsoid_allscale)
{
  const zsw::Scalar err_epsilon = 0.05;
  const zsw::Scalar sample_r = 0.025;
  zsw::mesh::TriMesh ori_mesh;
  if(!OpenMesh::IO::read_mesh(ori_mesh, "/home/wegatron/workspace/geometry/data/ellipsoid.obj")) {
    std::cerr << "[ERROR] can't read ori mesh" << std::endl;
    abort();
  }

  zsw::mesh::TriMesh deformed_mesh;
  if(!OpenMesh::IO::read_mesh(deformed_mesh, "/home/wegatron/workspace/geometry/data/ellipsoid_scale.obj")) {
    std::cerr << "[ERROR] can't read deform mesh" << std::endl;
    abort();
  }

  std::shared_ptr<zsw::DisWeightFunc> dwf(new zsw::GaussDisWeightFunc());
  zsw::LocalVectorFieldDeformer deformer(ori_mesh, deformed_mesh, sample_r, 100, dwf);
  deform(deformer, ori_mesh, "/home/wegatron/tmp/deform_test/lvfd_ellipsoid_scale", err_epsilon, sample_r);
}

BOOST_AUTO_TEST_CASE(ellipsoid_deform)
{
  const zsw::Scalar err_epsilon = 0.05;
  const zsw::Scalar sample_r = 0.025;
  zsw::mesh::TriMesh ori_mesh;
  if(!OpenMesh::IO::read_mesh(ori_mesh, "/home/wegatron/workspace/geometry/data/ellipsoid.obj")) {
    std::cerr << "[ERROR] can't read ori mesh" << std::endl;
    abort();
  }

  zsw::mesh::TriMesh deformed_mesh;
  if(!OpenMesh::IO::read_mesh(deformed_mesh, "/home/wegatron/workspace/geometry/data/ellipsoid_deformed.1.obj")) {
    std::cerr << "[ERROR] can't read deform mesh" << std::endl;
    abort();
  }

  std::shared_ptr<zsw::DisWeightFunc> dwf(new zsw::GaussDisWeightFunc());
  zsw::LocalVectorFieldDeformer deformer(ori_mesh, deformed_mesh, sample_r, 100, nullptr);
  deform(deformer, ori_mesh, "/home/wegatron/tmp/deform_test/lvfd_ellipsoid_deform", err_epsilon, sample_r);
}

BOOST_AUTO_TEST_SUITE_END()

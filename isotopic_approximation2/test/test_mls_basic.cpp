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

BOOST_AUTO_TEST_SUITE(test_calc_deform)

BOOST_AUTO_TEST_CASE(plane_scale_deform)
{
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> vs, dvs;

  for(size_t i=0; i<=100; ++i) {
    for(size_t j=0; j<=100; ++j) {
      Eigen::Matrix<zsw::Scalar,3,1> v, dv;
      v << i*0.01, j*0.01, 0;
      dv << i*0.01, j*0.005, 0;
      vs.push_back(v); dvs.push_back(dv);
    }
  }

  Eigen::Matrix<zsw::Scalar,3,1> vt; vt << 0, 0.5, 1;
  zsw::LocalTranslateDeformFunc def;
  Eigen::Matrix<zsw::Scalar,3,1> dvt;
  def.calcDeformedPos(vs, dvs, vt, dvt);
  std::cout << "deformed dvt:" << dvt.transpose() << std::endl;
}

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

BOOST_AUTO_TEST_CASE(ellipsoid_case)
{
  std::vector<zsw::Vector3s> ori_samples, deformed_samples, ori_shell_vs, deformed_shell_vs;
  const zsw::Scalar err_epsilon = 0.05;
  const zsw::Scalar sample_r = 0.025;
  sampleMesh("/home/wegatron/workspace/geometry/data/ellipsoid.obj",
             "/home/wegatron/workspace/geometry/data/ellipsoid_opt.obj",
             ori_samples, deformed_samples,
             ori_shell_vs, err_epsilon,
             sample_r);
  // writeout sample points
  zsw::writePoints("/home/wegatron/tmp/deform_test/ori_samples.vtk", ori_samples);
  zsw::writePoints("/home/wegatron/tmp/deform_test/deformed_samples.vtk", deformed_samples);
  zsw::writePoints("/home/wegatron/tmp/deform_test/ori_shell_vs.vtk", ori_shell_vs);
  // build kdTree
  zsw::KdTreeWarper ori_kdt; ori_kdt.buildTree(ori_samples[0].data(), ori_samples.size());
  //zsw::KdTreeWarper deformed_kdt; deformed_kdt.buildTree(deformed_samples[0].data(), deformed_samples.size());
  std::vector<zsw::KdTreeNode> vs_in_r;
  std::vector<zsw::Vector3s> vs, dvs;
  std::vector<zsw::Scalar> errs(ori_shell_vs.size(), 0);
  deformed_shell_vs.resize(ori_shell_vs.size());
  zsw::LocalTranslateDeformFunc df;
  // for each vertex query the near r pts
  for(size_t i=0; i<ori_shell_vs.size(); ++i) {
    vs_in_r.clear();
    ori_kdt.findWithinR(ori_shell_vs[i], 0.1, std::back_inserter(vs_in_r));
    vs.clear(); dvs.clear();
    for(const zsw::KdTreeNode &node : vs_in_r) {
      vs.push_back(ori_samples[node.index_]);
      dvs.push_back(deformed_samples[node.index_]);
    }
    df.calcDeformedPos(vs, dvs, ori_shell_vs[i], deformed_shell_vs[i]);
    zsw::Vector3s expected = ori_shell_vs[i]; expected[0] = expected[0]*0.25;
    errs[i] = (deformed_shell_vs[i] - expected).norm();
  }
  // writeout
  std::ofstream ofs("/home/wegatron/tmp/deform_test/deformed_shell_vs.vtk");
  std::vector<size_t> points_ind(ori_shell_vs.size(), 0);
  for(size_t i=0; i<ori_shell_vs.size(); ++i) { points_ind[i]=i; }
  point2vtk(ofs, deformed_shell_vs[0].data(), deformed_shell_vs.size(), points_ind.data(), points_ind.size());
  point_data(ofs, errs.begin(), errs.size(), "err");
  //zsw::writePoints("/home/wegatron/tmp/deform_test/deformed_shell_vs.vtk", deformed_shell_vs);
}

BOOST_AUTO_TEST_SUITE_END()

#include <fstream>
#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include<zswlib/mesh/mesh_type.h>
#include <zswlib/mesh/vtk.h>
#include <zswlib/zsw_log.h>
#include <zswlib/zsw_clock_c11.h>
#include <zswlib/const_val.h>
#include <zswlib/error_ctrl.h>
#include "../surface_generator.h"
#include "../triangulation2.h"
#include "../debug.h"


using namespace std;

bool adjEdge(const zsw::Tet &tet, const size_t v0, const size_t v1)
{
  size_t vcnt=0;
  for(size_t i=0; i<4; ++i) {
    if(tet.vid_[i] == v0 || tet.vid_[i]==v1) {      ++vcnt;    }
  }
  return vcnt!=2;
}

void test_init(const std::string &file_path, const std::string &output_prefix,
               const zsw::Scalar thick, const zsw::Scalar r)
{
  // input obj mesh
  zsw::mesh::TriMesh input_mesh;
  if(!OpenMesh::IO::read_mesh(input_mesh, file_path)) {
    std::cerr << "[ERROR] can't read mesh!" << std::endl;
  }
  std::cerr << "input mesh vertex num:" << input_mesh.n_vertices() << std::endl;

  // gen bi and bo points
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bo_points;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bi_points;
  zsw::genPoints(thick, input_mesh, bo_points, bi_points);
  // init triangulation

  std::cerr << "start triangulation init!!!" << std::endl;
  zsw::common::ClockC11 clock;
  zsw::Triangulation tr(0.3,"/home/wegatron/tmp/");
  CALL_FUNC(tr.construct(0.25, r, 0.3, bo_points, bi_points), abort());

  std::cerr << "init triangulation time cost: " << clock.time() << std::endl;
  // output triangulation
  std::function<bool(const zsw::Tet&)> ignore_bbox
    = std::bind(&zsw::Triangulation::ignoreWithPtType, &tr, std::placeholders::_1, zsw::BBOX_POINT);
  std::function<bool(const zsw::Tet&)> ignore_self_out
    = std::bind(&zsw::Triangulation::ignoreOnlyWithPtType, &tr, std::placeholders::_1, zsw::OUTER_POINT);
  // std::function<bool(const zsw::Tet&)> adj_edge
  //   = std::bind(adjEdge, std::placeholders::_1, 67, 56);

  // std::function<bool(const zsw::Tet&)> adj_edge2
  //   = std::bind(adjEdge, std::placeholders::_1, 50, 65);

  tr.writeTetMesh(output_prefix+"_tol.vtk", {ignore_bbox, ignore_self_out});
  std::vector<zsw::Vertex> &vertices = tr.getVertices();
  std::vector<zsw::Tet> &tets = tr.getTets();
  std::vector<size_t> sliver_tet_ids;
  for(size_t i=0; i<tets.size(); ++i) {
    const zsw::Tet &tet = tets[i];
    const Eigen::Matrix<zsw::Scalar,3,1> &v0 = vertices[tet.vid_[0]].pt_;
    const Eigen::Matrix<zsw::Scalar,3,1> &v1 = vertices[tet.vid_[1]].pt_;
    const Eigen::Matrix<zsw::Scalar,3,1> &v2 = vertices[tet.vid_[2]].pt_;
    const Eigen::Matrix<zsw::Scalar,3,1> &v3 = vertices[tet.vid_[3]].pt_;
    Eigen::Matrix<zsw::Scalar,3,1> va = v1 - v0;
    Eigen::Matrix<zsw::Scalar,3,1> vb = v2 - v0;
    Eigen::Matrix<zsw::Scalar,3,1> vc = v3 - v0;
    Eigen::Matrix<zsw::Scalar,3,1> vn = va.cross(vb); vn.normalize();
    if(fabs(vc.dot(vn)) < 10*zsw::const_val::eps) {
      tr.writeTet("/home/wegatron/tmp/simp_tol/debug/slivertet_"+std::to_string(i)+".vtk", i);
    }
  }

}

void test_mutualTessellation()
{
  // input obj mesh
  const string file_path="/home/wegatron/workspace/geometry/data/cube_2.obj";
  zsw::mesh::TriMesh input_mesh;
  if(!OpenMesh::IO::read_mesh(input_mesh, file_path)) {
    std::cerr << "[ERROR] can't read mesh!" << std::endl;
  }

  // gen bi and bo points
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bo_points;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bi_points;
  zsw::genPoints(0.5, input_mesh, bo_points, bi_points);
  // init triangulation
  zsw::Triangulation tr(0.3,"/home/wegatron/tmp/");
  CALL_FUNC(tr.construct(0.25, 0.1, 0.3, bo_points, bi_points), abort());

  tr.mutualTessellation();
  // output triangulation
  std::function<bool(const zsw::Tet&)> ignore_bbox
    = std::bind(&zsw::Triangulation::ignoreWithPtType, &tr, std::placeholders::_1, zsw::BBOX_POINT);
  std::function<bool(const zsw::Tet&)> ignore_out
    = std::bind(&zsw::Triangulation::ignoreWithPtType, &tr, std::placeholders::_1, zsw::OUTER_POINT);

  tr.writeTetMesh("/home/wegatron/tmp/cube_mutual.vtk", {ignore_bbox, ignore_out});
}

int main(int argc, char *argv[])
{
  // test_init("/home/wegatron/workspace/geometry/data/sphere.stl",
  //           "/home/wegatron/tmp/test_triangulation/tr", 0.1, 0.3);
  test_init("/home/wegatron/workspace/geometry/data/fandisk.obj",
            "/home/wegatron/tmp/test_triangulation/tr", 0.008, 0.004);
  // test_init2();
  //test_mutualTessellation();
  return 0;
}

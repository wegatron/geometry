#include <fstream>
#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include<zswlib/mesh/mesh_type.h>
#include <zswlib/mesh/vtk.h>
#include <zswlib/zsw_log.h>
#include <zswlib/zsw_clock_c11.h>
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
  zsw::Triangulation tr(r, bo_points, bi_points);
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
  // tr.writeTetMesh(output_prefix+"_adj.vtk", {adj_edge});
  // tr.writeTetMesh(output_prefix+"_adj2.vtk", {adj_edge2});
  //zsw::writeJudgePoints(output_prefix, tr.getJpts());
  // output edges
  // const vector<zsw::Edge>& edges = triangulation.getEdges();
  // for(const zsw::Edge &e : edges) {
  //   std::cerr << "e: " << e.vid_[0] << " " << e.vid_[1] << std::endl;
  // }
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
  zsw::Triangulation triangulation(0.1, bo_points, bi_points);
  triangulation.mutualTessellation();
  // output triangulation
  std::function<bool(const zsw::Tet&)> ignore_bbox
    = std::bind(&zsw::Triangulation::ignoreWithPtType, &triangulation, std::placeholders::_1, zsw::BBOX_POINT);
  std::function<bool(const zsw::Tet&)> ignore_out
    = std::bind(&zsw::Triangulation::ignoreWithPtType, &triangulation, std::placeholders::_1, zsw::OUTER_POINT);

  triangulation.writeTetMesh("/home/wegatron/tmp/cube_mutual.vtk", {ignore_bbox, ignore_out});
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

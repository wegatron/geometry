#include <fstream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include<zswlib/mesh/mesh_type.h>
#include <zswlib/mesh/vtk.h>
#include "../surface_generator.h"
#include "../triangulation2.h"

using namespace std;

void test0()
{
  // input obj mesh
  const string file_path="/home/wegatron/workspace/geometry/data/cube.obj";
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
  // output triangulation
  triangulation.writeTetMesh("/home/wegatron/tmp/cube1.vtk", 0);
  // output judge points
  const string jp_file="/home/wegatron/tmp/cube_jp.obj";
  ofstream ofs(jp_file, std::ofstream::out) ;
  const vector<zsw::Tet>& tets = triangulation.getTets();
  for(const zsw::Tet& tet : tets) {
    for(const zsw::JudgePoint &jp : tet.jpts_) {
      ofs << "v " << jp.pt_[0] << " " << jp.pt_[1] << " " << jp.pt_[2] << std::endl;
    }
  }
  // output edges
  const vector<zsw::Edge>& edges = triangulation.getEdges();
  for(const zsw::Edge &e : edges) {
    std::cerr << "e: " << e.vid_[0] << " " << e.vid_[1] << std::endl;
  }
}

void test01()
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
  // output triangulation
  triangulation.writeTetMesh("/home/wegatron/tmp/cube01.vtk", zsw::BBOX_POINT);
  // output judge points
  const string jp_file="/home/wegatron/tmp/cube_jp01.obj";
  ofstream ofs(jp_file, std::ofstream::out) ;
  const vector<zsw::Tet>& tets = triangulation.getTets();
  for(const zsw::Tet& tet : tets) {
    for(const zsw::JudgePoint &jp : tet.jpts_) {
      ofs << "v " << jp.pt_[0] << " " << jp.pt_[1] << " " << jp.pt_[2] << std::endl;
    }
  }
  // output edges
  const vector<zsw::Edge>& edges = triangulation.getEdges();
  for(const zsw::Edge &e : edges) {
    std::cerr << "e: " << e.vid_[0] << " " << e.vid_[1] << std::endl;
  }
}

void test1()
{
  // input obj mesh
  const string file_path="/home/wegatron/workspace/geometry/data/cube.obj";
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
  // output triangulation
  triangulation.writeTetMesh("/home/wegatron/tmp/cube_bbox.vtk", zsw::BBOX_POINT);
}


int main(int argc, char *argv[])
{
  test0();
  test1();
  test01();
  return 0;
}

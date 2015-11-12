#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include "../triangulation2.h"
#include "../surface_generator.h"

using namespace std;

void test(const std::string &file_path, const string &output_prefix, const zsw::Scalar thick_dis)
{
  zsw::mesh::TriMesh input_mesh;
  if(!OpenMesh::IO::read_mesh(input_mesh, file_path)) {
    std::cerr << "[ERROR] can't read mesh!" << std::endl;
  }

  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bo_points;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bi_points;
  zsw::genPoints(thick_dis, input_mesh, bo_points, bi_points);

  // judge points is 10 times dense as the basic mesh points
  zsw::Triangulation tr(thick_dis/10, bo_points, bi_points);
  std::function<bool(const zsw::Tet&)> ignore_bbox
    = std::bind(&zsw::Triangulation::ignoreWithPtType, &tr, std::placeholders::_1, zsw::BBOX_POINT);
  std::function<bool(const zsw::Tet&)> ignore_self_out
    = std::bind(&zsw::Triangulation::ignoreOnlyWithPtType, &tr, std::placeholders::_1, zsw::OUTER_POINT);
  std::function<bool(const zsw::Tet&)> ignore_self_in
    = std::bind(&zsw::Triangulation::ignoreOnlyWithPtType, &tr, std::placeholders::_1, zsw::INNER_POINT);

  tr.writeTetMesh(output_prefix+"_ori.vtk", {ignore_bbox, ignore_self_out, ignore_self_in});
  tr.writeTetMesh(output_prefix+"_ori_with_bbox.vtk", {});
  tr.simpTolerance();
  //tr.mutualTessellation();
  tr.writeTetMesh(output_prefix+"_simp_tol.vtk", {ignore_bbox, ignore_self_out, ignore_self_in});
  tr.writeTetMesh(output_prefix+"_simp_tol_with_bbox.vtk", {});
}

int main(int argc, char *argv[])
{
  test("/home/wegatron/workspace/geometry/data/sphere.stl", "/home/wegatron/tmp/simp_tol/sphere", 0.1);
  return 0;
}

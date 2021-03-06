#include <fstream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include <zswlib/error_ctrl.h>
#include "../triangulation2.h"
#include "../surface_generator.h"

using namespace std;

void test(const std::string &file_path, const string &output_prefix, const zsw::Scalar thick_dis, const zsw::Scalar sample_r)
{
  zsw::mesh::TriMesh input_mesh;
  if(!OpenMesh::IO::read_mesh(input_mesh, file_path)) {
    std::cerr << "[ERROR] can't read mesh!" << std::endl;
  }

  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bo_points;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bi_points;
  zsw::genPoints(thick_dis, input_mesh, bo_points, bi_points);

  zsw::Triangulation tr0(0.3, output_prefix+"_tmp/");
  CALL_FUNC(tr0.construct(0.1, sample_r, 0.3, bo_points, bi_points), abort());

  tr0.saveStatus(output_prefix+"_tmp/tmp_status", zsw::BC_STAGE);

  zsw::Triangulation tr(0.3, output_prefix+"_tmp/");
  zsw::Status status;
  tr.loadStatus(output_prefix+"_tmp/tmp_status", status);

  std::function<bool(const zsw::Tet&)> ignore_bbox
    = std::bind(&zsw::Triangulation::ignoreWithPtType, &tr, std::placeholders::_1, zsw::BBOX_POINT);
  std::function<bool(const zsw::Tet&)> ignore_out
    = std::bind(&zsw::Triangulation::ignoreWithPtType, &tr, std::placeholders::_1, zsw::OUTER_POINT);
  std::function<bool(const zsw::Tet&)> ignore_self_out
    = std::bind(&zsw::Triangulation::ignoreOnlyWithPtType, &tr, std::placeholders::_1, zsw::OUTER_POINT);
  std::function<bool(const zsw::Tet&)> ignore_self_in
    = std::bind(&zsw::Triangulation::ignoreOnlyWithPtType, &tr, std::placeholders::_1, zsw::INNER_POINT);

  std::function<bool(const zsw::Tet&)> ignore_not_with_zero_point
    = std::bind(&zsw::Triangulation::ignoreNotWithPtType, &tr, std::placeholders::_1, zsw::ZERO_POINT);

  tr.writeTetMesh(output_prefix+"tol_ori.vtk", {ignore_bbox, ignore_self_out, ignore_self_in});
  tr.writeTetMesh(output_prefix+"tol_ori_all.vtk", {});
  tr.writeTetMesh(output_prefix+"tol_in_ori.vtk", {ignore_bbox, ignore_out});

  tr.simpTolerance();
  tr.writeTetMesh(output_prefix+"_simp_tol_before_mt.vtk", {ignore_bbox, ignore_self_out, ignore_self_in});
  // tr.mutualTessellation();
  // tr.writeTetMesh(output_prefix+"_simp_tol_after_mt.vtk", {ignore_not_with_zero_point});
  // tr.writeSurface(output_prefix+"_simp_tol_after_zero_mt.obj", zsw::ZERO_POINT);
  // tr.writeSurface(output_prefix+"_simp_tol_after_inner_mt.obj", zsw::INNER_POINT);
  // tr.writeSurface(output_prefix+"_simp_tol_after_outer_mt.obj", zsw::OUTER_POINT);
}

int main(int argc, char *argv[])
{
  test("/home/wegatron/workspace/geometry/data/sphere.stl",
       "/home/wegatron/tmp/approximate/sphere/sphere", 0.1, 0.03);
  return 0;
}

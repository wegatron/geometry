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
    abort();
  }

  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bo_points;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bi_points;
  zsw::genPoints(thick_dis, input_mesh, bo_points, bi_points);

  zsw::Triangulation tr;
  CALL_FUNC(tr.construct(sample_r, bo_points, bi_points), abort());

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
  tr.writeTetMesh(output_prefix+"tol_in_ori.vtk", {ignore_bbox, ignore_out});
  tr.writeTetMesh(output_prefix+"tol_total_ori.vtk", {});

  //tr.simpTolerance();
  tr.writeTetMesh(output_prefix+"_simp_tol_before_mt.vtk", {ignore_bbox, ignore_self_out, ignore_self_in});
  tr.writeTetMesh(output_prefix+"tol_total_before_mt.vtk", {});
  tr.writeTetMeshAdjVs("/home/wegatron/tmp/tet_adj80.vtk", {80});
  tr.writeTetMeshAdjVs("/home/wegatron/tmp/tet_adj124.vtk", {124});
  CALL_FUNC(tr.isGood(), abort());

  const std::vector<zsw::Tet> &tets=tr.getTets();
  size_t tcnt=0;
  for(size_t ti=0; ti<tets.size(); ++ti) {
    if(!tets[ti].valid_) { continue; }
    size_t cnt=0;
    for(size_t vid : tets[ti].vid_) {      if(vid==78 || vid==120) { ++cnt; }    }
    if(cnt==2) {      tr.writeTet("/home/wegatron/tmp/e78_120_"+std::to_string(tcnt++)+".vtk", ti);    }
  }

  tr.mutualTessellation();
  tr.writeTetMesh(output_prefix+"tol_total_after_mt.vtk", {});
  tr.writeTetMesh(output_prefix+"_simp_tol_after_mt.vtk", {ignore_bbox, ignore_self_out, ignore_self_in});
  tr.writeSurface2(output_prefix+"_zero_surf_before_simp.obj", zsw::ZERO_POINT);

  CALL_FUNC(tr.isGood(), abort());
  tr.simpZeroSurface();
  tr.writeTetMesh(output_prefix+"tol_total_last.vtk", {});
  tr.writeTetMesh(output_prefix+"_after_simp_zero.vtk", {ignore_bbox, ignore_self_out, ignore_self_in});
  tr.writeSurface2(output_prefix+"_simped_zero.obj", zsw::ZERO_POINT);

  CALL_FUNC(tr.isGood(), abort());
  std::cerr << "tris fine after simpZeroSurface" << std::endl;
}

int main(int argc, char *argv[])
{
  test("/home/wegatron/workspace/geometry/data/sphere.stl", "/home/wegatron/tmp/approximate/sphere/sphere", 0.1, 0.01);
  return 0;
}

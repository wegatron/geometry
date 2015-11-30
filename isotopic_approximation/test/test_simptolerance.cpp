#include <fstream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
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

  // judge points is 10 times dense as the basic mesh points
  zsw::Triangulation tr(sample_r, bo_points, bi_points);
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
  // std::cout << "input two vid:" << std::endl;
  // vector<size_t> vids(2,0);
  // std::cin >> vids[0] >> vids[1];
  // tr.writeTetMeshAdjVs("/home/wegatron/tmp/simp_tol/debug/adj.vtk", vids);
  // tr.writeBoundTris("/home/wegatron/tmp/simp_tol/debug/bound_tri.obj", vids[0], vids[1]);
  // while(1) {
  //   tr.testCollapseDebug(vids[0], vids[1]);
  // }
#if 0
  {
    std::ofstream ofs_in(output_prefix+"jp_in.obj", std::ofstream::out) ;
    std::ofstream ofs_out(output_prefix+"jp_out.obj", std::ofstream::out) ;
    const vector<zsw::Tet>& tets = tr.getTets();
    for(const zsw::Tet &tet : tets) {
      for(const zsw::JudgePoint &jp : tet.jpts_) {
        if(jp.val_exp_ > 0 ) { ofs_out << "v " << jp.pt_[0] << " " << jp.pt_[1] << " " << jp.pt_[2] << std::endl; }
        else { ofs_in << "v " << jp.pt_[0] << " " << jp.pt_[1] << " " << jp.pt_[2] << std::endl; }
      }
    }
  }
#endif

  tr.simpTolerance();

  tr.writeTetMesh(output_prefix+"_simp_tol_before_mt.vtk", {ignore_bbox, ignore_self_out, ignore_self_in});
  // tr.writeTetMeshAdjV(output_prefix+"_simp_tol_before_mt_adjv4", 4);
  tr.mutualTessellation();
  tr.writeTetMesh(output_prefix+"_simp_tol_after_mt.vtk", {ignore_not_with_zero_point});
  tr.writeSurface(output_prefix+"_simp_tol_after_zero_mt.obj", zsw::ZERO_POINT);
  tr.writeSurface(output_prefix+"_simp_tol_after_inner_mt.obj", zsw::INNER_POINT);
  tr.writeSurface(output_prefix+"_simp_tol_after_outer_mt.obj", zsw::OUTER_POINT);
}

int main(int argc, char *argv[])
{
  test(argv[1], argv[2], atof(argv[3]), atof(argv[4]));
  //test("/home/wegatron/workspace/geometry/data/sphere.stl", "/home/wegatron/tmp/simp_tol/sphere/sphere", 0.1, 0.01);
  //test("/home/wegatron/workspace/geometry/data/fandisk.obj", "/home/wegatron/tmp/simp_tol/fandisk/fandisk", 0.006, 0.003);
  return 0;
}

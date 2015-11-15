#include <fstream>
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
  std::function<bool(const zsw::Tet&)> ignore_out
    = std::bind(&zsw::Triangulation::ignoreWithPtType, &tr, std::placeholders::_1, zsw::OUTER_POINT);
  std::function<bool(const zsw::Tet&)> ignore_self_out
    = std::bind(&zsw::Triangulation::ignoreOnlyWithPtType, &tr, std::placeholders::_1, zsw::OUTER_POINT);
  std::function<bool(const zsw::Tet&)> ignore_self_in
    = std::bind(&zsw::Triangulation::ignoreOnlyWithPtType, &tr, std::placeholders::_1, zsw::INNER_POINT);

  std::function<bool(const zsw::Tet&)> ignore_not_with_zero_point
    = std::bind(&zsw::Triangulation::ignoreNotWithPtType, &tr, std::placeholders::_1, zsw::ZERO_POINT);

  // {
  //   const string out_jp_file=output_prefix+"out_jp.obj";
  //   std::ofstream ofs(out_jp_file, std::ofstream::out) ;
  //   const vector<zsw::Tet>& tets = tr.getTets();
  //   for(const zsw::Tet& tet : tets) {
  //     for(const zsw::JudgePoint &jp : tet.jpts_) {
  //       if(jp.val_exp_<0.5) { continue; }
  //       ofs << "v " << jp.pt_[0] << " " << jp.pt_[1] << " " << jp.pt_[2] << std::endl;
  //     }
  //   }
  // }

  tr.writeTetMesh(output_prefix+"tol_ori.vtk", {ignore_bbox, ignore_self_out, ignore_self_in});
  tr.writeTetMesh(output_prefix+"tol_in_ori.vtk", {ignore_bbox, ignore_out});
  //tr.writeTetMesh(output_prefix+"_tol_zero.vtk", {ignore_bbox, ignore_out});
  //tr.writeTetMesh(output_prefix+"_ori_with_bbox.vtk", {});
  tr.simpTolerance();

  {
    const string in_jp_file=output_prefix+"in_jp.obj";
    std::ofstream ofs(in_jp_file, std::ofstream::out) ;
    const vector<zsw::Tet>& tets = tr.getTets();
    for(const zsw::Tet &tet : tets) {
      for(const zsw::JudgePoint &jp : tet.jpts_) {
        ofs << "v " << jp.pt_[0] << " " << jp.pt_[1] << " " << jp.pt_[2] << std::endl;
      }
    }
    ofs.close();
  }

  // {
  //   size_t wtet_id=-1;
  //   const string in_jp_file=output_prefix+"in_jp.obj";
  //   std::ofstream ofs(in_jp_file, std::ofstream::out) ;
  //   const vector<zsw::Tet>& tets = tr.getTets();
  //   size_t jp_id=0;
  //   for( size_t t_id=0; t_id<tets.size(); ++t_id) {
  //     const zsw::Tet &tet=tets[t_id];
  //     if(!tet.valid_) { continue; }
  //     std::cerr << "tet[" << t_id << "] jpts size:" << tet.jpts_.size() << std::endl;
  //     for(const zsw::JudgePoint &jp : tet.jpts_) {
  //       if(jp.val_exp_>0.5) { continue; }
  //       if(fabs(jp.val_exp_-jp.val_cur_) >1.0) {
  //         std::cerr << "[ERROR] !!!!!!!!!!!!!!!!!!!!" << std::endl;
  //       }
  //       if(jp_id==70100) { std::cerr << "70100->tid:" << t_id << std::endl; wtet_id=t_id;  }
  //       ++jp_id;
  //       ofs << "v " << jp.pt_[0] << " " << jp.pt_[1] << " " << jp.pt_[2] << std::endl;
  //     }
  //   }
  //   std::cerr << "wtet_id:" << wtet_id << std::endl;
  //   tr.writeTet(output_prefix+"jpt_debug_tet.vtk", wtet_id);
  // }

  tr.writeTetMesh(output_prefix+"_simp_tol_before_mt.vtk", {ignore_bbox, ignore_self_out, ignore_self_in});
  // tr.writeTetMeshAdjV(output_prefix+"_simp_tol_before_mt_adjv4", 4);
  tr.mutualTessellation();
  //tr.writeTetMesh(output_prefix+"_simp_tol_after_mt.vtk", {ignore_not_with_zero_point});
  tr.writeSurface(output_prefix+"_simp_tol_after_mt.obj", zsw::ZERO_POINT);
  // tr.writeTetMeshAdjV(output_prefix+"_simp_tol_after_mt_adjv4", 4);
  //tr.writeTetMesh(output_prefix+"_simp_zero.vtk", {ignore_bbox, ignore_out});
  //tr.writeTetMesh(output_prefix+"_simp_tol_with_bbox.vtk", {});
}

int main(int argc, char *argv[])
{
  test("/home/wegatron/workspace/geometry/data/sphere.stl", "/home/wegatron/tmp/simp_tol/sphere", 0.1);
  return 0;
}

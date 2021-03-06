#include <fstream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include <zswlib/error_ctrl.h>
#include <zswlib/zsw_clock_c11.h>
#include "../triangulation2.h"
#include "../surface_generator.h"

using namespace std;

void test(const std::string &file_path, const string &output_prefix, const zsw::Scalar normal_cond_scale,
          const zsw::Scalar thick_dis, const zsw::Scalar sample_r, const zsw::Scalar flat_threshold)
{
  std::ofstream ofs;
  OPEN_STREAM(output_prefix+"_tmp/run_time", ofs, std::ofstream::out, abort());
  zsw::common::ClockC11 zsw_clock;
  zsw::mesh::TriMesh input_mesh;
  if(!OpenMesh::IO::read_mesh(input_mesh, file_path)) {
    std::cerr << "[ERROR] can't read mesh!" << std::endl;
    abort();
  }

  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bo_points;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bi_points;
  zsw::genPoints(thick_dis, input_mesh, bo_points, bi_points);

  zsw::Triangulation tr(normal_cond_scale, output_prefix+"_tmp/");
  CALL_FUNC(tr.construct(flat_threshold, sample_r, thick_dis, bo_points, bi_points), abort());

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
  tr.saveStatus(output_prefix+"_tmp/status_before_simp_tol", zsw::BC_STAGE);
  ofs << "construct triangulation time:" << zsw_clock.time() << std::endl;
  tr.simpTolerance();
  tr.writeTetMesh(output_prefix+"_simp_tol_before_mt.vtk", {ignore_bbox, ignore_self_out, ignore_self_in});
  tr.writeTetMesh(output_prefix+"tol_total_before_mt.vtk", {});
  assert(tr.isGood()==0);
  ofs<< "simp_tolerance_time_cost:" << zsw_clock.time() << std::endl;
  tr.mutualTessellation();
  tr.writeTetMesh(output_prefix+"tol_total_after_mt.vtk", {});
  tr.writeTetMesh(output_prefix+"_simp_tol_after_mt.vtk", {ignore_bbox, ignore_self_out, ignore_self_in});
  tr.writeSurface2(output_prefix+"_zero_surf_before_simp.obj", zsw::ZERO_POINT);
  assert(tr.isGood()==0);
  ofs<< "simp_mutuall_tessellation_time_cost:" << zsw_clock.time() << std::endl;

  tr.simpZeroSurface();
  tr.writeTetMesh(output_prefix+"tol_total_last.vtk", {});
  tr.writeTetMesh(output_prefix+"_after_simp_zero.vtk", {ignore_bbox, ignore_self_out, ignore_self_in});
  tr.writeSurface(output_prefix+"_simped_zero.obj", zsw::ZERO_POINT);
  ofs << "simp_zero_surface_time_cost:" << zsw_clock.time() << std::endl;
  ofs << "total_time_cost:" << zsw_clock.totalTime() << std::endl;
  assert(tr.isGood()==0);
}

int main(int argc, char *argv[])
{
  test(argv[1], argv[2], atof(argv[3]), atof(argv[4]), atof(argv[5]), atof(argv[6]));
  return 0;
}

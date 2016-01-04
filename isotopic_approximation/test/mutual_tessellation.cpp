#include <fstream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include <zswlib/error_ctrl.h>
#include "../triangulation2.h"
#include "../surface_generator.h"

using namespace std;

void test(const string &filepath, const string &output_prefix, const zsw::Scalar thick,
          const zsw::Scalar sample_r, const zsw::Scalar flat_threshold)
{
  zsw::mesh::TriMesh in_mesh;
  if(!OpenMesh::IO::read_mesh(in_mesh, filepath)) {
    std::cerr << "[ERROR] can't read mesh!" << std::endl;
    abort();
  }

  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bo_points;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bi_points;
  zsw::genPoints(thick, in_mesh, bo_points, bi_points);
  zsw::Triangulation tr(0.3,"/home/wegatron/tmp/");
  CALL_FUNC(tr.construct(flat_threshold, sample_r, 0.3, bo_points, bi_points), abort());
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
  tr.writeTetMesh(output_prefix+"debug_all.vtk", {});
  tr.mutualTessellation();
  tr.writeSurface2(output_prefix+"debug_mutual.obj", zsw::ZERO_POINT);
}

int main(int argc, char *argv[])
{
  if(atoi(argv[1]) & 1) {
    test("/home/wegatron/workspace/geometry/data/cylinder_smoothed.ply", "/home/wegatron/tmp/mutual_tessellation/cylinder/cylinder_"+std::string(argv[2])+"_", 0.02, 0.01, atof(argv[2]));
  }
  if(atoi(argv[1])&2) {
    test("/home/wegatron/workspace/geometry/data/fandisk_smoothed.ply", "/home/wegatron/tmp/mutual_tessellation/fandisk/fandisk_"+std::string(argv[2])+"_", 0.004, 0.02, atof(argv[2]));
  }
  if(atoi(argv[1])&4) {
    test("/home/wegatron/workspace/geometry/data/fertility.stl", "/home/wegatron/tmp/mutual_tessellation/fertility/fertility_"+std::string(argv[2])+"_", 0.5, 0.2, atof(argv[2]));
  }
  if(atoi(argv[1])&8) {
    test("/home/wegatron/workspace/geometry/data/bunny.obj", "/home/wegatron/tmp/mutual_tessellation/bunny/bunny_"+std::string(argv[2])+"_", 0.002, 0.0008, atof(argv[2]));
  }
  return 0;
}

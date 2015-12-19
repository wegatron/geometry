#include <fstream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include <zswlib/error_ctrl.h>
#include "../triangulation2.h"
#include "../surface_generator.h"

using namespace std;

void debugSimpJpts(const std::string &filepath, const std::string &output_prefix,
               const zsw::Scalar thick_dis, const zsw::Scalar sample_r,
               const size_t v0, const size_t v1)
{
  zsw::mesh::TriMesh input_mesh;
  if(!OpenMesh::IO::read_mesh(input_mesh, filepath)) {
    std::cerr << "[ERROR] can't read mesh!" << std::endl;
  }

  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bo_points;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bi_points;
  zsw::genPoints(thick_dis, input_mesh, bo_points, bi_points);

  // judge points is 10 times dense as the basic mesh points
  zsw::Triangulation tr(0.3, output_prefix+"_tmp/");
  CALL_FUNC(tr.construct(0.25, sample_r, bo_points, bi_points), abort());

  std::function<bool(const zsw::Tet&)> ignore_bbox
    = std::bind(&zsw::Triangulation::ignoreWithPtType, &tr, std::placeholders::_1, zsw::BBOX_POINT);
  std::function<bool(const zsw::Tet&)> ignore_out
    = std::bind(&zsw::Triangulation::ignoreWithPtType, &tr, std::placeholders::_1, zsw::OUTER_POINT);
  std::function<bool(const zsw::Tet&)> ignore_self_out
    = std::bind(&zsw::Triangulation::ignoreOnlyWithPtType, &tr, std::placeholders::_1, zsw::OUTER_POINT);
  std::function<bool(const zsw::Tet&)> ignore_self_in
    = std::bind(&zsw::Triangulation::ignoreOnlyWithPtType, &tr, std::placeholders::_1, zsw::INNER_POINT);
  tr.writeTetMesh(output_prefix+"tol.vtk", {ignore_bbox, ignore_self_out, ignore_self_in});
  tr.writeTetMesh(output_prefix+"tol_in.vtk", {ignore_bbox, ignore_out});

  {
    const vector<zsw::Vertex> &vertices = tr.getVertices();
    const vector<zsw::Edge> &edges = tr.getEdges();
    size_t target_eid=-1;
    for(size_t e_id : vertices[v0].edge_ids_) {
      if(edges[e_id].valid_ &&
         (edges[e_id].vid_[0]==v1 || edges[e_id].vid_[1]==v1)) {        target_eid=e_id; break;      }
    }
    std::set<size_t> eids_set;
    tr.tryCollapseBoundaryEdge(target_eid, eids_set);
  }
}

int main(int argc, char *argv[])
{
  debugSimpJpts(argv[1], argv[2], atof(argv[3]), atof(argv[4]), atoi(argv[5]), atoi(argv[6]));
  return 0;
}

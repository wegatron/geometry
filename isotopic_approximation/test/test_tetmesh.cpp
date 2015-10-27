#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include "../triangulation.h"
#include "../surface_generator.h"

using namespace std;


int main(int argc, char *argv[])
{
  zsw::mesh::TriMesh input_mesh;
  std::string file_path="/home/wegatron/workspace/geometry/data/sphere.stl";
  if(!OpenMesh::IO::read_mesh(input_mesh, file_path)) {
    std::cerr << "[ERROR] can't read mesh!" << std::endl;
  }
  zsw::SurfaceGenerator sfg;
  vector<zsw::Point> bz_points, bo_points, bi_points;
  sfg.genPoints(0.2, input_mesh, bz_points, bo_points, bi_points);
  zsw::TetMesh tm(bz_points, bo_points, bi_points, 0.02);

  tm.writeVtk("/home/wegatron/tmp/tmp.vtk");

  const vector<zsw::TetMesh::Tet> &tets = tm.getTets();
  const vector<zsw::TetMesh::Vertex> &vertices = tm.getVertices();
  std::cerr << "vertices_[34].tet_ids_.size()=" << vertices[34].tet_ids_.size() << std::endl;
  std::cerr << "==========================" << std::endl;
  size_t ccnt=0;
  for(int i=0; i<tets.size(); ++i) {
    if(!tets[i].valid_) { continue; }
    assert(!(vertices[tets[i].vind0_].pt_type_==vertices[tets[i].vind1_].pt_type_ &&
             vertices[tets[i].vind1_].pt_type_==vertices[tets[i].vind2_].pt_type_ &&
             vertices[tets[i].vind2_].pt_type_==vertices[tets[i].vind3_].pt_type_));
    if(tets[i].vind0_ !=34 && tets[i].vind1_ !=34 && tets[i].vind2_ !=34 && tets[i].vind3_ !=34) { continue; }
    ++ccnt;
  }
  std::cerr << "size expected:" << ccnt << std::endl;
  return 0;
}

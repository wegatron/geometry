#define BOOST_AUTO_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include "../triangulation.h"
#include "../surface_generator.h"

using namespace std;

BOOST_AUTO_TEST_SUITE(collapse_edge_with_fake_kernel_region)

void calcFv(const zsw::TetMesh::Edge &edge, const vector<zsw::TetMesh::Vertex> &vertices,
            const vector<zsw::TetMesh::Tet> &tets,
            set<size_t> &fv)
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
}

BOOST_AUTO_TEST_CASE(collapse_edge_fkr0)
{
  zsw::mesh::TriMesh input_mesh;
  std::string file_path="/home/wegatron/workspace/geometry/data/sphere.stl";
  if(!OpenMesh::IO::read_mesh(input_mesh, file_path)) {
    BOOST_FAIL("[ERROR] can't read mesh!");
  }
  zsw::SurfaceGenerator sfg;
  vector<zsw::Point> bz_points, bo_points, bi_points;
  sfg.genPoints(0.2, input_mesh, bz_points, bo_points, bi_points);
  zsw::TetMesh tm(bz_points, bo_points, bi_points, 0.02);

  const std::vector<zsw::TetMesh::Vertex> &vertices = tm.getVertices();
  const std::vector<zsw::TetMesh::Tet> &tets = tm.getTets();
  const std::vector<zsw::TetMesh::Edge> &edges = tm.getEdges();

  vector<size_t> invalid_tet_ids;
  vector<size_t> updated_tet_ids;
  const size_t vind0=67;
  const size_t vind1=62;
  for(size_t tet_id : vertices[vind1].tet_ids_) {
    if(tets[tet_id].vind0_==vind0 || tets[tet_id].vind1_==vind0
       || tets[tet_id].vind2_==vind0 || tets[tet_id].vind3_==vind0) {
      invalid_tet_ids.push_back(tet_id);
    } else {      updated_tet_ids.push_back(tet_id);    }
  }

  tm.writeVtk("/home/wegatron/tmp/tmp_before.vtk");
  const zsw::Point target_pt((vertices[vind1].pt_[0]+vertices[vind0].pt_[0])/2, (vertices[vind1].pt_[1]+vertices[vind0].pt_[1])/2, (vertices[vind1].pt_[2]+vertices[vind0].pt_[2])/2);

  BOOST_CHECK_MESSAGE(tm.testCollapseEdge(vind0,vind1),
                      "No edge "+ std::to_string(vind0) + "-" + std::to_string(vind1) +" to collapse!");

  const std::vector<zsw::TetMesh::Vertex> &vertices2 = tm.getVertices();
  // check vertex vind1's father is vind0'
  BOOST_CHECK_MESSAGE(vertices[vind1].father_==vind0, "vertices[" + std::to_string(vind1)+"].father_="
                      +std::to_string(vertices[vind1].father_)+" which should be "+std::to_string(vind0));
  // check vertex vind0's pt is the center of the collapsed edge
  BOOST_CHECK(vertices[vind0].pt_ == target_pt);
  // check tets
  for(size_t tet_id : invalid_tet_ids) {    BOOST_CHECK(!tets[tet_id].valid_);   }
  for(size_t tet_id : updated_tet_ids) {
    BOOST_CHECK_MESSAGE(
    tets[tet_id].vind0_==vind0 || tets[tet_id].vind1_==vind0 || tets[tet_id].vind2_==vind0 || tets[tet_id].vind3_==vind0,
                        std::to_string(tets[tet_id].vind0_)+","+std::to_string(tets[tet_id].vind0_)+","
    +std::to_string(tets[tet_id].vind0_)+","+std::to_string(tets[tet_id].vind0_));
  }
  // check edges' vind
  size_t tmp_cnt=0;
  for(size_t eid : vertices[vind1].edge_ids_) {
    if(edges[eid].vind0_==vind1 || edges[eid].vind1_==vind1) { BOOST_CHECK(!edges[eid].valid_); ++tmp_cnt; }
    else { BOOST_CHECK(edges[eid].vind0_==vind0 || edges[eid].vind1_==vind0); }
  }
  BOOST_CHECK(tmp_cnt==1);
  /// check edges' fv and fv_cnt
  // check each edge linked by adjacent vertex's fv is right

  set<size_t> vids;
  for(size_t tet_id : updated_tet_ids) {
    if(tets[tet_id].vind0_!=vind0) { vids.insert(tets[tet_id].vind0_); }
    if(tets[tet_id].vind1_!=vind0) { vids.insert(tets[tet_id].vind1_); }
    if(tets[tet_id].vind2_!=vind0) { vids.insert(tets[tet_id].vind2_); }
    if(tets[tet_id].vind3_!=vind0) { vids.insert(tets[tet_id].vind3_); }
  }

  for(size_t vid : vids) {
    for(size_t eid : vertices[vid].edge_ids_) {
      if(!edges[eid].valid_) { continue; }
      set<size_t> fv;
      calcFv(edges[eid], vertices, tets, fv);
      BOOST_CHECK(fv.size() == edges[eid].fv_cnt_);
      bool flag=true;
      size_t tfv_cnt=0;
      for(size_t tfv : edges[eid].fv_) {
        if(vertices[tfv].father_!=-1) { continue; }
        if(fv.find(tfv)==fv.end()) { flag=false; break; }
        ++tfv_cnt;
      }
      BOOST_CHECK(flag && tfv_cnt==edges[eid].fv_cnt_);
    }
  }

  tm.writeVtk("/home/wegatron/tmp/tmp_after.vtk");
}

BOOST_AUTO_TEST_SUITE_END()

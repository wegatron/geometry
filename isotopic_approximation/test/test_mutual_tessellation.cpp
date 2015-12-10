#define BOOST_AUTO_TEST_MAIN

#include <boost/test/included/unit_test.hpp>

#include <fstream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>
#include "../triangulation2.h"
#include "../surface_generator.h"

using namespace std;

BOOST_AUTO_TEST_SUITE(test_tessellation)

bool pairComp(const std::pair<size_t, size_t> &a, const std::pair<size_t, size_t> &b)
{
  if(a.first < b.first) { return true; }
  if(a.first==b.first && a.second<b.second) { return true; }
  return false;
}

BOOST_AUTO_TEST_CASE(test_sphere)
{
  zsw::mesh::TriMesh in_mesh;
  if(!OpenMesh::IO::read_mesh(in_mesh, "/home/wegatron/workspace/geometry/data/sphere.obj")) {
    std::cerr << "[ERROR] can't read mesh!" << std::endl;
    abort();
  }

  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bo_points;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bi_points;
  zsw::genPoints(0.1, in_mesh, bo_points, bi_points);
  zsw::Triangulation tr(0.05, bo_points, bi_points);
  if(!tr.isGoodTriangulation()) { std::cerr << "invalid triangulation!" << std::endl; abort(); }

  {
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
    tr.writeTetMesh("/home/wegatron/tmp/debug_tol_ori.vtk", {ignore_bbox, ignore_self_out, ignore_self_in});
    tr.writeTetMesh("/home/wegatron/tmp/debug_all.vtk", {});
  }

  tr.mutualTessellation();
  const std::vector<zsw::Tet> &tets = tr.getTets();
  const std::vector<zsw::Edge> &edges = tr.getEdges();
  const std::vector<zsw::Vertex> &vertices = tr.getVertices();
  std::set<std::pair<size_t,size_t>,zsw::PairCompFunc> e_set_tr(pairComp), e_set_check(pairComp);
  size_t ev_cnt=0;
  for(const zsw::Edge &e : edges) {
    if(!e.valid_) { continue; }
    ++ev_cnt;
    std::pair<size_t,size_t> ep;
    ep.first=e.vid_[0]; ep.second=e.vid_[1];
    if(ep.first>ep.second) { std::swap(ep.first, ep.second); }
    e_set_tr.insert(ep);
  }
  BOOST_CHECK(e_set_tr.size()==ev_cnt);

  for(const zsw::Tet &tet : tets) {
    if(!tet.valid_) { continue; }
    size_t vids[4]; std::copy(tet.vid_, tet.vid_+4, vids); sort(vids, vids+4);
    e_set_check.insert(std::pair<size_t,size_t>(vids[0], vids[1]));
    e_set_check.insert(std::pair<size_t,size_t>(vids[0], vids[2]));
    e_set_check.insert(std::pair<size_t,size_t>(vids[0], vids[3]));
    e_set_check.insert(std::pair<size_t,size_t>(vids[1], vids[2]));
    e_set_check.insert(std::pair<size_t,size_t>(vids[1], vids[3]));
    e_set_check.insert(std::pair<size_t,size_t>(vids[2], vids[3]));
  }

  auto it_tr = e_set_tr.begin();
  auto it_check = e_set_check.begin();
  while(it_tr!=e_set_tr.end() && it_check!=e_set_check.end()) {
    if(it_tr->first != it_check->first || it_tr->second!=it_check->second) {
      std::cerr << "it_tr:" << it_tr->first << " " << it_tr->second << std::endl;
      std::cerr << "it_check:" << it_check->first << " " << it_check->second << std::endl;
      std::cerr << "it_tr_type:" << vertices[it_tr->first].pt_type_ << " " << vertices[it_tr->second].pt_type_ << std::endl;
      std::cerr << "it_check_type:" << vertices[it_check->first].pt_type_ << " " << vertices[it_check->second].pt_type_ << std::endl;
      std::abort();
    }
    ++it_tr; ++it_check;
  }
  BOOST_CHECK(e_set_check.size()==ev_cnt);
}

BOOST_AUTO_TEST_SUITE_END()

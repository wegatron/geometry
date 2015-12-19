#define BOOST_AUTO_TEST_MAIN

#include <fstream>
#include <boost/test/included/unit_test.hpp>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include<zswlib/mesh/mesh_type.h>
#include <zswlib/mesh/vtk.h>
#include <zswlib/error_ctrl.h>
#include "../surface_generator.h"

#include "../triangulation2.h"

using namespace std;

BOOST_AUTO_TEST_SUITE(link_condition)

BOOST_AUTO_TEST_CASE(fake_triangulation)
{
  zsw::Triangulation tr(0.3,"/home/wegatron/tmp/");
  std::vector<zsw::Vertex> &vertices=tr.getVertices();
  std::vector<zsw::Tet> &tets=tr.getTets();
  std::vector<zsw::Edge> &edges=tr.getEdges();
  Eigen::Matrix<zsw::Scalar,3,1> v[5];
  v[0]<<0,0,0; v[1]<<1,0,0; v[2]<<0,1,1; v[3]<<0.5,0.5,1; v[4]=(v[1]+v[2])/2.0;
  vertices.push_back({true, zsw::OUTER_POINT, v[0], {0,1}, {}});  vertices.push_back({true, zsw::OUTER_POINT, v[1], {0,2},{}});
  vertices.push_back({true, zsw::OUTER_POINT, v[2], {1,2}, {}});  vertices.push_back({true, zsw::OUTER_POINT, v[3], {0,1,2}, {}});
  vertices.push_back({true, zsw::OUTER_POINT, v[4], {0,1,2}, {}});

  tets.push_back({true,{0,1,3,4}}); tets.push_back({true,{0,2,3,4}}); tets.push_back({true,{1,2,3,4}});
  edges.push_back({true, 0,1}); edges.push_back({true, 0,2});
  edges.push_back({true, 0,3}); edges.push_back({true, 0,4});
  edges.push_back({true, 1,2}); edges.push_back({true, 1,3});
  edges.push_back({true, 1,4}); edges.push_back({true, 2,3});
  edges.push_back({true, 2,4}); edges.push_back({true, 3,4});

  BOOST_CHECK(!tr.testLinkCondition(edges[0]));
  BOOST_CHECK(!tr.testLinkCondition(edges[1]));
  BOOST_CHECK(tr.testLinkCondition(edges[2]));
  BOOST_CHECK(tr.testLinkCondition(edges[3]));
  BOOST_CHECK(!tr.testLinkCondition(edges[4]));
  BOOST_CHECK(tr.testLinkCondition(edges[5]));
  BOOST_CHECK(tr.testLinkCondition(edges[6]));
  BOOST_CHECK(tr.testLinkCondition(edges[7]));
  BOOST_CHECK(tr.testLinkCondition(edges[8]));
  BOOST_CHECK(tr.testLinkCondition(edges[9]));
}

BOOST_AUTO_TEST_CASE(real_triangulation)
{
  const string file_path="/home/wegatron/workspace/geometry/data/cube_2.obj";
  zsw::mesh::TriMesh input_mesh;
  if(!OpenMesh::IO::read_mesh(input_mesh, file_path)) {
    std::cerr << "[ERROR] can't read mesh!" << std::endl;
  }

  // gen bi and bo points
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bo_points;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bi_points;
  zsw::genPoints(0.5, input_mesh, bo_points, bi_points);
  // init triangulation
  zsw::Triangulation tr(0.3,"/home/wegatron/tmp/");
  CALL_FUNC(tr.construct(0.25, 0.1, bo_points, bi_points), abort());
  const std::vector<zsw::Edge> &edges=tr.getEdges();
  const std::vector<zsw::Vertex> &vertices=tr.getVertices();
  const std::vector<zsw::Tet> &tets=tr.getTets();
  for(const zsw::Edge &e : edges) {
    if(!tr.testLinkCondition(e)) {
      std::cerr << "link condition fail with edge : " << e.vid_[0] << ", " << e.vid_[1] << std::endl;
      size_t i=0;
      for(size_t tid : vertices[e.vid_[0]].tet_ids_) {
        tr.writeTet("/home/wegatron/tmp/test_linkcondition/cube_linkcond_single"
                               +std::to_string(i++)+".vtk", tid);
      }
      std::set<size_t> fv; // vertex construct a face with edge e
      std::set<size_t> adj_v0; // vertex link e.vid_[0]
      std::set<size_t> adj_v1; // vertex link e.vid_[1]
      for(size_t tid : vertices[e.vid_[0]].tet_ids_) {
        bool isfv=false;
        for(size_t vid : tets[tid].vid_) {
          if(vid == e.vid_[1]) { isfv=true; }
          adj_v0.insert(vid);
        }
        if(isfv) {    for(size_t vid : tets[tid].vid_) {      fv.insert(vid);    }    }
      }
      fv.erase(e.vid_[0]); fv.erase(e.vid_[1]);
      std::cerr << "fv : ";
      for(size_t fv_id : fv) {
        std::cerr << fv_id << " ";
      }
      std::cerr<< std::endl;
      adj_v0.erase(e.vid_[0]); adj_v0.erase(e.vid_[1]);

      for(size_t tid : vertices[e.vid_[1]].tet_ids_) {
        for(size_t vid : tets[tid].vid_) {
          adj_v1.insert(vid);
        }
      }
      adj_v1.erase(e.vid_[0]); adj_v1.erase(e.vid_[1]);

      std::vector<size_t> cv;
      std::set<size_t>::iterator it0=adj_v0.begin();
      std::set<size_t>::iterator it1=adj_v1.begin();
      while(it0!=adj_v0.end() && it1!=adj_v1.end()) {
        if(*it0 == *it1) { cv.push_back(*it0); ++it0; ++it1; }
        else if(*it0>*it1) { ++it1; }
        else { ++it0; }
      }
      std::cerr << "cv : ";
      for(size_t cv_id : cv) {
        std::cerr << cv_id << " ";
      }
      std::cerr<< std::endl;
      break;
    }
  }
  tr.writeTetMesh("/home/wegatron/tmp/test_linkcondition/cube_link_cond.vtk", {});
}

BOOST_AUTO_TEST_SUITE_END()

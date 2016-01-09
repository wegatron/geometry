#include <iostream>
#include <fstream>
#include <unordered_set>
#include <zswlib/error_ctrl.h>
#include "isotopic_approximation.h"
#include "basic_op.h"
#include "bound_sphere.h"
#include "sampling.h"

using namespace std;

namespace zsw{

  void Approximation::init(const zsw::Scalar &surf_sample_r, const zsw::Scalar &tet_sample_r,
                           const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bi_vertices,
                           const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bo_vertices)
  {
    surf_sample_r_=surf_sample_r;
    tet_sample_r_=tet_sample_r;
    Eigen::Matrix<zsw::Scalar,3,2> bbox;
    calcBBOX(bo_vertices, bbox);
    zsw::Scalar scale = 0.5*(bbox.block<3,1>(0,0) - bbox.block<3,1>(0,1)).norm();
    Eigen::Matrix<zsw::Scalar,3,1> transform = 0.5*(bbox.block<3,1>(0,0)+bbox.block<3,1>(0,1));
    zsw::BoundSphere bs("/home/wegatron/workspace/geometry/data/bound_sphere.obj", scale, transform);
    const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bs_vertices = bs.getVertices();
    std::vector<std::pair<Point, VertexInfo>> vertices;
    size_t vid=0;
    for(const Eigen::Matrix<zsw::Scalar,3,1> &v : bs_vertices) {
      vertices.push_back({Point(v[0], v[1], v[2]), VertexInfo(vid++, zsw::BBOX_POINT, v, 0.0)});
    }
    for(const Eigen::Matrix<zsw::Scalar,3,1> &v : bi_vertices) {
      vertices.push_back({Point(v[0],v[1],v[2]), VertexInfo(vid++, zsw::INNER_POINT, v, 0.0)});
    }
    for(const Eigen::Matrix<zsw::Scalar,3,1> &v : bo_vertices) {
      vertices.push_back({Point(v[0],v[1],v[2]), VertexInfo(vid++, zsw::OUTER_POINT, v, 0.0)});
    }
    tw_.reset(new zsw::TriangulationWapper(vertices));
    createJudgePoints();
  }

  void Approximation::createJudgePoints()
  {
    const TTds &tds=tw_->getTds();
    Eigen::Matrix<zsw::Scalar,3,4> bi_tri_points;
    Eigen::Matrix<zsw::Scalar,3,4> bo_tri_points;
    for(TTds::Cell_iterator cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) {
      size_t bi_cnt=0;
      size_t bo_cnt=0;
      for(size_t i=0; i<4; ++i) {
        if(cit->vertex(i)->info().pt_type_==zsw::INNER_POINT) {
          bi_tri_points(0,bi_cnt)=cit->vertex(i)->point()[0];
          bi_tri_points(1,bi_cnt)=cit->vertex(i)->point()[1];
          bi_tri_points(2,bi_cnt)=cit->vertex(i)->point()[2];
          ++bi_cnt;
        } else if(cit->vertex(i)->info().pt_type_==zsw::OUTER_POINT) {
          bi_tri_points(0,bo_cnt)=cit->vertex(i)->point()[0];
          bi_tri_points(1,bo_cnt)=cit->vertex(i)->point()[1];
          bi_tri_points(2,bo_cnt)=cit->vertex(i)->point()[2];
          ++bo_cnt;
        }
      }
      if(bi_cnt==3 && bo_cnt==1) { sampleTriangle(bi_tri_points.block<3,3>(0,0), surf_sample_r_, bi_jpts_);}
      else if(bo_cnt==3 && bi_cnt==1) { sampleTriangle(bo_tri_points.block<3,3>(0,0), surf_sample_r_, bo_jpts_); }
    }
    jpts_ptr_bi_.reset(new zsw::Flann<zsw::Scalar>(bi_jpts_[0].data(), bi_jpts_.size()));
    jpts_ptr_bo_.reset(new zsw::Flann<zsw::Scalar>(bo_jpts_[0].data(), bo_jpts_.size()));
    jpts_.reserve(bi_jpts_.size()+bo_jpts_.size());
    for(const Eigen::Matrix<zsw::Scalar,3,1> &jpt : bi_jpts_) { jpts_.push_back({jpt,-1,-1}); }
    for(const Eigen::Matrix<zsw::Scalar,3,1> &jpt : bo_jpts_) { jpts_.push_back({jpt,1,1}); }
  }

  void Approximation::simpTolerance()
  {
    TTds &tds=tw_->getTds();
    std::unordered_map<std::string, TTds::Edge> edge_map;
    for(TTds::Edge_iterator eit=tds.edges_begin();
        eit!=tds.edges_end(); ++eit) {
      if(tw_->isBoundaryEdge(*eit)) {
        std::pair<size_t,size_t> key(eit->first->vertex(eit->second)->info().index_,
                                     eit->first->vertex(eit->second)->info().index_);
        if(key.first>key.second) { swap(key.first, key.second); }
        std::string key_str=std::to_string(key.first) + "," + std::to_string(key.second);
        if(edge_map.find(key_str)==edge_map.end()) { edge_map.insert(std::make_pair(key_str, *eit)); }
#if 0
        else {
          CGAL::Container_from_circulator<TTds::Cell_circulator> cells0(tds.incident_cells(*eit));
          std::cerr << "cell0:" << std::endl;
          for(auto cell : cells0) {
            std::cerr << cell.vertex(0)->info().index_ << " " <<
              cell.vertex(1)->info().index_ << " " <<
              cell.vertex(2)->info().index_ << " " <<
              cell.vertex(3)->info().index_ << std::endl;
          }

          CGAL::Container_from_circulator<TTds::Cell_circulator> cells1(tds.incident_cells(edge_map[key_str]));
          std::cerr << "cell1:" << std::endl;
          for(auto cell : cells1) {
            std::cerr << cell.vertex(0)->info().index_ << " " <<
              cell.vertex(1)->info().index_ << " " <<
              cell.vertex(2)->info().index_ << " " <<
              cell.vertex(3)->info().index_ << std::endl;
          }
          std::cerr << "---------------------" << std::endl;
        }
#endif
      }
    }
    while(!edge_map.empty()) {
      TTds::Edge &e = edge_map.begin()->second; edge_map.erase(edge_map.begin());
      if(tds.is_edge(e.first, e.second, e.third)) { continue; }
      tryCollapseBoundaryEdge(e, edge_map);
    }
  }

  void Approximation::writeJudgePoints(const std::string &filepath) const
  {
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
        << jpts_.size() << " float" << std::endl;
    for(const JudgePoint &jpt : jpts_) {
      ofs << jpt.pt_.transpose() << std::endl;
    }
    ofs << "CELLS " << jpts_.size() << " " << jpts_.size()*2 << std::endl;
    for(size_t i=0; i<jpts_.size(); ++i) {      ofs << "1 " << i << std::endl;    }
    ofs << "CELL_TYPES " << jpts_.size() << std::endl;
    for(size_t i=0; i<jpts_.size(); ++i) { ofs << "1" << std::endl; }
  }

  void Approximation::tryCollapseBoundaryEdge(TTds::Edge &e,
                                              std::unordered_map<std::string,TTds::Edge> &edge_map)
  {
    if(!tw_->isSatisfyLinkCondition(e)) { return; }
    std::vector<Fhd> bound_tris;
    std::vector<Vhd> opposite_vs;
    tw_->calcBoundTris(e, bound_tris, opposite_vs);
    Eigen::Matrix<zsw::Scalar,3,2> bbox;
    calcFhdBBox(bound_tris, bbox);
    std::vector<JudgePoint*> jpts_in_bbox;
    for(JudgePoint &jpt : jpts_) {
      if(jpt.pt_[0]<bbox(0,0) || jpt.pt_[1]<bbox(1,0) || jpt.pt_[2]<bbox(2,0) ||
         jpt.pt_[0]>bbox(0,1) || jpt.pt_[1]>bbox(1,1) || jpt.pt_[2]>bbox(2,1)) { continue; }
      jpts_in_bbox.push_back(&jpt);
    }
    // candicate merge points in kernel region
    KernelRegionJudger krj;
    constructKernelRegionJudger(bound_tris, opposite_vs, krj);
    std::vector<JudgePoint*> candicate_point;
    for(JudgePoint * jpt_ptr : jpts_in_bbox) {
      if(krj.judge(jpt_ptr->pt_)) { candicate_point.push_back(jpt_ptr); }
    }
    // sort jpt by error
    sort(candicate_point.begin(), candicate_point.end(), [](const JudgePoint *a, const JudgePoint *b){
        return fabs(a->val_cur_-a->val_exp_) > fabs(b->val_cur_-b->val_exp_);
      });
    JudgePoint *merge_pt=nullptr;
    std::vector<VertexUpdateData> vup;
    std::vector<JudgePointUpdateData> jup;
    for(JudgePoint * jpt_ptr : candicate_point) {
      vup.clear(); jup.clear();
      if(isSatisfyErrorBound(bound_tris, jpts_in_bbox, jpt_ptr->pt_, vup, &jup)) { merge_pt=jpt_ptr; break; }
    }
    Vhd vhd=e.first->vertex(0);
    tw_->collapseEdge(e, vhd, merge_pt->pt_);
    updateVertex(vup);
    updateJudgePoint(jup);
    boundaryEdgeBack(vhd, edge_map);
  }

  void Approximation::tryCollapseZeroEdge(TTds::Edge &e,
                                          std::unordered_map<std::string,TTds::Edge> &edge_map)
  {
    TTds &tds=tw_->getTds();
    if(!tw_->isSatisfyLinkCondition(e)) { return; }
    std::vector<Fhd> bound_tris;
    std::vector<Vhd> opposite_vs;
    tw_->calcBoundTris(e, bound_tris, opposite_vs);
    std::vector<JudgePoint*> jpts_in_bbox;
    calcJptsInBbox(bound_tris, jpts_in_bbox);
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> sample_points;
    sampleIncidentCells(e, sample_points);
    KernelRegionJudger krj;
    constructKernelRegionJudger(bound_tris, opposite_vs, krj);
    const Eigen::Matrix<zsw::Scalar,3,1> *merge_pt=nullptr;
    std::vector<VertexUpdateData> vup;
    for(const Eigen::Matrix<zsw::Scalar,3,1> &pt : sample_points) {
      vup.clear(); // can't parallel
      if(krj.judge(pt) && isSatisfyErrorBound(bound_tris, jpts_in_bbox, pt, vup, nullptr)) { merge_pt=&pt; break; }
    }
    if(merge_pt==nullptr) { return; }
    Vhd vhd=e.first->vertex(0);
    tw_->collapseEdge(e, vhd, *merge_pt);
    updateVertex(vup);
    zeroEdgeBack(vhd, edge_map);
  }

  void Approximation::mutuallTessellation()
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  }

  void Approximation::simpZeroSurface()
  {
    TTds &tds = tw_->getTds();
    std::unordered_map<std::string, TTds::Edge> edge_map;
    for(TTds::Edge_iterator eit=tds.edges_begin();
        eit!=tds.edges_end(); ++eit) {
      if(eit->first->vertex(0)->info().pt_type_!=zsw::ZERO_POINT
         || eit->first->vertex(1)->info().pt_type_!=zsw::ZERO_POINT) { continue; }
      size_t key[2] = {eit->first->vertex(0)->info().index_, eit->first->vertex(1)->info().index_};
      if(key[0]>key[1]) { swap(key[0],key[1]); }
      std::string key_str(std::to_string(key[0])+","+std::to_string(key[1]));
      if(edge_map.find(key_str)==edge_map.end()) { edge_map.insert(std::make_pair(key_str, *eit)); }
    }
    while(!edge_map.empty()) {
      TTds::Edge &e=edge_map.begin()->second; edge_map.erase(edge_map.begin());
      if(tds.is_edge(e.first, e.second, e.third)) { continue; }
      tryCollapseZeroEdge(e, edge_map);
    }
  }

  void Approximation::writeZeroSurface(const std::string &filepath) const
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  }

  void Approximation::writeTetMesh(const std::string &filepath,
                                   std::vector<std::function<bool(const TTds::Cell_handle)>> ignore_tet_funcs) const
  {
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    const DelaunayTriangulation &delaunay=tw_->getDelaunay();
    const TTds &tds=tw_->getTds();
    ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
        << tds.number_of_vertices() << " float" << std::endl;

    CGAL::Unique_hash_map<TTds::Vertex_handle,size_t> v_map;
    size_t v_index=0;
    for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
      ofs << *vit << std::endl;
      v_map[vit] = v_index++;
    }

    stringstream ss;
    size_t valid_cells_number=0;
    for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) {
      bool ignore=false;
      for(auto igf : ignore_tet_funcs) {        if(igf(cit)) { ignore=true; break; }      }
      if(ignore) { continue; }
      ss << "4 " << v_map[cit->vertex(0)] << " " <<
        v_map[cit->vertex(1)] << " " <<
        v_map[cit->vertex(2)] << " " <<
        v_map[cit->vertex(3)] << std::endl;
      ++valid_cells_number;
    }
    ofs << "CELLS "<< valid_cells_number << " " << valid_cells_number*5 <<std::endl;
    ofs << ss.str();
    ofs << "CELL_TYPES " << valid_cells_number << std::endl;
    for(size_t i=0; i<valid_cells_number; ++i) {      ofs << "10" << std::endl;    }
    ofs.close();
  }

  void calcFhdBBox(const std::vector<Fhd> &bound_tris, Eigen::Matrix<zsw::Scalar,3,2> &bbox)
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  }

}
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_set>
#include <zswlib/error_ctrl.h>
#include "isotopic_approximation.h"
#include "basic_op.h"
#include "bound_sphere.h"
#include "sampling.h"
#include "smoother.h"

using namespace std;

namespace zsw{

  void Approximation::simpZeroSurface()
  {
    std::unordered_map<std::string, TTds::Edge> z_map;
    const TTds &tds = tw_->getTds();
    for(TTds::Edge_iterator eit=tds.edges_begin();
        eit!=tds.edges_end(); ++eit) {
      if(!tw_->isZeroEdge(*eit)) { continue; }
      std::string key_str=edge2key(*eit);
      z_map.insert(std::make_pair(key_str, *eit));
    }
    size_t z_c_step=0;
    size_t try_z_c_step=0;
    while(!z_map.empty()) {
      TTds::Edge e=z_map.begin()->second; z_map.erase(z_map.begin());
      if(!tw_->isZeroEdge(e)) { continue; }
      if(tryCollapseBZEdge(e, z_map)) {
        if(++z_c_step%50==0) { std::cout << "[INFO] zero edge collapsed " << z_c_step << std::endl; }
      }
      if(++try_z_c_step%100==0) { std::cout << "[INFO] try zero edge collapsed " << try_z_c_step << std::endl; }
    }
    NZSWLOG("zsw_info") << "zero edge try collapse:" << try_z_c_step << std::endl;
    NZSWLOG("zsw_info") << "zero edge collapsed total:" << z_c_step << std::endl;
    NZSWLOG("zsw_info") << "zero edge collapse suc:" << z_c_step*1.0/try_z_c_step << std::endl;
  }

  bool Approximation::tryCollapseBZEdge(TTds::Edge &e,
                                        std::unordered_map<std::string,TTds::Edge> &z_map,
                                        bool is_bz_back)
  {
    if(!tw_->isSatisfyLinkCondition(e)) { return false; }
    std::vector<VertexTriple> bound_tris;
    std::vector<Vhd> opposite_vs;
    tw_->calcBoundTris(e, bound_tris, opposite_vs);
    std::vector<const JudgePoint*> jpts_in_bbox;
    calcJptsInBbox(&bound_tris[0].first, 3*bound_tris.size(), jpts_in_bbox);
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> sample_points;
    sampleAdjCells(e, sample_points);
    KernelRegionJudger krj;
    constructKernelRegionJudger(bound_tris, opposite_vs, krj);
#if 1
    std::vector<std::vector<VertexUpdateData>> tmp_vups;
    const Eigen::Matrix<zsw::Scalar,3,1> *merge_pt=nullptr;
#pragma omp parallel for
    for(size_t i=0; i<sample_points.size(); ++i) {
      if(merge_pt!=nullptr) { continue; }
      const Eigen::Matrix<zsw::Scalar,3,1> &pt=sample_points[i];
      if(krj.judge(pt) &&
         isTetsSatisfyNormalCondition(bound_tris, pt, zsw::ZERO_POINT)
         && isSatisfyErrorBound(bound_tris, jpts_in_bbox, pt, 0, tmp_vups[i], nullptr)) {
        merge_pt=&sample_points[i];
      }
    }
    if(merge_pt==nullptr) { return false; }
    std::vector<VertexUpdateData> &vup=tmp_vups[merge_pt-&sample_points[0]];
#else
    const Eigen::Matrix<zsw::Scalar,3,1> *merge_pt=nullptr;
    std::vector<VertexUpdateData> vup;
    for(const Eigen::Matrix<zsw::Scalar,3,1> &pt : sample_points) {
      vup.clear(); // can't parallel
      if(krj.judge(pt) &&
         isTetsSatisfyNormalCondition(bound_tris, pt, zsw::ZERO_POINT)
         && isSatisfyErrorBound(bound_tris, jpts_in_bbox, pt, 0, vup, nullptr)) { merge_pt=&pt; break; }
    }
    if(merge_pt==nullptr) { return false; }
#endif
    //Vhd vhd=e.first->vertex(e.second);
    Vhd vhd=(e.first->vertex(e.second)->info().pt_type_==zsw::ZERO_POINT) ? e.first->vertex(e.second) : e.first->vertex(e.third);
    tw_->collapseEdge(e, vhd, *merge_pt);
    //updateVertex(vup);
    if(is_bz_back) { bzEdgeBack(vhd, z_map); }
    else { zeroEdgeBack(vhd, z_map); }
    return true;
  }

  void Approximation::zeroEdgeBack(Vhd vhd, std::unordered_map<std::string,TTds::Edge> &edge_map) const
  {
    const TTds &tds=tw_->getTds();
    std::vector<TTds::Edge> edges;
    tds.incident_edges(vhd, std::back_inserter(edges));
    for(const TTds::Edge &e : edges) {
      if(!tw_->isZeroEdge(e)) { continue; }
      std::string key_str=edge2key(e);
      edge_map.insert(std::make_pair(key_str, e));
    }
  }
}

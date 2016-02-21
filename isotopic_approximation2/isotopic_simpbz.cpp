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

  void Approximation::simpBZEdges()
  {
    std::unordered_map<std::string, TTds::Edge> bz_map;
    const TTds &tds=tw_->getTds();
    for(auto eit=tds.edges_begin(); eit!=tds.edges_end(); ++eit) {
      if(!tw_->isBZEdge(*eit)) { continue; }
      std::string key_str=edge2key(*eit);
      bz_map.insert(std::make_pair(key_str, *eit));
    }

    size_t zb_c_step=0;
    while(!bz_map.empty()) {
      TTds::Edge e=bz_map.begin()->second; bz_map.erase(bz_map.begin());
      if(!tw_->isBZEdge(e) && !tw_->isZeroEdge(e)) { continue; }
      if(tryCollapseBZEdge(e, bz_map)) {
        std::cout << "[INFO] all edge collapsed " << zb_c_step << std::endl;
        ++zb_c_step;
      }
    }
    std::cout << "[INFO] zb collapsed total:" << zb_c_step << std::endl;
  }

  bool Approximation::tryCollapseBZEdge(TTds::Edge &e, std::unordered_map<std::string,TTds::Edge> &bz_map)
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
    const Eigen::Matrix<zsw::Scalar,3,1> *merge_pt=nullptr;
    std::vector<VertexUpdateData> vup;
    for(const Eigen::Matrix<zsw::Scalar,3,1> &pt : sample_points) {
      vup.clear(); // can't parallel
      if(krj.judge(pt) && isSatisfyErrorBound(bound_tris, jpts_in_bbox, pt, 0, vup, nullptr)) { merge_pt=&pt; break; }
    }
    if(merge_pt==nullptr) { return false; }
    Vhd vhd=(e.first->vertex(e.second)->info().pt_type_==zsw::ZERO_POINT) ? e.first->vertex(e.second) : e.first->vertex(e.third);
    //Vhd vhd=e.first->vertex(e.second); vhd->info().pt_type_==zsw::ZERO_POINT;
#if 0
    static size_t cnt=0;
    ++cnt;
    if(cnt==58) {
    writeAdjcentCells("/home/wegatron/tmp/adj_cell.vtk", e);
    writeJudgePoints("/home/wegatron/tmp/jpts_in_bbox.vtk", jpts_in_bbox);
  }
#endif
    tw_->collapseEdge(e, vhd, *merge_pt);
    //updateVertex(vup);
    bzEdgeBack(vhd, bz_map);
    // if(z_map!=nullptr) { zeroEdgeBack(vhd, *z_map); }
    return true;
  }

  void Approximation::bzEdgeBack(Vhd vhd, std::unordered_map<std::string, TTds::Edge> &edge_map) const
  {
    const TTds &tds=tw_->getTds();
    std::vector<TTds::Edge> edges;
    tds.incident_edges(vhd, std::back_inserter(edges));
    for(const TTds::Edge &e : edges) {
      if(!tw_->isBZEdge(e) || !tw_->isZeroEdge(e)) { continue; }
      std::string key_str=edge2key(e);
      edge_map.insert(std::make_pair(key_str, e));
    }
  }

}

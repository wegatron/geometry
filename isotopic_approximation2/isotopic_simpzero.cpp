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
    std::queue<std::pair<TTds::Edge, size_t>> z_q;
    const TTds &tds = tw_->getTds();
    for(TTds::Edge_iterator eit=tds.edges_begin();
        eit!=tds.edges_end(); ++eit) {
      if(tw_->isZeroEdge(*eit)) { z_q.push(std::make_pair(*eit, 0)); }
    }
    tw_->resetVertexLastUpdate();
    size_t zero_cnt = countZeroPoints();
    size_t z_step = 0;
    size_t z_step_suc = 0;
    //print_sp_size_= true;
    while(!z_q.empty()) {
      TTds::Edge e = z_q.front().first;
      size_t last_update = z_q.front().second;
      z_q.pop();
      if(!tw_->isZeroEdge(e)) { continue; }
      if(e.first->vertex(e.second)->info().last_update_ > last_update ||
         e.first->vertex(e.third)->info().last_update_ > last_update) { continue; }
      if(++z_step % 100 == 0) { std::cout << "[INFO] try zero edge collapsed " << z_step << std::endl; /* print_sp_size_ = true; */}
      if(tryCollapseBZEdge(e, z_q, z_step_suc, false)) {
        --zero_cnt;
        if(++z_step_suc %50 == 0) {
          std::cout << "zero edge collapsed " << z_step_suc << std::endl;
          NZSWLOG("zsw_info") << "Zero count " << zero_cnt << std::endl;
          NZSWLOG("time&pt_count") << tclock_.time() << "  " << zero_cnt << std::endl;
          writeTetMesh(tmp_outdir_+"simp_z"+to_string(z_step_suc)+".vtk", {zsw::ignore_out, zsw::ignore_bbox});
        }
      }
    }
    NZSWLOG("zsw_info") << "zero edge try collapse:" << z_step << std::endl;
    NZSWLOG("zsw_info") << "zero edge collapsed total:" << z_step_suc << std::endl;
    NZSWLOG("zsw_info") << "zero edge collapse suc:" << z_step_suc*1.0/z_step << std::endl;
  }

  bool Approximation::tryCollapseBZEdge(TTds::Edge &e,
                                        std::queue<std::pair<TTds::Edge, size_t>> &z_q,
                                        size_t cur_update,
                                        bool is_bz_back)
  {
    if(!tw_->isSatisfyLinkCondition(e)) { return false; }
    std::vector<VertexTriple> bound_tris;
    std::vector<Vhd> opposite_vs;
    tw_->calcBoundTris(e, bound_tris, opposite_vs);
    std::vector<const JudgePoint*> jpts_in_bbox;
    calcJptsInBbox(&bound_tris[0].first, 3*bound_tris.size(), jpts_in_bbox, true);
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> sample_points;
    sampleAdjCells(e, sample_points);
    KernelRegionJudger krj;
    constructKernelRegionJudger(bound_tris, opposite_vs, krj);
    const Eigen::Matrix<zsw::Scalar,3,1> *merge_pt=nullptr;
    bz_krj_need_judge_cnt_+=sample_points.size();
    std::vector<bool> in_krj(sample_points.size(), false);
    {
      BLOCK_TIME_ANALYSIS("krj_z");
#pragma omp parallel for
      for(size_t i=0; i<sample_points.size(); ++i) {
        in_krj[i] = krj.judge(sample_points[i]);
      }
    }
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> krj_points;
    size_t nj_pt = 0;
    for(size_t i=0; i<sample_points.size(); ++i) {      if(in_krj[i]) { krj_points.push_back(sample_points[i]); ++nj_pt; }    }
    size_t erj_pt = 0;
    for(const Eigen::Matrix<zsw::Scalar,3,1> &pt : krj_points) {
      if(isTetsSatisfyNormalCondition(bound_tris, pt, zsw::ZERO_POINT)) {
        ++erj_pt;
        bz_judge_pt_cnt_+=jpts_in_bbox.size();
        if(isSatisfyErrorBound(bound_tris, jpts_in_bbox, pt, 0, nullptr)) { merge_pt=&pt; break; }
      }
    }
    bz_normal_cond_judge_cnt_ += nj_pt;
    bz_error_bound_judge_cnt_ += erj_pt;
    // if(print_sp_size_) {
    //   print_sp_size_ = false;
    //   NZSWLOG("zsw_info") << "all pt=" << sample_points.size() << " nj_pt=" << nj_pt << " erj_pt=" << erj_pt << std::endl;
    //   PRINT_COST("krj_z");
    // }
    if(merge_pt==nullptr) { return false; }
    Vhd vhd=(e.first->vertex(e.second)->info().pt_type_==zsw::ZERO_POINT) ? e.first->vertex(e.second) : e.first->vertex(e.third);
    tw_->collapseEdge(e, vhd, *merge_pt);
    if(is_bz_back) { bzEdgeBack(vhd, z_q, cur_update); }
    else { zeroEdgeBack(vhd, z_q, cur_update); }
    return true;
  }

  void Approximation::zeroEdgeBack(Vhd vhd, std::queue<std::pair<TTds::Edge, size_t>> &z_q, size_t cur_update) const
  {
    const TTds &tds=tw_->getTds();
    std::vector<TTds::Edge> edges;
    tds.incident_edges(vhd, std::back_inserter(edges));
    for(const TTds::Edge &e : edges) {
      if(tw_->isZeroEdge(e)) { z_q.push(std::make_pair(e, cur_update)); }
    }
  }
}

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
  void Approximation::simpTolerance()
  {
    const TTds &tds=tw_->getTds();
    std::unordered_map<std::string, TTds::Edge> edge_map;
    for(TTds::Edge_iterator eit=tds.edges_begin();
        eit!=tds.edges_end(); ++eit) {
      if(!tw_->isBoundaryEdge(*eit)) { continue;  }
      std::string key_str=edge2key(*eit);
      edge_map.insert(std::make_pair(key_str, *eit));
    }
    NZSWLOG("zsw_info") << "edge size:" << edge_map.size() << std::endl;
    size_t b_c_step=0;
    size_t try_b_c_step=0;
    while(!edge_map.empty()) {
      TTds::Edge e = edge_map.begin()->second; edge_map.erase(edge_map.begin());
      if(!tw_->isBoundaryEdge(e)) { continue; }
      if(tryCollapseBoundaryEdge(e, edge_map)) {
        if(++b_c_step%50==0) {
          NZSWLOG("zsw_info") << "boundary collapsed:" << b_c_step << std::endl;
          writeTetMesh(tmp_outdir_+"sim_tol_"+std::to_string(b_c_step/50)
                       +".vtk", {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
        }
      }
      if(++try_b_c_step%100==0) {
        NZSWLOG("zsw_info") << "try boundary collapsed:" << try_b_c_step << std::endl;
      }
    }
    NZSWLOG("zsw_info") << "boundary collapsed total:" << b_c_step << std::endl;
    NZSWLOG("zsw_info") << "boundary collapsed suc:" << b_c_step*1.0/try_b_c_step << std::endl;
#if 0
    size_t normal_cond_debug_i=0;
    for(auto chd=tds.cells_begin(); chd!=tds.cells_end(); ++chd) {
      if(tw_->isTolCell(chd)) {
        Eigen::Matrix<zsw::Scalar,3,4> tri_pts;
        tri_pts<<
          chd->vertex(0)->point()[0], chd->vertex(1)->point()[0], chd->vertex(2)->point()[0], chd->vertex(3)->point()[0],
          chd->vertex(0)->point()[1], chd->vertex(1)->point()[1], chd->vertex(2)->point()[1], chd->vertex(3)->point()[1],
          chd->vertex(0)->point()[2], chd->vertex(1)->point()[2], chd->vertex(2)->point()[2], chd->vertex(3)->point()[2];
        Eigen::Matrix<zsw::Scalar,4,1> val;
        for(size_t i=0; i<4; ++i) {
          if(chd->vertex(i)->info().pt_type_==zsw::INNER_POINT){ val[i]=-1; }
          else { val[i]=1; }
        }
        Eigen::Matrix<zsw::Scalar,3,1> bc=0.25*(
                                                tri_pts.block<3,1>(0,0)+tri_pts.block<3,1>(0,1)+
                                                tri_pts.block<3,1>(0,2)+tri_pts.block<3,1>(0,3));
        const static Eigen::Matrix<zsw::Scalar,1,4> tmp_v=Eigen::Matrix<zsw::Scalar,1,4>::Ones()*(1-normal_cond_scale_);
        Eigen::Matrix<zsw::Scalar,3,4> scaled_tri_pts=normal_cond_scale_*tri_pts+bc*tmp_v;
        std::string filepath(tmp_outdir_+"normal_cond/nbc_"+std::to_string(normal_cond_debug_i++)+".vtk");
        if(!normalCondition(val, scaled_tri_pts, tri_pts, inner_jpts_, outer_jpts_,
                            inner_kdtree_, outer_kdtree_,
                            true, &filepath)) {
          std::cerr << "SimpTolerance result normal cond failed!!!" << std::endl;
          abort();
        }
      }
    }
#endif
  }

  bool Approximation::tryCollapseBoundaryEdge(TTds::Edge &e,
                                              std::unordered_map<std::string,TTds::Edge> &edge_map)
  {
    const TTds &tds = tw_->getTds();
    if(!tw_->isSatisfyLinkCondition(e)) {      return false;    }
    std::vector<VertexTriple> bound_tris;
    std::vector<Vhd> opposite_vs;
    tw_->calcBoundTris(e, bound_tris, opposite_vs);
    std::vector<const JudgePoint*> jpts_in_bbox;
    calcJptsInBbox(&bound_tris[0].first, 3*bound_tris.size(), jpts_in_bbox);
    // candicate merge points in kernel region
    KernelRegionJudger krj;
    constructKernelRegionJudger(bound_tris, opposite_vs, krj);
    std::function<bool(zsw::Scalar v)> judge_func;
    if(e.first->vertex(e.second)->info().pt_type_==zsw::INNER_POINT) {
      judge_func=[](zsw::Scalar v){ return v<0; };
    } else { judge_func=[](zsw::Scalar v){ return v>0; }; }
    std::vector<Plane> adj_zero_support_planes;
    tw_->calcAdjZeroSupportPlanes(e, adj_zero_support_planes);
    std::vector<bool> is_candicate(jpts_in_bbox.size(), false);
    std::vector<std::pair<const JudgePoint*, zsw::Scalar>> tmp_points(jpts_in_bbox.size(), std::make_pair(nullptr,0.0));
#pragma omp parallel for
    for(size_t i=0; i<jpts_in_bbox.size(); ++i) {
      tmp_points[i].first=jpts_in_bbox[i];
      tmp_points[i].second=0.0;
      if(!judge_func(jpts_in_bbox[i]->val_exp_)) { continue; }
      if(krj.judge(jpts_in_bbox[i]->pt_)) {
        is_candicate[i]=true;
        for(const Plane &plane : adj_zero_support_planes) {
          if(plane.normal_.dot(jpts_in_bbox[i]->pt_ - plane.v0_) <0) { is_candicate[i]=false; break; }
          zsw::Scalar tmp=plane.normal_.dot(jpts_in_bbox[i]->pt_ - plane.v0_);
          tmp_points[i].second+=tmp*tmp;
        }
      }
    }
    std::vector<std::pair<const JudgePoint*, zsw::Scalar>> candicate_points;
    for(size_t i=0; i<is_candicate.size(); ++i) { if(is_candicate[i]) { candicate_points.push_back(tmp_points[i]); } }
#if 0
    // test adj info and jpts selection is right
    writeAdjcentCells("/home/wegatron/tmp/adj_cell.vtk", e);
    writeJudgePoints("/home/wegatron/tmp/jpts_in_bbox.vtk", jpts_in_bbox);
    writeJudgePoints("/home/wegatron/tmp/candicate_points.vtk", candicate_points);
    abort();
#endif

    // sort jpt by error
    sort(candicate_points.begin(), candicate_points.end(), [](const std::pair<const JudgePoint*, zsw::Scalar> &a,
                                                              const std::pair<const JudgePoint*, zsw::Scalar> &b){
           return a.second<b.second;
         });
    const JudgePoint *merge_pt=nullptr;
    std::vector<VertexUpdateData> vup;
    std::vector<JudgePointUpdateData> jup;
    const zsw::Scalar v_pt=(e.first->vertex(e.second)->info().pt_type_==zsw::INNER_POINT) ? -1 : 1;
    for(std::pair<const JudgePoint*, zsw::Scalar> &cd_pt : candicate_points) {
      const JudgePoint *jpt_ptr=cd_pt.first;
      vup.clear(); jup.clear();
      if(isTetsSatisfyNormalCondition(bound_tris, jpt_ptr->pt_, e.first->vertex(e.second)->info().pt_type_)
         && isSatisfyErrorBound(bound_tris, jpts_in_bbox, jpt_ptr->pt_, v_pt, vup, &jup))
        { merge_pt=jpt_ptr; break; }
    }
    if(merge_pt==nullptr) { return false; }
    Vhd vhd=e.first->vertex(e.second);
    tw_->collapseEdge(e, vhd, merge_pt->pt_);
    std::for_each(jup.begin(), jup.end(), [](const JudgePointUpdateData &dt){dt.jpt->val_cur_=dt.val_cur_;});
    //updateVertex(vup);
    boundaryEdgeBack(vhd, edge_map);
#if 0
    if(!checkNormalCondition()) {
      std::cerr << "merge_pt:" << merge_pt->pt_.transpose() << std::endl;
      std::cerr << "bound_tris:" << std::endl;
      for(VertexTriple &vt : bound_tris) {
        std::cerr << vt.first->point()[0] << " " << vt.first->point()[1] << " " << vt.first->point()[2] << std::endl;
        std::cerr << vt.second->point()[0] << " " << vt.second->point()[1] << " " << vt.second->point()[2] << std::endl;
        std::cerr << vt.third->point()[0] << " " << vt.third->point()[1] << " " << vt.third->point()[2] << std::endl;
      }
      abort();
    }
#endif
    return true;
  }

  void Approximation::boundaryEdgeBack(Vhd vhd, std::unordered_map<std::string, TTds::Edge> &edge_map) const
  {
    const TTds &tds=tw_->getTds();
    std::vector<TTds::Edge> edges;
    tds.incident_edges(vhd, std::back_inserter(edges));
    for(const TTds::Edge &e : edges) {
      if(!tw_->isBoundaryEdge(e)) { continue; }
      std::string key_str=edge2key(e);
      edge_map.insert(std::make_pair(key_str, e));
    }
  }

  bool Approximation::isTetsSatisfyNormalCondition(const std::vector<VertexTriple> &bound_tris,
                                                   const Eigen::Matrix<zsw::Scalar,3,1> &pt,
                                                   const PointType point_type) const
  {
    assert(point_type==zsw::INNER_POINT || point_type==zsw::OUTER_POINT || point_type==zsw::ZERO_POINT);
    Eigen::Matrix<zsw::Scalar,4,1> val;
    Eigen::Matrix<zsw::Scalar,3,4> tri_pts;
    tri_pts.block<3,1>(0,0)=pt;
    const static Eigen::Matrix<zsw::Scalar,1,4> tmp_v=Eigen::Matrix<zsw::Scalar,1,4>::Ones()*(1-normal_cond_scale_);
    for(VertexTriple vt : bound_tris) {
      if(!isConstructTZCell(point_type, vt.first->info().pt_type_, vt.second->info().pt_type_, vt.third->info().pt_type_, val))
        { continue; }
      tri_pts(0,1)=vt.first->point()[0]; tri_pts(1,1)=vt.first->point()[1]; tri_pts(2,1)=vt.first->point()[2];
      tri_pts(0,2)=vt.second->point()[0]; tri_pts(1,2)=vt.second->point()[1]; tri_pts(2,2)=vt.second->point()[2];
      tri_pts(0,3)=vt.third->point()[0]; tri_pts(1,3)=vt.third->point()[1]; tri_pts(2,3)=vt.third->point()[2];
      Eigen::Matrix<zsw::Scalar,3,1> bc=0.25*(
                                              tri_pts.block<3,1>(0,0)+tri_pts.block<3,1>(0,1)+
                                              tri_pts.block<3,1>(0,2)+tri_pts.block<3,1>(0,3));
      Eigen::Matrix<zsw::Scalar,3,4> scaled_tri_pts=normal_cond_scale_*tri_pts+bc*tmp_v;
      if(!normalCondition(val, scaled_tri_pts, tri_pts, inner_jpts_, outer_jpts_, inner_kdtree_, outer_kdtree_)) { return false;}
    }
    return true;
  }
}

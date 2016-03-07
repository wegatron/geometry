#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_set>
#include <zswlib/error_ctrl.h>
#include "isotopic_approximation.h"
#include "isotopic_debug.h"
#include "basic_op.h"
#include "bound_sphere.h"
#include "sampling.h"
#include "smoother.h"
#include <zswlib/zsw_clock_c11.h>

using namespace std;

namespace zsw{
  void Approximation::refine(const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bs_jpts)
  {
#if 0
    testKdtree();
#endif
    for(const Eigen::Matrix<zsw::Scalar,3,1> &in_jpt : inner_jpts_) { jpts_.push_back({in_jpt, -1, 1}); }
    for(const Eigen::Matrix<zsw::Scalar,3,1> &out_jpt : outer_jpts_) { jpts_.push_back({out_jpt, 1, 1}); }
    inner_kdtree_.buildTree(inner_jpts_[0].data(), inner_jpts_.size());
    outer_kdtree_.buildTree(outer_jpts_[0].data(), outer_jpts_.size());
    std::vector<std::pair<Point, VertexInfo>> init_vertices;
    init_vertices.reserve(bs_jpts.size());
    for(size_t ind=0; ind<bs_jpts.size(); ++ind) {
      const Eigen::Matrix<zsw::Scalar,3,1> &pt=bs_jpts[ind];
      init_vertices.push_back(std::make_pair(Point(pt[0],pt[1],pt[2]),VertexInfo(ind, zsw::BBOX_POINT, pt, 0)));
    }
    tw_.reset(new TriangulationWapper(init_vertices));
    std::priority_queue<std::pair<zsw::Scalar,JudgePoint*>,
                        std::vector<std::pair<zsw::Scalar,JudgePoint*>>,
                        ErrorMaxComparison> err_queue;
    for(size_t i=0; i<jpts_.size();++i) {
      zsw::Scalar err=fabs(jpts_[i].val_cur_-jpts_[i].val_exp_);
      if(err>1) { err_queue.push(std::make_pair(err, &jpts_[i])); }
    }
    TTds &tds=tw_->getTds();
    size_t add_pt_for_err=0;
    size_t add_pt_for_normal=0;
    do{
      while(!err_queue.empty()) {
        std::pair<zsw::Scalar,JudgePoint*> jpt_info=err_queue.top(); err_queue.pop();
        zsw::Scalar real_err=fabs(jpt_info.second->val_cur_-jpt_info.second->val_exp_);
        if(fabs(real_err-jpt_info.first) > zsw::const_val::eps) { continue; } // have already updated
        PointType pt_type= (jpt_info.second->val_exp_<0) ? zsw::INNER_POINT : zsw::OUTER_POINT;
        VertexInfo vertex_info(-1, pt_type, jpt_info.second->pt_, 0.0);
        std::vector<Chd> chds;
        tw_->addPointInDelaunaySafe(jpt_info.second->pt_, vertex_info, chds);
        if(++add_pt_for_err%100==0) { NZSWLOG("zsw_info") << "add_pt_for_err:" << add_pt_for_err << std::endl;  }
        updateJptsInCells(chds, err_queue);
      }
      // find a tet viloate the normal condition
      bool viloate_normal_cond=false;
      std::vector<Chd> tmp_chds;
      for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) {
        if(!cit->info().satisfy_normal_cond_ && !checkUpNormalCondition(cit, tmp_chds)) {
          if(++add_pt_for_normal%50==0) { NZSWLOG("zsw_info") << "add_pt_for_normal:" << add_pt_for_normal << std::endl; }
          viloate_normal_cond=true;
          break;
        }
      }
      if(viloate_normal_cond) { updateJptsInCells(tmp_chds, err_queue); }
      else { break; }
    } while(true);
#if 0
    checkNormalcondition();
#endif
    NZSWLOG("zsw_info") << "refined add pt cnt:" << tw_->getTds().number_of_vertices()-bs_jpts.size()-1 << std::endl;
  }

  void Approximation::upgradeRefine(std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bs_jpts)
  {
    // scale inner_jpts_, outer_jpts_ bs_jpts in x direction by 0.25
    for(size_t i=0; i<inner_jpts_.size(); ++i) { inner_jpts_[i][0]*=0.25; }
    for(size_t i=0; i<outer_jpts_.size(); ++i) { outer_jpts_[i][0]*=0.25; }
    for(size_t i=0; i<bs_jpts.size(); ++i) { bs_jpts[i][0]*=0.25; }
    // deformed refine
    refine(bs_jpts);
    //writeTetMesh("/home/wegatron/tmp/deformed_refine_res.vtk", {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
    // scale back jpts_, inner_jpts_, outer_jpts_ in x direction
    for(size_t i=0; i<inner_jpts_.size(); ++i) { inner_jpts_[i][0]=inner_jpts_[i][0]*4; }
    for(size_t i=0; i<outer_jpts_.size(); ++i) { outer_jpts_[i][0]=outer_jpts_[i][0]*4; }
    // get deform back the tetmesh
    TTds &tds=tw_->getTds();
    for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
      if(vit->info().pt_type_==zsw::INVALID_POINT) { continue; }
      vit->set_point(Point(vit->point()[0]*4, vit->point()[1], vit->point()[2]));
    }
  }

  void Approximation::refine2(std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &ori_bs_jpts,
                              std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_inner_jpts,
                              std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_outer_jpts,
                              std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_bs_jpts)
  {
    for(const Eigen::Matrix<zsw::Scalar,3,1> &in_jpt : inner_jpts_) { jpts_.push_back({in_jpt, -1, 1}); }
    for(const Eigen::Matrix<zsw::Scalar,3,1> &out_jpt : outer_jpts_) { jpts_.push_back({out_jpt, 1, 1}); }
    inner_kdtree_.buildTree(inner_jpts_[0].data(), inner_jpts_.size());
    outer_kdtree_.buildTree(outer_jpts_[0].data(), outer_jpts_.size());
    std::vector<std::pair<Point, VertexInfo>> init_vertices;
    init_vertices.reserve(ori_bs_jpts.size());
    for(size_t ind=0; ind<ori_bs_jpts.size(); ++ind) {
      const Eigen::Matrix<zsw::Scalar,3,1> &pt=deformed_bs_jpts[ind];
      init_vertices.push_back(std::make_pair(Point(pt[0],pt[1],pt[2]),VertexInfo(ind, zsw::BBOX_POINT, ori_bs_jpts[ind], 0)));
    }
    tw_.reset(new TriangulationWapper(init_vertices));
    std::priority_queue<std::pair<zsw::Scalar,JudgePoint*>,
                        std::vector<std::pair<zsw::Scalar,JudgePoint*>>,
                        ErrorMaxComparison> err_queue;
    for(size_t i=0; i<jpts_.size();++i) {
      zsw::Scalar err=fabs(jpts_[i].val_cur_-jpts_[i].val_exp_);
      if(err>1) { err_queue.push(std::make_pair(err, &jpts_[i])); }
    }
    TTds &tds=tw_->getTds();
    size_t add_pt_for_err=0;
    while(!err_queue.empty()) {
      std::pair<zsw::Scalar,JudgePoint*> jpt_info=err_queue.top(); err_queue.pop();
      zsw::Scalar real_err=fabs(jpt_info.second->val_cur_-jpt_info.second->val_exp_);
      if(fabs(real_err-jpt_info.first) > zsw::const_val::eps) { continue; } // have already updated
      PointType pt_type= (jpt_info.second->val_exp_<0) ? zsw::INNER_POINT : zsw::OUTER_POINT;
      size_t deformed_jpt_index = std::distance(&jpts_[0], jpt_info.second);
      Eigen::Matrix<zsw::Scalar,3,1> deformed_pt = (deformed_jpt_index<inner_jpts_.size()) ? deformed_inner_jpts[deformed_jpt_index] : deformed_outer_jpts[deformed_jpt_index-inner_jpts_.size()];
      VertexInfo vertex_info(-1, pt_type, jpt_info.second->pt_, 0.0);
      std::vector<Chd> chds;
      tw_->addPointInDelaunaySafe(deformed_pt, vertex_info, chds);
      if(++add_pt_for_err%100==0) { NZSWLOG("zsw_info") << "add_pt_for_err:" << add_pt_for_err << std::endl; }
      updateJptsInCells2(chds, err_queue);
    }
    //std::cout << __FILE__ << __LINE__ << std::endl;
    // // find a tet viloate the normal condition
    // bool viloate_normal_cond=false;
    // std::vector<Chd> tmp_chds;
    // for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) {
    //   if(!cit->info().satisfy_normal_cond_ && !checkUpNormalCondition(cit, tmp_chds)) {
    //     if(++add_pt_for_normal%50==0) { NZSWLOG("zsw_info") << "add_pt_for_normal:" << add_pt_for_normal << std::endl; }
    //     viloate_normal_cond=true;
    //     break;
    //   }
    // }
    // if(viloate_normal_cond) { updateJptsInCells(tmp_chds, err_queue); }
    // else { break; }
    // } while(true);
    writeTetMesh(tmp_outdir_ +"before_remove_bbox_vertex.vtk", {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
    //tw_->removeBBoxPts();
    writeTetMesh(tmp_outdir_ +"before_swap_vertex.vtk", {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
    for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
      if(vit->info().pt_type_ != zsw::INVALID_POINT) {
        vit->set_point(Point(vit->info().pos_ori_[0], vit->info().pos_ori_[1], vit->info().pos_ori_[2]));
      }
    }
    NZSWLOG("zsw_info") << "refined add pt cnt:" << tw_->getTds().number_of_vertices()-ori_bs_jpts.size()-1 << std::endl;
  }

  bool Approximation::checkUpNormalCondition(Chd chd, std::vector<Chd> &chds)
  {
    if(!tw_->isTolCell(chd)) { return true; }
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
    if(normalCondition(val, scaled_tri_pts, tri_pts, inner_jpts_, outer_jpts_, inner_kdtree_, outer_kdtree_)) {
      chd->info().satisfy_normal_cond_=true;
      return true;
    }
#if 0
    static size_t ind=0;
    std::string filepath=tmp_outdir_+"normal_cond_debug/nc_"+std::to_string(ind)+".vtk";
    if(!normalCondition(val, scaled_tri_pts, tri_pts, inner_jpts_, outer_jpts_, inner_kdtree_, outer_kdtree_,
                        true, &filepath)) {
      std::cerr << "Refine result normal condition failed on " << ind << std::endl;
    }
    writeTetMesh(tmp_outdir_+"before_add_"+std::to_string(ind)+".vtk",  {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
    ind++;
#endif
    // insert point into delaunay triangulation
    Eigen::Matrix<zsw::Scalar,3,1> center;
    calcCircumcenter(tri_pts, center);
    // std::cerr << "center " << ind-1 << ":" << center.transpose() << std::endl;
    // find the minest
    std::vector<size_t> in_indices;
    std::vector<zsw::Scalar> in_dist;
    inner_kdtree_.queryNearest(center, in_indices, in_dist);
    std::vector<size_t> out_indices;
    std::vector<zsw::Scalar> out_dist;
    outer_kdtree_.queryNearest(center, out_indices, out_dist);
    size_t jpt_ind = (in_dist[0]<out_dist[0]) ? in_indices[0] : out_indices[0]+inner_jpts_.size();

    // std::cerr << "nearest pt:" << jpts_[jpt_ind].pt_ << std::endl;
    // if(ind-1!=25) {    return true; }
    // add and update
    PointType pt_type=(jpts_[jpt_ind].val_exp_<0) ? zsw::INNER_POINT : zsw::OUTER_POINT;
    VertexInfo vertex_info(-1, pt_type, jpts_[jpt_ind].pt_, 0.0);
    //writeTetMesh("/home/wegatron/tmp/before_add_pt.vtk",  {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
    tw_->addPointInDelaunaySafe(jpts_[jpt_ind].pt_, vertex_info, chds);
    //writeTetMesh("/home/wegatron/tmp/after_add_pt.vtk",  {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
    return false;
  }

  void calcCircumcenter(const Eigen::Matrix<zsw::Scalar,3,4> &tri_pts, Eigen::Matrix<zsw::Scalar,3,1> &center)
  {
    Eigen::Matrix<zsw::Scalar,3,1> n0,n1,n2;
    Eigen::Matrix<zsw::Scalar,3,3> A;
    n0=A.block<3,1>(0,0)=tri_pts.block<3,1>(0,1)-tri_pts.block<3,1>(0,0);
    n1=A.block<3,1>(0,1)=tri_pts.block<3,1>(0,2)-tri_pts.block<3,1>(0,0);
    n2=A.block<3,1>(0,2)=tri_pts.block<3,1>(0,3)-tri_pts.block<3,1>(0,0);
    Eigen::Matrix<zsw::Scalar,3,1> b;
    b[0] = 0.5*(A.block<3,1>(0,0)).dot(tri_pts.block<3,1>(0,1)+tri_pts.block<3,1>(0,0));
    b[1] = 0.5*(A.block<3,1>(0,1)).dot(tri_pts.block<3,1>(0,2)+tri_pts.block<3,1>(0,0));
    b[2] = 0.5*(A.block<3,1>(0,2)).dot(tri_pts.block<3,1>(0,3)+tri_pts.block<3,1>(0,0));
    Eigen::Matrix<zsw::Scalar,3,3> Atr = A.transpose();
    Eigen::PartialPivLU<Eigen::Matrix<zsw::Scalar,3,3>> pplu; pplu.compute(Atr);
    center = pplu.solve(b);
#if 0
    zsw::Scalar dis[4];
    for(size_t i=0; i<4; ++i) {
      dis[i]=(center-tri_pts.block<3,1>(0,i)).norm();
    }
    for(size_t i=1; i<4; ++i) {
      if(fabs(dis[i]-dis[i-1]) > zsw::const_val::eps) {
        std::cerr << Atr << std::endl;
        std::cerr << "---------------" << std::endl;
        std::cerr << n0.transpose() << std::endl;
        std::cerr << n1.transpose() << std::endl;
        std::cerr << n2.transpose() << std::endl;
        std::cerr << "err:" << (A*center-b).norm() << std::endl;
        std::cerr << fabs(dis[i]-dis[i-1]) << std::endl;
        std::cerr << "--" << n0.dot(center)-b[0]<< std::endl;
        std::cerr << "--" << n1.dot(center)-b[1] << std::endl;
        std::cerr << "--" << n2.dot(center)-b[2] << std::endl;
        abort();
      }
    }
#endif
  }
}

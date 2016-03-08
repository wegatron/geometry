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
    for(const Eigen::Matrix<zsw::Scalar,3,1> &in_jpt : inner_jpts_) { jpts_.push_back({in_jpt, in_jpt, -1, 1}); }
    for(const Eigen::Matrix<zsw::Scalar,3,1> &out_jpt : outer_jpts_) { jpts_.push_back({out_jpt, out_jpt, 1, 1}); }
    inner_kdtree_.buildTree(inner_jpts_[0].data(), inner_jpts_.size());
    outer_kdtree_.buildTree(outer_jpts_[0].data(), outer_jpts_.size());
    std::vector<std::pair<Point, VertexInfo>> init_vertices;
    init_vertices.reserve(bs_jpts.size());
    for(size_t ind=0; ind<bs_jpts.size(); ++ind) {
      const Eigen::Matrix<zsw::Scalar,3,1> &pt=bs_jpts[ind];
      init_vertices.push_back(std::make_pair(Point(pt[0],pt[1],pt[2]),VertexInfo(ind, zsw::BBOX_POINT, pt)));
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
        VertexInfo vertex_info(-1, pt_type, jpt_info.second->pt_cur_);
        std::vector<Chd> chds;
        tw_->addPointInDelaunay(jpt_info.second->pt_cur_, vertex_info, chds);
        if(++add_pt_for_err%100==0) { NZSWLOG("zsw_info") << "add_pt_for_err:" << add_pt_for_err << std::endl;  }
        updateJptsInCells(chds, err_queue);
      }
      // find a tet viloate the normal condition
      bool viloate_normal_cond=false;
      std::vector<Chd> tmp_chds;
      for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) {
        if(!cit->info().satisfy_normal_cond_ && !checkUpNormalCondition(cit, tmp_chds, true)) {
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

  void Approximation::refineD(std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &ori_bs_jpts,
                              std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_inner_jpts,
                              std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_outer_jpts,
                              std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_bs_jpts)
  {
    inner_kdtree_.buildTree(inner_jpts_[0].data(), inner_jpts_.size());
    outer_kdtree_.buildTree(outer_jpts_[0].data(), outer_jpts_.size());

    jpts_.reserve(inner_jpts_.size()+outer_jpts_.size());
    size_t pt_size=inner_jpts_.size();
    for(size_t i=0; i<pt_size; ++i) { jpts_.push_back({deformed_inner_jpts[i], inner_jpts_[i], -1, 1}); }
    pt_size=outer_jpts_.size();
    for(size_t i=0; i<pt_size; ++i) { jpts_.push_back({deformed_outer_jpts[i], outer_jpts_[i], 1, 1}); }

    std::vector<std::pair<Point, VertexInfo>> init_vertices;
    init_vertices.reserve(ori_bs_jpts.size());
    for(size_t ind=0; ind<ori_bs_jpts.size(); ++ind) {
      const Eigen::Matrix<zsw::Scalar,3,1> &pt=deformed_bs_jpts[ind];
      init_vertices.push_back(std::make_pair(Point(pt[0],pt[1],pt[2]),VertexInfo(ind, zsw::BBOX_POINT, ori_bs_jpts[ind])));
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
    size_t add_pt_for_err = 0;
    size_t add_pt_for_normal = 0;
    do{
      while(!err_queue.empty()) {
        std::pair<zsw::Scalar, JudgePoint*> jpt_info=err_queue.top(); err_queue.pop();
        zsw::Scalar real_err=fabs(jpt_info.second->val_cur_-jpt_info.second->val_exp_);
        if(fabs(real_err-jpt_info.first) > zsw::const_val::eps) { continue; } // have already updated
        PointType pt_type= (jpt_info.second->val_exp_<0) ? zsw::INNER_POINT : zsw::OUTER_POINT;
        VertexInfo vertex_info(-1, pt_type, jpt_info.second->pt_cur_);
        std::vector<Chd> chds;
        tw_->addPointInDelaunay(jpt_info.second->pt_c_, vertex_info, chds);
        if(++add_pt_for_err%100==0) {
          NZSWLOG("zsw_info") << "add_pt_for_err:" << add_pt_for_err << std::endl;
          // writeTetMesh(tmp_outdir_+"simp_tol_ori_"+ std::to_string(add_pt_for_err) +".vtk", {zsw::ignore_bbox}, nullptr, false);
          // writeTetMesh(tmp_outdir_+"simp_tol_deformed_"+ std::to_string(add_pt_for_err) +".vtk", {zsw::ignore_bbox}, nullptr, true);
        }
        updateJptsInCellsD(chds, err_queue);
      }

      // find a tet viloate the normal condition
      bool viloate_normal_cond=false;
      std::vector<Chd> tmp_chds;
      for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) {
        if(!cit->info().satisfy_normal_cond_ && !checkUpNormalCondition(cit, tmp_chds, false)) {
          if(++add_pt_for_normal%50==0) {
            NZSWLOG("zsw_info") << "add_pt_for_normal:" << add_pt_for_normal << std::endl;
            // writeTetMesh(tmp_outdir_+"simp_tol_ori_n_"+ std::to_string(add_pt_for_normal) +".vtk", {zsw::ignore_bbox}, nullptr, false);
            // writeTetMesh(tmp_outdir_+"simp_tol_deformed_n_"+ std::to_string(add_pt_for_normal) +".vtk", {zsw::ignore_bbox}, nullptr, true);
          }
          viloate_normal_cond=true;
          break;
        }
      }
      if(viloate_normal_cond) { updateJptsInCellsD(tmp_chds, err_queue); }
      else { break; }
    } while(true);
    writeTetMesh(tmp_outdir_ +"before_swap_vertex.vtk", {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
    tw_->swapVertex();
    NZSWLOG("zsw_info") << "refined add pt cnt:" << tw_->getTds().number_of_vertices()-ori_bs_jpts.size()-1 << std::endl;
  }

  bool Approximation::checkUpNormalCondition(Chd chd, std::vector<Chd> &chds, bool using_cur_pts)
  {
    if(!tw_->isTolCell(chd)) { return true; }
    Eigen::Matrix<zsw::Scalar,3,4> tri_pts;
    if(using_cur_pts) {
      tri_pts<<
        chd->vertex(0)->point()[0], chd->vertex(1)->point()[0], chd->vertex(2)->point()[0], chd->vertex(3)->point()[0],
        chd->vertex(0)->point()[1], chd->vertex(1)->point()[1], chd->vertex(2)->point()[1], chd->vertex(3)->point()[1],
        chd->vertex(0)->point()[2], chd->vertex(1)->point()[2], chd->vertex(2)->point()[2], chd->vertex(3)->point()[2];
    } else {
      tri_pts.block<3,1>(0,0) = chd->vertex(0)->info().pos_c_;
      tri_pts.block<3,1>(0,1) = chd->vertex(1)->info().pos_c_;
      tri_pts.block<3,1>(0,2) = chd->vertex(2)->info().pos_c_;
      tri_pts.block<3,1>(0,3) = chd->vertex(3)->info().pos_c_;
    }
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
    writeTetMesh(tmp_outdir_+"before_add_"+std::to_string(ind)+".vtk",  {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out}, nullptr, false);
    abort();
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
    VertexInfo vertex_info(-1, pt_type, jpts_[jpt_ind].pt_cur_);
    size_t nv = tw_->getTds().number_of_vertices();
    //writeTetMesh("/home/wegatron/tmp/before_add_pt.vtk",  {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
    if(using_cur_pts) {  tw_->addPointInDelaunay(jpts_[jpt_ind].pt_cur_, vertex_info, chds);
    } else {  tw_->addPointInDelaunay(jpts_[jpt_ind].pt_c_, vertex_info, chds);  }
    //writeTetMesh("/home/wegatron/tmp/after_add_pt.vtk",  {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
    return tw_->getTds().number_of_vertices() == nv;
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

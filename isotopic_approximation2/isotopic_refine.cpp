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
  void Approximation::refine()
  {
#if 0
    for(size_t i=0; i<1000; ++i) {
      testKdtree(outer_kdtree_, outer_jpts_);
    }
    for(size_t i=0; i<1000; ++i) {
      testKdtree(inner_kdtree_, inner_jpts_);
    }
    std::cerr << "test_pass!!!" << std::endl;
#endif
    std::priority_queue<std::pair<zsw::Scalar,JudgePoint*>,
                        std::vector<std::pair<zsw::Scalar,JudgePoint*>>,
                        ErrorMaxComparison> err_queue;
    std::cout << "[INFO] " << "refine start!" << std::endl;
    const TTds &tds=tw_->getTds();
    bool add_pt_flag=false;
    size_t debug_check_normalcond_time=0;
    do{
      std::cerr << "cell_size:" << tds.number_of_cells() << std::endl;
      std::cerr << "vertices_size:" << tds.number_of_vertices() << std::endl;
      for(size_t i=0; i<jpts_.size();++i) {
        zsw::Scalar err=fabs(jpts_[i].val_cur_-jpts_[i].val_exp_);
        if(err>1) { err_queue.push(std::make_pair(err, &jpts_[i])); }
      }
      if(err_queue.empty()) { break; }
      size_t debug_pt_size=0;
      // zsw::common::ClockC11 clock;
      while(!err_queue.empty()) {
        std::pair<zsw::Scalar,JudgePoint*> jpt_info=err_queue.top(); err_queue.pop();
        zsw::Scalar real_err=fabs(jpt_info.second->val_cur_-jpt_info.second->val_exp_);
        if(fabs(real_err-jpt_info.first) > zsw::const_val::eps) { continue; } // have already updated

        if(++debug_pt_size%100==0) {          std::cerr << "pt_size:" << debug_pt_size << std::endl;        }

        PointType pt_type= (jpt_info.second->val_exp_<0) ? zsw::INNER_POINT : zsw::OUTER_POINT;
        VertexInfo vertex_info(-1, pt_type, jpt_info.second->pt_, 0.0);
        std::vector<Chd> chds;

        // if(debug_pt_size%100==0) { std::cerr << "find pt time cost:" << clock.time() << std::endl;        }

        tw_->addPointInDelaunaySafe(jpt_info.second->pt_, vertex_info, chds);

        // if(debug_pt_size%100==0) {          std::cerr << "add pt into delaunay cost:" << clock.time() << std::endl;        }

        std::vector<std::vector<JudgePoint*>> updated_jpts(chds.size());
#pragma omp parallel for
        for(size_t i=0; i<chds.size(); ++i) {
          if(tw_->isValidCell(chds[i])) { updateJptsInCell(chds[i], &updated_jpts[i]); }
        }
        std::set<JudgePoint*> up_jpt_set;
        for(std::vector<JudgePoint*> up_jpt : updated_jpts) {
          for(JudgePoint * jpt : up_jpt) {            up_jpt_set.insert(jpt);          }
        }
        for(JudgePoint * jpt : up_jpt_set) {
          zsw::Scalar tmp_err=fabs(jpt->val_cur_-jpt->val_exp_);
          if(tmp_err>1.0) { err_queue.push(std::make_pair(tmp_err, jpt)); }
        }
      }

      TTds tmp_tds=tw_->getTds();
      mutuallTessellation(&tmp_tds);
      // if(isZeroTetExist(tmp_tds)) { std::cerr << "Exist Zero tet!!!" << std::endl; }
      writeTetMesh("/home/wegatron/tmp/check_normal_cond_"+std::to_string(debug_check_normalcond_time)+".vtk",
                   {zsw::ignore_bbox, zsw::ignore_out}, &tmp_tds);
      ++debug_check_normalcond_time;
      std::cerr << "start checkup normal cond!!!" << std::endl;
      add_pt_flag=false;
      std::queue<Chd> chds_queue;
      for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) { chds_queue.push(cit); }
      static size_t add_pt_num=0;
      while(!chds_queue.empty()) {
        Chd chd = chds_queue.front(); chds_queue.pop();
        if(!tw_->isTolCell(chd)) { continue; }
        if(!checkUpNormalCondition(chd, chds_queue)) {
          if(++add_pt_num%100==0) { std::cerr << "add pt_num=" << add_pt_num << std::endl; }
          add_pt_flag=true;
        }
      }
      if(!add_pt_flag) { break; }
      std::vector<Chd> chds;
      for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) {
        if(tw_->isValidCell(cit))  { chds.push_back(cit); }
      }
#pragma omp parallel for
      for(size_t i=0; i<chds.size(); ++i) {
        updateJptsInCell(chds[i],nullptr);
      }
      //break;
      // std::cerr << "normal cond check end!!!" << std::endl;
    }while(1);
#if 0
    size_t ind=0;
    for(auto chd=tds.cells_begin(); chd!=tds.cells_end(); ++chd) {
      if(!tw_->isTolCell(chd)) { continue; }
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
      if(!normalCondition(val, scaled_tri_pts, tri_pts, inner_jpts_, outer_jpts_, inner_kdtree_, outer_kdtree_)) {
        std::string filepath=tmp_outdir_+"normal_cond/nc_"+std::to_string(ind)+".vtk";
        if(!normalCondition(val, scaled_tri_pts, tri_pts, inner_jpts_, outer_jpts_, inner_kdtree_, outer_kdtree_,
                            true, &filepath)) {
          std::cerr << "Refine result normal condition failed on " << ind << std::endl;
          //abort();
        }
        ++ind;
      }
    }
#endif
  }

  JudgePoint * Approximation::findMaxErrorJpt(zsw::Scalar &error)
  {
    JudgePoint *ret=nullptr;
    error=1.0;
    for(JudgePoint &jpt : jpts_) {
      zsw::Scalar tmp_error=fabs(jpt.val_cur_-jpt.val_exp_);
      if(tmp_error > error) {
        error=tmp_error;
        ret=&jpt;
      }
    }
    return ret;
  }

  bool Approximation::checkUpNormalCondition(Chd chd, std::queue<Chd> &chds_queue)
  {
    assert(tw_->isTolCell(chd));
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
    std::vector<Chd> chds;
    //writeTetMesh("/home/wegatron/tmp/before_add_pt.vtk",  {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
    tw_->addPointInDelaunaySafe(jpts_[jpt_ind].pt_, vertex_info, chds);
    //writeTetMesh("/home/wegatron/tmp/after_add_pt.vtk",  {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
    for(Chd chd : chds) { chds_queue.push(chd); }
    return false;
  }

  bool Approximation::checkNormalCondition() const
  {
    const TTds &tds=tw_->getTds();
    size_t normal_cond_debug_i=0;
    for(auto chd=tds.cells_begin(); chd!=tds.cells_end(); ++chd) {
      if(!tw_->isTolCell(chd)) { continue; }
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
      std::string filepath(tmp_outdir_+"normal_cond/check_"+std::to_string(normal_cond_debug_i++)+".vtk");
      if(!normalCondition(val, scaled_tri_pts, tri_pts, inner_jpts_, outer_jpts_,
                          inner_kdtree_, outer_kdtree_,
                          true, &filepath)) {
        std::cerr << "normal cond failed!!!" << std::endl;
        std::cerr << "normal_cond_debug_i=" << normal_cond_debug_i-1 << std::endl;
        return false;
      }
    }
    return true;
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
#if 1
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

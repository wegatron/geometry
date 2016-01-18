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
  void Approximation::refine()
  {
    std::priority_queue<std::pair<zsw::Scalar,JudgePoint*>,
                        std::vector<std::pair<zsw::Scalar,JudgePoint*>>,
                        ErrorMaxComparison> err_queue;
    std::cout << "[INFO] " << "refine start!" << std::endl;
    size_t pre=0, cur=1;
    std::unordered_set<std::string> cell_key_set[2];
    tw_->initCellKeySet(cell_key_set[0]);
    const TTds &tds=tw_->getTds();
    bool add_pt_flag=false;
    do{
      std::cerr << "key_set_size:" << cell_key_set[pre].size() << std::endl;
      std::cerr << "cell_size:" << tds.number_of_cells() << std::endl;
      std::cerr << "vertices_size:" << tds.number_of_vertices() << std::endl;
      for(size_t i=0; i<jpts_.size();++i) {
        zsw::Scalar err=fabs(jpts_[i].val_cur_-jpts_[i].val_exp_);
        if(err>1) { err_queue.push(std::make_pair(err, &jpts_[i])); }
      }
      if(err_queue.empty()) { break; }
      while(!err_queue.empty()) {
        std::pair<zsw::Scalar,JudgePoint*> jpt_info=err_queue.top(); err_queue.pop();
        zsw::Scalar real_err=fabs(jpt_info.second->val_cur_-jpt_info.second->val_exp_);
        if(fabs(real_err-jpt_info.first) > zsw::const_val::eps) { continue; } // have already updated
        PointType pt_type= (jpt_info.second->val_exp_<0) ? zsw::INNER_POINT : zsw::OUTER_POINT;
        VertexInfo vertex_info(-1, pt_type, jpt_info.second->pt_, 0.0);
        std::vector<Chd> chds;
        tw_->addPointInDelaunaySafe(jpt_info.second->pt_, vertex_info, chds, &cell_key_set[pre], &cell_key_set[cur]);
        swap(pre,cur);
        for(Chd chd : chds) { if(tw_->isValidCell(chd)) {updateJptsInCell(chd, &err_queue);} }
      }
      add_pt_flag=false;
      std::queue<Chd> chds_queue;
      for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) { chds_queue.push(cit); }
      while(!chds_queue.empty()) {
        Chd chd = chds_queue.front(); chds_queue.pop();
        if(!tw_->isTolCell(chd)) { continue; }
        cell_key_set[cur].clear();
        if(!checkUpNormalCondition(chd, chds_queue, &cell_key_set[pre], &cell_key_set[cur])) { add_pt_flag=true; }
        swap(pre,cur);
      }
      if(!add_pt_flag) { break; }
      for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) {
        if(tw_->isValidCell(cit)) { updateJptsInCell(cit,nullptr); }
      }
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
      std::string filepath=tmp_outdir_+"normal_cond/nc_"+std::to_string(ind++)+".vtk";
      if(!normalCondition(val, scaled_tri_pts, tri_pts, inner_jpts_, outer_jpts_, inner_kdtree_ptr_, outer_kdtree_ptr_,
                          true, &filepath)) {
        std::cerr << "Refine result normal condition failed!!!" << std::endl;
        abort();
      }
    }
#endif
  }

  bool Approximation::checkUpNormalCondition(Chd chd, std::queue<Chd> &chds_queue,
                                             std::unordered_set<std::string> *cell_key_set_pre,
                                             std::unordered_set<std::string> *cell_key_set_cur)
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
    if(normalCondition(val, scaled_tri_pts, tri_pts, inner_jpts_, outer_jpts_, inner_kdtree_ptr_, outer_kdtree_ptr_)) {
      return true;
    }

    // static size_t debug_step=0;
    // std::cerr << "tri_pts:\n" << tri_pts << std::endl;
    //  writeTetMesh("/home/wegatron/tmp/upnorm_"+std::to_string(debug_step++)+".vtk",
    //               {ignore_self_in, ignore_bbox, ignore_self_out});

    // insert point into delaunay triangulation
    Eigen::Matrix<zsw::Scalar,3,1> center;
    calcCircumcenter(tri_pts, center);
    // find the minest
    std::vector<size_t> in_indices;
    std::vector<zsw::Scalar> in_dist;
    inner_kdtree_ptr_->queryNearest(center, in_indices, in_dist);
    std::vector<size_t> out_indices;
    std::vector<zsw::Scalar> out_dist;
    outer_kdtree_ptr_->queryNearest(center, out_indices, out_dist);
    size_t jpt_ind = (in_dist[0]<out_dist[0]) ? in_indices[0] : out_indices[0]+inner_jpts_.size();
    // check if the point is already in
    if(fabs(jpts_[jpt_ind].val_cur_-jpts_[jpt_ind].val_exp_)<zsw::const_val::eps) {
      // already in
#if 0
      static int ind=0;
      std::vector<Eigen::Matrix<zsw::Scalar,3,1>> pts;
      std::vector<Eigen::Matrix<zsw::Scalar,3,4>> cells;
      cells.push_back(tri_pts);
      pts.push_back(center); pts.push_back(jpts_[jpt_ind].pt_);
      writeCellsAndPoints(tmp_outdir_+"normal_cond/not_st_nc_"+std::to_string(ind++)+".vtk", cells, pts);
#endif
      return true;
    }
    // add and update
    PointType pt_type=(jpts_[jpt_ind].val_exp_<0) ? zsw::INNER_POINT : zsw::OUTER_POINT;
    VertexInfo vertex_info(-1, pt_type, jpts_[jpt_ind].pt_, 0.0);
    std::vector<Chd> chds;
    tw_->addPointInDelaunaySafe(jpts_[jpt_ind].pt_, vertex_info, chds, cell_key_set_pre, cell_key_set_cur);
    //for(Chd chd : chds) { updateJptsInCell(chd, nullptr); }
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
                          inner_kdtree_ptr_, outer_kdtree_ptr_,
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

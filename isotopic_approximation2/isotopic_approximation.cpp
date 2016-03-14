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
#include <zswlib/zsw_clock_c11.h>

using namespace std;

namespace zsw{

  void Approximation::init(const zsw::Scalar err_epsilon,
                           const zsw::Scalar tri_sample_r,
                           const zsw::Scalar tet_sample_r,
                           std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &inner_jpts,
                           std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &outer_jpts,
                           std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bs_jpts)
  {
    tri_sample_r_=tri_sample_r;
    tet_sample_r_=tet_sample_r;
    inner_jpts_=inner_jpts; outer_jpts_=outer_jpts;
    NZSWLOG("zsw_info") << "inner judge point size:" << inner_jpts_.size() << std::endl;
    NZSWLOG("zsw_info") << "outer judge point size:" << outer_jpts_.size() << std::endl;
    jpts_.reserve(inner_jpts_.size()+outer_jpts_.size());
    clock_.clearCur();
    refine(bs_jpts);
    NZSWLOG("zsw_info") << "refine time cost" << clock_.time() << std::endl;
    NZSWLOG("zsw_info") << "refine complete, init finished!" << std::endl;
    NZSWLOG("zsw_info") << "vertex size:" << tw_->getTds().number_of_vertices() << std::endl;
    NZSWLOG("zsw_info") << "cell size:" << tw_->getTds().number_of_cells() << std::endl;
  }

  void Approximation::initD(const zsw::Scalar err_epsilon,
                            const zsw::Scalar tri_sample_r,
                            const zsw::Scalar tet_sample_r,
                            std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &ori_inner_jpts,
                            std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &ori_outer_jpts,
                            std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &ori_bs_jpts,
                            std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_inner_jpts,
                            std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_outer_jpts,
                            std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_bs_jpts)
  {
    tri_sample_r_=tri_sample_r;
    tet_sample_r_=tet_sample_r;
    inner_jpts_ = ori_inner_jpts; outer_jpts_ = ori_outer_jpts;
    NZSWLOG("zsw_info") << "inner judge point size:" << inner_jpts_.size() << std::endl;
    NZSWLOG("zsw_info") << "outer judge point size:" << outer_jpts_.size() << std::endl;
    jpts_.reserve(inner_jpts_.size()+outer_jpts_.size());
    refineD(ori_bs_jpts, deformed_inner_jpts, deformed_outer_jpts, deformed_bs_jpts);
    NZSWLOG("zsw_info") << "refine time cost" << clock_.time() << std::endl;
    NZSWLOG("zsw_info") << "refine complete, init finished!" << std::endl;
    NZSWLOG("zsw_info") << "vertex size:" << tw_->getTds().number_of_vertices() << std::endl;
    NZSWLOG("zsw_info") << "cell size:" << tw_->getTds().number_of_cells() << std::endl;
  }

  void  Approximation::simp(const std::string &tmp_output_dir)
  {
    writeTetMesh(tmp_output_dir+"before_simp_tol.vtk", {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
    writeTetMesh(tmp_output_dir+"before_simp_tol_deform.vtk", {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out},
                 nullptr, false);
    TTds tmp_tds=tw_->getTds();
    mutuallTessellation();
    writeTetMesh(tmp_output_dir+"after_refine_zero_surf.vtk", {zsw::ignore_bbox, zsw::ignore_out});
    writeTetMesh(tmp_output_dir+"after_mt_tol.vtk", {zsw::ignore_bbox, zsw::ignore_self_out, zsw::ignore_self_in});
    tw_->setTds(tmp_tds);
    clock_.clearCur();
    simpTolerance();
    NZSWLOG("zsw_info") << "simp_tol time cost:" << clock_.time() << std::endl;
    writeTetMesh(tmp_output_dir+"after_simp_tol.vtk", {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
    writeTetMesh(tmp_output_dir+"after_simp_deform_tol.vtk", {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out},
                 nullptr, false);
    mutuallTessellation();
    NZSWLOG("zsw_info") << "mutuall tessellation time cost:" << clock_.time() << std::endl;
    writeTetMesh(tmp_output_dir+"after_simp_tol_zero_surf.vtk", {zsw::ignore_bbox, zsw::ignore_out});
    simpZeroSurface();

    NZSWLOG("bz_info") << "bz_krj_need_judge_cnt:" << bz_krj_need_judge_cnt_ << std::endl;
    NZSWLOG("bz_info") << "bz_normal_cond_judge_cnt:" << bz_normal_cond_judge_cnt_ << std::endl;
    NZSWLOG("bz_info") << "bz_error_bound_judge_cnt:" << bz_error_bound_judge_cnt_ << std::endl;
    NZSWLOG("bz_info") << "bz_judge_pt_cnt:" << bz_judge_pt_cnt_ << std::endl;

    NZSWLOG("zsw_info") << "simp_zero surf time cost:" << clock_.time() << std::endl;
    writeTetMesh(tmp_output_dir+"simped_zero_surf.vtk", {zsw::ignore_bbox, zsw::ignore_out});
    simpBZEdges();

    NZSWLOG("bz_info") << "bz_krj_need_judge_cnt:" << bz_krj_need_judge_cnt_ << std::endl;
    NZSWLOG("bz_info") << "bz_normal_cond_judge_cnt:" << bz_normal_cond_judge_cnt_ << std::endl;
    NZSWLOG("bz_info") << "bz_error_bound_judge_cnt:" << bz_error_bound_judge_cnt_ << std::endl;
    NZSWLOG("bz_info") << "bz_judge_pt_cnt:" << bz_judge_pt_cnt_ << std::endl;

    NZSWLOG("zsw_info") << "simp_bz time cost:" << clock_.time() << std::endl;
    NZSWLOG("zsw_info") << "total time cost:" << clock_.totalTime() << std::endl;
    NZSWLOG("zsw_info") << "final point count:" << countZeroPoints() << std::endl;
  }

  size_t Approximation::countZeroPoints() const
  {
    size_t cnt=0;
    const TTds &tds=tw_->getTds();
    for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
      if(vit->info().pt_type_==zsw::ZERO_POINT) { ++cnt; }
    }
    return cnt;
  }

  void Approximation::updateJptsInCell(Chd chd,
                                       std::vector<JudgePoint*> *updated_jpts,
                                       bool using_cur_pts)
  {
    if(!(chd->vertex(0)->info().pt_type_!=zsw::INVALID_POINT && chd->vertex(1)->info().pt_type_!=zsw::INVALID_POINT
         && chd->vertex(2)->info().pt_type_!=zsw::INVALID_POINT && chd->vertex(3)->info().pt_type_!=zsw::INVALID_POINT))
      { return; }
    // calc jpts in bbox
    std::vector<const JudgePoint*> jpts_in_bbox;
    Vhd vhds[4] = {chd->vertex(0), chd->vertex(1), chd->vertex(2), chd->vertex(3)};
    calcJptsInBbox(vhds,4,jpts_in_bbox, using_cur_pts);
    Eigen::Matrix<zsw::Scalar,4,4> A;
    if(using_cur_pts) {
      A <<
        chd->vertex(0)->point()[0], chd->vertex(1)->point()[0], chd->vertex(2)->point()[0], chd->vertex(3)->point()[0],
        chd->vertex(0)->point()[1], chd->vertex(1)->point()[1], chd->vertex(2)->point()[1], chd->vertex(3)->point()[1],
        chd->vertex(0)->point()[2], chd->vertex(1)->point()[2], chd->vertex(2)->point()[2], chd->vertex(3)->point()[2],
        1,1,1,1;
    } else {
      A.block<3,1>(0,0) = chd->vertex(0)->info().pos_c_;
      A.block<3,1>(0,1) = chd->vertex(1)->info().pos_c_;
      A.block<3,1>(0,2) = chd->vertex(2)->info().pos_c_;
      A.block<3,1>(0,3) = chd->vertex(3)->info().pos_c_;
      A(3,0) = A(3,1) = A(3,2) = A(3,3) = 1.0;
    }
    Eigen::PartialPivLU<Eigen::Matrix<zsw::Scalar,4,4>> pplu; pplu.compute(A);
    for(const JudgePoint * jpt : jpts_in_bbox) {
      Eigen::Matrix<zsw::Scalar,4,1> x, b;
      b.block<3,1>(0,0) = jpt->pt_cur_; b[3] = 1.0;
      x=pplu.solve(b);
      if(x[0]<-zsw::const_val::eps || x[1]<-zsw::const_val::eps || x[2]<-zsw::const_val::eps
         || x[3]<-zsw::const_val::eps) { continue; } // jpt not in this cell
      assert((A*x-b).norm()<zsw::const_val::eps);
      Eigen::Matrix<zsw::Scalar,4,1> val;
      for(size_t vi=0; vi<4; ++vi) {
        if(chd->vertex(vi)->info().pt_type_==zsw::INNER_POINT) { val[vi]=-1; }
        else if(chd->vertex(vi)->info().pt_type_==zsw::ZERO_POINT){ val[vi]=0; }
        else { val[vi]=1; }
      }
      JudgePoint *tmp_jpt = const_cast<JudgePoint*>(jpt);
      tmp_jpt->val_cur_=val.dot(x);
      if(updated_jpts!=nullptr) {      updated_jpts->push_back(tmp_jpt); }
      // if(err_queue!=nullptr) {
      //   zsw::Scalar err=fabs(tmp_jpt->val_cur_-tmp_jpt->val_exp_);
      //   if(err>1) { err_queue->push(std::make_pair(err, tmp_jpt)); }
      // }
    }
  }

  void Approximation::updateJptsInCells(const std::vector<Chd> &chds,     std::priority_queue<std::pair<zsw::Scalar,JudgePoint*>,
                                        std::vector<std::pair<zsw::Scalar,JudgePoint*>>,
                                        ErrorMaxComparison> &err_queue)
  {
    FUNCTION_TIME_ANALYSIS();
    std::vector<std::vector<JudgePoint*>> updated_jpts(chds.size());
#pragma omp parallel for
    for(size_t i=0; i<chds.size(); ++i) {
      if(tw_->isValidCell(chds[i])) {
        updateJptsInCell(chds[i], &updated_jpts[i], true);
        chds[i]->info().satisfy_normal_cond_=false;
      }
    }
    std::set<JudgePoint*> up_jpt_set;
    for(std::vector<JudgePoint*> up_jpt : updated_jpts) {
      for(JudgePoint * jpt : up_jpt) {            up_jpt_set.insert(jpt);          }
    }
    for(JudgePoint * jpt : up_jpt_set) {
      zsw::Scalar tmp_err=fabs(jpt->val_cur_-jpt->val_exp_);
      if(tmp_err > 1.0 - alpha_) { err_queue.push(std::make_pair(tmp_err, jpt)); }
    }
  }

  void Approximation::updateJptsInCellsD(const std::vector<Chd> &chds,     std::priority_queue<std::pair<zsw::Scalar,JudgePoint*>,
                                         std::vector<std::pair<zsw::Scalar,JudgePoint*>>,
                                         ErrorMaxComparison> &err_queue)
  {
    FUNCTION_TIME_ANALYSIS();
    std::vector<std::vector<JudgePoint*>> updated_jpts(chds.size());
#pragma omp parallel for
    for(size_t i=0; i<chds.size(); ++i) {
      if(tw_->isValidCell(chds[i])) {
        updateJptsInCell(chds[i], &updated_jpts[i], false);
        chds[i]->info().satisfy_normal_cond_=false;
      }
    }

    std::set<JudgePoint*> up_jpt_set;
    for(std::vector<JudgePoint*> up_jpt : updated_jpts) {
      for(JudgePoint * jpt : up_jpt) {            up_jpt_set.insert(jpt);          }
    }
    for(JudgePoint * jpt : up_jpt_set) {
      zsw::Scalar tmp_err=fabs(jpt->val_cur_-jpt->val_exp_);
      if(tmp_err > 1.0 - alpha_) { err_queue.push(std::make_pair(tmp_err, jpt)); }
    }
  }

  bool Approximation::isSatisfyErrorBound(const std::vector<VertexTriple> &bound_tris,
                                          const std::vector<const JudgePoint*> &jpts_in_bbox,
                                          const Eigen::Matrix<zsw::Scalar,3,1> &merge_pt,
                                          const zsw::Scalar v_pt,
                                          std::vector<JudgePointUpdateData> * jup)
  {
    FUNCTION_TIME_ANALYSIS();
    ++bz_judge_pt_cnt_;
    std::vector<bool> is_updated(jpts_in_bbox.size(), false);
    size_t false_cnt=0;
    for(const VertexTriple &vt : bound_tris) {
      assert(vt.first->info().pt_type_!=zsw::INVALID_POINT);
      assert(vt.second->info().pt_type_!=zsw::INVALID_POINT);
      assert(vt.third->info().pt_type_!=zsw::INVALID_POINT);

      Eigen::Matrix<zsw::Scalar,4,4> A;
      A(0,0)=vt.first->point()[0]; A(1,0)=vt.first->point()[1];  A(2,0)=vt.first->point()[2]; A(3,0)=1;
      A(0,1)=vt.second->point()[0]; A(1,1)=vt.second->point()[1];  A(2,1)=vt.second->point()[2]; A(3,1)=1;
      A(0,2)=vt.third->point()[0]; A(1,2)=vt.third->point()[1];  A(2,2)=vt.third->point()[2]; A(3,2)=1;
      A.block<3,1>(0,3)=merge_pt; A(3,3)=1;
      Eigen::Matrix<zsw::Scalar,4,1> val;
      if(vt.first->info().pt_type_==zsw::INNER_POINT) { val[0]=-1; }
      else if(vt.first->info().pt_type_==zsw::ZERO_POINT) { val[0]=0; }
      else { val[0]=1; }
      if(vt.second->info().pt_type_==zsw::INNER_POINT) { val[1]=-1; }
      else if(vt.second->info().pt_type_==zsw::ZERO_POINT) { val[1]=0; }
      else { val[1]=1; }
      if(vt.third->info().pt_type_==zsw::INNER_POINT) { val[2]=-1; }
      else if(vt.third->info().pt_type_==zsw::ZERO_POINT) { val[2]=0; }
      else { val[2]=1; }
      val[3]=v_pt;
      Eigen::PartialPivLU<Eigen::Matrix<zsw::Scalar,4,4>> pplu; pplu.compute(A);
      std::vector<zsw::Scalar> n_vals(jpts_in_bbox.size(), -100);
#pragma omp parallel for
      for(size_t jpt_i=0; jpt_i<jpts_in_bbox.size(); ++jpt_i) {
        if(false_cnt!=0 || is_updated[jpt_i]) { continue; }
        const JudgePoint *jpt=jpts_in_bbox[jpt_i];
        Eigen::Matrix<zsw::Scalar,4,1> b; b.block<3,1>(0,0)=jpt->pt_cur_; b[3]=1;
        Eigen::Matrix<zsw::Scalar,4,1> x=pplu.solve(b);
        if(x[0]<-zsw::const_val::eps || x[1]<-zsw::const_val::eps
           || x[2]<-zsw::const_val::eps || x[3]<-zsw::const_val::eps) { continue; }
        is_updated[jpt_i]=true;
        zsw::Scalar tmp_val = val.dot(x);
#pragma omp atomic
        false_cnt += (fabs(tmp_val-jpt->val_exp_)>1-alpha_);
        if(jup!=nullptr) { n_vals[jpt_i]=tmp_val; }
      }
      if(false_cnt!=0) { break; }
      if(jup!=nullptr) {
        for(size_t i=0; i<jpts_in_bbox.size(); ++i) {
          if(n_vals[i]<-10.0) { continue; }
          jup->push_back({const_cast<JudgePoint*>(jpts_in_bbox[i]), n_vals[i]});
        }
      }
    }
    return false_cnt==0;
  }

  void Approximation::sampleAdjCells(const TTds::Edge &e, std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &sample_points) const
  {
    FUNCTION_TIME_ANALYSIS();
    Vhd vhds[2]={e.first->vertex(e.second), e.first->vertex(e.third)};
    std::vector<Chd> cells;
    const TTds &tds=tw_->getTds();
    tds.incident_cells(e.first->vertex(e.second), std::back_inserter(cells));
    std::map<std::string,Chd> cell_map;
    for(Chd chd : cells) {
      std::string key_str = cell2key(chd);
      cell_map[key_str]=chd;
    }
    cells.clear();
    tds.incident_cells(e.first->vertex(e.third), std::back_inserter(cells));
    for(Chd chd : cells) {
      std::string key_str = cell2key(chd);
      cell_map[key_str]=chd;
    }
    for(auto it=cell_map.begin(); it!=cell_map.end(); ++it) {
      Eigen::Matrix<zsw::Scalar,3,1> v0,v1,v2,v3;
      v0<< it->second->vertex(0)->point()[0],it->second->vertex(0)->point()[1],it->second->vertex(0)->point()[2];
      v1<< it->second->vertex(1)->point()[0],it->second->vertex(1)->point()[1],it->second->vertex(1)->point()[2];
      v2<< it->second->vertex(2)->point()[0],it->second->vertex(2)->point()[1],it->second->vertex(2)->point()[2];
      v3<< it->second->vertex(3)->point()[0],it->second->vertex(3)->point()[1],it->second->vertex(3)->point()[2];
      if(sample_tet_by_ref_) {
        Eigen::Matrix<zsw::Scalar,3,1> rv0 = it->second->vertex(0)->info().pos_c_;
        Eigen::Matrix<zsw::Scalar,3,1> rv1 = it->second->vertex(1)->info().pos_c_;
        Eigen::Matrix<zsw::Scalar,3,1> rv2 = it->second->vertex(2)->info().pos_c_;
        Eigen::Matrix<zsw::Scalar,3,1> rv3 = it->second->vertex(3)->info().pos_c_;
        zsw::Scalar ref_tet_sample_r = tet_sample_r_ * g_scale_;
        sampleTetRefTet(v0, v1, v2, v3, rv0, rv1, rv2, rv3, ref_tet_sample_r, sample_points);
      } else { sampleTet(v0,v1,v2,v3,tet_sample_r_, sample_points); }
    }
#if 0
    // cout adj cells
    writeAdjcentCells(tmp_outdir_+"adj_cells.vtk", e);
    writePoints(tmp_outdir_+"sp.vtk", sample_points);
    abort();
#endif
  }

  // void sampleKrj(const zsw::KernelRegionJudger &krj, std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &sample_points) const
  // {
  //   std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  // }

  void Approximation::calcJptsInBbox(Vhd *vhd_ptr, const size_t n, std::vector<const JudgePoint*> &jpts_in_bbox,
                                     bool using_cur_pts) const
  {
    Eigen::Matrix<zsw::Scalar,3,2> bbox;
    if(using_cur_pts) {
      calcVerticesBbox(vhd_ptr, n, bbox);
    } else {
      calcVerticesBboxD(vhd_ptr, n, bbox);
    }
    for(const JudgePoint &jpt : jpts_) {
      if(jpt.pt_cur_[0]<bbox(0,0) || jpt.pt_cur_[1]<bbox(1,0) || jpt.pt_cur_[2]<bbox(2,0) ||
         jpt.pt_cur_[0]>bbox(0,1) || jpt.pt_cur_[1]>bbox(1,1) || jpt.pt_cur_[2]>bbox(2,1)) { continue; }
      jpts_in_bbox.push_back(&jpt);
    }
  }

  void Approximation::constructKernelRegionJudger(const std::vector<VertexTriple> &bound_tris,
                                                  std::vector<Vhd> &opposite_vs, KernelRegionJudger &krj) const
  {
    for(size_t fi=0; fi<bound_tris.size(); ++fi) {
      Eigen::Matrix<zsw::Scalar,3,1> v0;
      Eigen::Matrix<zsw::Scalar,3,1> v[3];
      v0[0]=opposite_vs[fi]->point()[0]; v0[1]=opposite_vs[fi]->point()[1]; v0[2]=opposite_vs[fi]->point()[2];
      v[0]<< bound_tris[fi].first->point()[0], bound_tris[fi].first->point()[1], bound_tris[fi].first->point()[2];
      v[1]<< bound_tris[fi].second->point()[0], bound_tris[fi].second->point()[1], bound_tris[fi].second->point()[2];
      v[2]<< bound_tris[fi].third->point()[0], bound_tris[fi].third->point()[1], bound_tris[fi].third->point()[2];
      krj.addConstraint(v[0],v[1],v[2],v0);
    }
  }

  void calcVerticesBbox(Vhd *vhd_ptr, const size_t n, Eigen::Matrix<zsw::Scalar,3,2> &bbox)
  {
    bbox(0,0)=bbox(0,1)=(*vhd_ptr)->point()[0];
    bbox(1,0)=bbox(1,1)=(*vhd_ptr)->point()[1];
    bbox(2,0)=bbox(2,1)=(*vhd_ptr)->point()[2];
    for(size_t i=1; i<n; ++i) {
      ++vhd_ptr;
      for(size_t di=0; di<3; ++di) {
        if((*vhd_ptr)->point()[di] < bbox(di,0)) { bbox(di,0)=(*vhd_ptr)->point()[di]; }
        else if((*vhd_ptr)->point()[di] > bbox(di,1)) { bbox(di,1)=(*vhd_ptr)->point()[di]; }
      }
    }
    Eigen::Matrix<zsw::Scalar,3,1> eps_mat; eps_mat<<zsw::const_val::eps, zsw::const_val::eps, zsw::const_val::eps;
    bbox.block<3,1>(0,0) = bbox.block<3,1>(0,0)-eps_mat;
    bbox.block<3,1>(0,1) = bbox.block<3,1>(0,1)+eps_mat;
  }

  void calcVerticesBboxD(Vhd *vhd_ptr, const size_t n, Eigen::Matrix<zsw::Scalar,3,2> &bbox)
  {
    bbox.block<3,1>(0,0) = (*vhd_ptr)->info().pos_c_;
    bbox.block<3,1>(0,1) = (*vhd_ptr)->info().pos_c_;
    for(size_t i=1; i<n; ++i) {
      ++vhd_ptr;
      for(size_t di=0; di<3; ++di) {
        if((*vhd_ptr)->info().pos_c_[di] < bbox(di,0)) { bbox(di,0)=(*vhd_ptr)->info().pos_c_[di]; }
        else if((*vhd_ptr)->info().pos_c_[di]> bbox(di,1)) { bbox(di,1)=(*vhd_ptr)->info().pos_c_[di]; }
      }
    }
    Eigen::Matrix<zsw::Scalar,3,1> eps_mat; eps_mat<<zsw::const_val::eps, zsw::const_val::eps, zsw::const_val::eps;
    bbox.block<3,1>(0,0) = bbox.block<3,1>(0,0)-eps_mat;
    bbox.block<3,1>(0,1) = bbox.block<3,1>(0,1)+eps_mat;
  }

  void Approximation::testTdsValid()
  {
    const TTds &tds=tw_->getTds();
    if(!tds.is_valid()) { abort(); }
    for(auto vit = tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
      std::list<TTds::Cell_handle> cells;
      tds.incident_cells(vit, std::back_inserter(cells));
    }
  }

  void Approximation::testCollapse()
  {
    const TTds &tds = tw_->getTds();
    // for(auto eit=tds.edges_begin(); eit!=tds.edges_end(); ++eit) {
    //   if(tw_->isBoundaryEdge(*eit)) {
    //     const Point &pt = eit->first->vertex(eit->second)->point();
    //     Eigen::Matrix<zsw::Scalar,3,1> merge_pt;
    //     merge_pt<< pt[0], pt[1], pt[2];
    //     Vhd vhd = eit->first->vertex(eit->second);
    //     TTds::Edge e = *eit;
    //     tw_->collapseEdge(e, vhd, merge_pt);
    //     break;
    //   }
    // }

    std::cerr << "number of vertices:" << tds.number_of_vertices() << std::endl;
    std::cerr << "num of cell:" << tds.number_of_cells() << std::endl;

    // for(auto eit=tds.edges_begin(); eit!=tds.edges_end(); ++eit) {
    //   Vhd vhds[2] = {eit->first->vertex(eit->second), eit->first->vertex(eit->third)};
    //   if((vhds[0]->info().pt_type_==zsw::OUTER_POINT && vhds[1]->info().pt_type_==zsw::OUTER_POINT)) {
    //     TTds::Edge edge=*eit;
    //     Point &pa = edge.first->vertex(edge.second)->point();
    //     Point &pb = edge.first->vertex(edge.third)->point();
    //     Point pt((pa[0]+pb[0])/2.0, (pa[1]+pb[1])/2.0, (pa[2]+pb[2])/2.0);
    //     tw_->insertInEdge(edge, pt, zsw::OUTER_POINT);
    //     break;
    //   }
    // }
    mutuallTessellation();

    std::cerr << "number of vertices:" << tds.number_of_vertices() << std::endl;
    std::cerr << "num of cell:" << tds.number_of_cells() << std::endl;
    testTdsValid();
  }
}

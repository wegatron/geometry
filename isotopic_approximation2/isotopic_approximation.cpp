#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_set>
#include <zswlib/error_ctrl.h>
#include "isotopic_approximation.h"
#include "basic_op.h"
#include "bound_sphere.h"
#include "sampling.h"

using namespace std;

namespace zsw{

  // void Approximation::init(const zsw::Scalar &surf_sample_r, const zsw::Scalar &tet_sample_r,
  //                          const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bi_vertices,
  //                          const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bo_vertices)
  // {
  //   surf_sample_r_=surf_sample_r;
  //   tet_sample_r_=tet_sample_r;
  //   Eigen::Matrix<zsw::Scalar,3,2> bbox;
  //   calcBBOX(bo_vertices, bbox);
  //   zsw::Scalar scale = 0.5*(bbox.block<3,1>(0,0) - bbox.block<3,1>(0,1)).norm();
  //   Eigen::Matrix<zsw::Scalar,3,1> transform = 0.5*(bbox.block<3,1>(0,0)+bbox.block<3,1>(0,1));
  //   zsw::BoundSphere bs("/home/wegatron/workspace/geometry/data/bound_sphere.obj", scale, transform);
  //   const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bs_vertices = bs.getVertices();
  //   std::vector<std::pair<Point, VertexInfo>> vertices;
  //   size_t vid=0;
  //   for(const Eigen::Matrix<zsw::Scalar,3,1> &v : bs_vertices) {
  //     vertices.push_back({Point(v[0], v[1], v[2]), VertexInfo(vid++, zsw::BBOX_POINT, v, 0.0)});
  //   }
  //   for(const Eigen::Matrix<zsw::Scalar,3,1> &v : bi_vertices) {
  //     vertices.push_back({Point(v[0],v[1],v[2]), VertexInfo(vid++, zsw::INNER_POINT, v, 0.0)});
  //   }
  //   for(const Eigen::Matrix<zsw::Scalar,3,1> &v : bo_vertices) {
  //     vertices.push_back({Point(v[0],v[1],v[2]), VertexInfo(vid++, zsw::OUTER_POINT, v, 0.0)});
  //   }
  //   tw_.reset(new zsw::TriangulationWapper(vertices));
  //   createJudgePoints();
  //   std::cout << "vertices :" << vertices.size() << std::endl;
  //   std::cout << "judgepoints:" << jpts_.size() << std::endl;
  // }

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
    inner_kdtree_ptr_.reset(new zsw::Flann<zsw::Scalar>(inner_jpts_[0].data(), inner_jpts_.size()));
    outer_kdtree_ptr_.reset(new zsw::Flann<zsw::Scalar>(outer_jpts_[0].data(), outer_jpts_.size()));
    jpts_.reserve(inner_jpts_.size()+outer_jpts_.size());
    for(const Eigen::Matrix<zsw::Scalar,3,1> &in_jpt : inner_jpts_) { jpts_.push_back({in_jpt, -1, 1}); }
    for(const Eigen::Matrix<zsw::Scalar,3,1> &out_jpt : outer_jpts_) { jpts_.push_back({out_jpt, 1, 1}); }
    std::vector<std::pair<Point, VertexInfo>> init_vertices;
    init_vertices.reserve(bs_jpts.size());
    for(size_t ind=0; ind<bs_jpts.size(); ++ind) {
      const Eigen::Matrix<zsw::Scalar,3,1> &pt=bs_jpts[ind];
      init_vertices.push_back(std::make_pair(Point(pt[0],pt[1],pt[2]),VertexInfo(ind, zsw::BBOX_POINT, pt, 0)));
    }
    tw_.reset(new TriangulationWapper(init_vertices));
    std::cout << "[INFO] inner judge point size:" << inner_jpts_.size() << std::endl;
    std::cout << "[INFO] outer judge point size:" << outer_jpts_.size() << std::endl;
    refine();
    std::cout << "[INFO] refine complete, init finished!" << std::endl;
    std::cout << "[INFO] vertex size:" << tw_->getTds().number_of_vertices() << std::endl;
    std::cout << "[INFO] cell size:" << tw_->getTds().number_of_cells() << std::endl;
  }

  void Approximation::refine()
  {
    std::priority_queue<std::pair<zsw::Scalar,JudgePoint*>,
                        std::vector<std::pair<zsw::Scalar,JudgePoint*>>,
                        ErrorMaxComparison> err_queue;
    for(size_t i=0; i<jpts_.size();++i) {
      zsw::Scalar err=fabs(jpts_[i].val_cur_-jpts_[i].val_exp_);
      if(err>1) { err_queue.push(std::make_pair(err, &jpts_[i])); }
    }

    std::cout << "[INFO] " << "refine start!" << std::endl;
    size_t pre=0, cur=1;
    std::unordered_set<std::string> cell_key_set[2];
    tw_->initCellKeySet(cell_key_set[0]);
    const TTds &tds=tw_->getTds();
    while(!err_queue.empty()) {
      std::pair<zsw::Scalar,JudgePoint*> jpt_info=err_queue.top(); err_queue.pop();
      zsw::Scalar real_err=fabs(jpt_info.second->val_cur_-jpt_info.second->val_exp_);
      if(fabs(real_err-jpt_info.first) > zsw::const_val::eps) { continue; } // have already updated
      PointType pt_type= (jpt_info.second->val_exp_<0) ? zsw::INNER_POINT : zsw::OUTER_POINT;
      VertexInfo vertex_info(-1, pt_type, jpt_info.second->pt_, 0.0);
      std::vector<Chd> chds;
      //tw_->addPointInDelaunay(jpt_info.second->pt_, vertex_info, chds);
      tw_->addPointInDelaunaySafe(jpt_info.second->pt_, vertex_info, chds, &cell_key_set[pre], &cell_key_set[cur]);
      swap(pre,cur);
      for(Chd chd : chds) { updateJptsInCell(chd, &err_queue); }

      // writeTetMesh("/home/wegatron/tmp/refine_debug_"+std::to_string(debug_step++)+".vtk", {ignore_self_in, ignore_bbox, ignore_self_out});
    }
    writeTetMesh("/home/wegatron/tmp/before_checkUpNormalCondition.vtk", {ignore_self_in, ignore_bbox, ignore_self_out});
    std::cout << "[INFO] " << "start checkUpNormalCondition!" << std::endl;
    std::queue<Chd> chds_queue;
    for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) { chds_queue.push(cit); }
    while(!chds_queue.empty()) {
      Chd chd = chds_queue.front(); chds_queue.pop();
      if(!tw_->isTolCell(chd)) { continue; }
      cell_key_set[cur].clear();
      checkUpNormalCondition(chd, chds_queue, &cell_key_set[pre], &cell_key_set[cur]);
      swap(pre,cur);
    }

    // for(const JudgePoint &jpt : jpts_) {
    //   assert(fabs(jpt.val_cur_-jpt.val_exp_)<1);
    // }
  }

  void Approximation::checkUpNormalCondition(Chd chd, std::queue<Chd> &chds_queue,
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
    if(normalCondition(val, tri_pts, scaled_tri_pts, inner_jpts_, outer_jpts_, inner_kdtree_ptr_, outer_kdtree_ptr_)) return;

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
    if(fabs(jpts_[jpt_ind].val_cur_-jpts_[jpt_ind].val_exp_)<zsw::const_val::eps) { return; }
    // add and update
    PointType pt_type=(jpts_[jpt_ind].val_exp_<0) ? zsw::INNER_POINT : zsw::OUTER_POINT;
    VertexInfo vertex_info(-1, pt_type, jpts_[jpt_ind].pt_, 0.0);
    std::vector<Chd> chds;
    tw_->addPointInDelaunaySafe(jpts_[jpt_ind].pt_, vertex_info, chds, cell_key_set_pre, cell_key_set_cur);
    for(Chd chd : chds) { updateJptsInCell(chd, nullptr); }
    for(Chd chd : chds) { chds_queue.push(chd); }
  }

  void Approximation::updateJptsInCell(Chd chd, std::priority_queue<std::pair<zsw::Scalar,JudgePoint*>,
                                       std::vector<std::pair<zsw::Scalar,JudgePoint*>>,
                                       ErrorMaxComparison> *err_queue)
  {
    assert(chd->vertex(0)->info().pt_type_!=zsw::INVALID_POINT && chd->vertex(0)->info().pt_type_!=zsw::INVALID_POINT
           && chd->vertex(0)->info().pt_type_!=zsw::INVALID_POINT && chd->vertex(0)->info().pt_type_!=zsw::INVALID_POINT);
    // calc jpts in bbox
    std::vector<const JudgePoint*> jpts_in_bbox;
    Vhd vhds[4] = {chd->vertex(0), chd->vertex(1), chd->vertex(2), chd->vertex(3)};
    calcJptsInBbox(vhds,4,jpts_in_bbox);
    for(const JudgePoint * jpt : jpts_in_bbox) {
      Eigen::Matrix<zsw::Scalar,4,4> A;
      Eigen::Matrix<zsw::Scalar,4,1> x, b;
      A <<
        chd->vertex(0)->point()[0], chd->vertex(1)->point()[0], chd->vertex(2)->point()[0], chd->vertex(3)->point()[0],
        chd->vertex(0)->point()[1], chd->vertex(1)->point()[1], chd->vertex(2)->point()[1], chd->vertex(3)->point()[1],
        chd->vertex(0)->point()[2], chd->vertex(1)->point()[2], chd->vertex(2)->point()[2], chd->vertex(3)->point()[2],
        1,1,1,1;
      b<< jpt->pt_[0],jpt->pt_[1],jpt->pt_[2],1;
      Eigen::PartialPivLU<Eigen::Matrix<zsw::Scalar,4,4>> pplu; pplu.compute(A);
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

      // if(jpt-&jpts_[0]==15175) {
      //   std::cerr << jpt->pt_.transpose() << std::endl;
      //   std::cerr << cell2str(chd) << std::endl;
      //   std::cerr << val.transpose() << std::endl;
      //   std::cerr << "A:\n" << A << std::endl;
      // }

      JudgePoint *tmp_jpt = const_cast<JudgePoint*>(jpt);
      tmp_jpt->val_cur_=val.dot(x);
      if(err_queue!=nullptr) {
        zsw::Scalar err=fabs(tmp_jpt->val_cur_-tmp_jpt->val_exp_);
        if(err>1) { err_queue->push(std::make_pair(err, tmp_jpt)); }
      }
    }
  }

  // void Approximation::createJudgePoints()
  // {
  //   const TTds &tds=tw_->getTds();
  //   Eigen::Matrix<zsw::Scalar,3,4> bi_tri_points;
  //   Eigen::Matrix<zsw::Scalar,3,4> bo_tri_points;
  //   for(TTds::Cell_iterator cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) {
  //     size_t bi_cnt=0;
  //     size_t bo_cnt=0;
  //     for(size_t i=0; i<4; ++i) {
  //       if(cit->vertex(i)->info().pt_type_==zsw::INNER_POINT) {
  //         bi_tri_points(0,bi_cnt)=cit->vertex(i)->point()[0];
  //         bi_tri_points(1,bi_cnt)=cit->vertex(i)->point()[1];
  //         bi_tri_points(2,bi_cnt)=cit->vertex(i)->point()[2];
  //         ++bi_cnt;
  //       } else if(cit->vertex(i)->info().pt_type_==zsw::OUTER_POINT) {
  //         bo_tri_points(0,bo_cnt)=cit->vertex(i)->point()[0];
  //         bo_tri_points(1,bo_cnt)=cit->vertex(i)->point()[1];
  //         bo_tri_points(2,bo_cnt)=cit->vertex(i)->point()[2];
  //         ++bo_cnt;
  //       }
  //     }
  //     if(bi_cnt==3 && bo_cnt==1) { sampleTriangle(bi_tri_points.block<3,3>(0,0), surf_sample_r_, bi_jpts_);}
  //     else if(bo_cnt==3 && bi_cnt==1) { sampleTriangle(bo_tri_points.block<3,3>(0,0), surf_sample_r_, bo_jpts_); }
  //   }
  //   jpts_ptr_bi_.reset(new zsw::Flann<zsw::Scalar>(bi_jpts_[0].data(), bi_jpts_.size()));
  //   jpts_ptr_bo_.reset(new zsw::Flann<zsw::Scalar>(bo_jpts_[0].data(), bo_jpts_.size()));
  //   jpts_.reserve(bi_jpts_.size()+bo_jpts_.size());
  //   for(const Eigen::Matrix<zsw::Scalar,3,1> &jpt : bi_jpts_) { jpts_.push_back({jpt,-1,-1}); }
  //   for(const Eigen::Matrix<zsw::Scalar,3,1> &jpt : bo_jpts_) { jpts_.push_back({jpt,1,1}); }
  // }

  void Approximation::simpTolerance()
  {
    const TTds &tds=tw_->getTds();
    std::unordered_map<std::string, TTds::Edge> edge_map;
    for(TTds::Edge_iterator eit=tds.edges_begin();
        eit!=tds.edges_end(); ++eit) {
      if(!tw_->isBoundaryEdge(*eit)) { continue;  }
      std::string key_str=edge2key(*eit);
      edge_map.insert(std::make_pair(key_str, *eit));
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
    std::cout << "[INFO] edge size:" << edge_map.size() << std::endl;
    size_t b_c_step=0;
    while(!edge_map.empty()) {
      TTds::Edge e = edge_map.begin()->second; edge_map.erase(edge_map.begin());
      if(!tw_->isBoundaryEdge(e)) { continue; }
      // std::cout << "try collapse edge:";
      // std::cout << e.first->vertex(e.second)->info().index_ << " " <<
      //   e.first->vertex(e.third)->info().index_ << std::endl;
      if(tryCollapseBoundaryEdge(e, edge_map)) {
        if(b_c_step++%50==0) {  std::cout << "[INFO] boundary collapsed:" << b_c_step << std::endl;  }
        // writeTetMesh("/home/wegatron/tmp/sim_tol_"+std::to_string(debug_step++)
        //              +".vtk", {zsw::ignore_bbox, zsw::ignore_self_in, zsw::ignore_self_out});
      }
    }
    std::cout << "[INFO] boundary collapsed total:" << b_c_step << std::endl;
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
        std::string filepath("/home/wegatron/tmp/normal_cond_debug/nd_"+std::to_string(normal_cond_debug_i++)+".vtk");
        normalCondition(val, tri_pts, scaled_tri_pts, inner_jpts_, outer_jpts_,
                        inner_kdtree_ptr_, outer_kdtree_ptr_,
                        true, &filepath);
      }
    }
#endif
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

  bool Approximation::isTolTetsSatisfyNormalCondition(const std::vector<VertexTriple> &bound_tris,
                                       const Eigen::Matrix<zsw::Scalar,3,1> &pt,
                                       const PointType point_type) const
  {
    assert(point_type!=zsw::BBOX_POINT && point_type!=zsw::INVALID_POINT);
    Eigen::Matrix<zsw::Scalar,4,1> val;
    Eigen::Matrix<zsw::Scalar,3,4> tri_pts;
    tri_pts.block<3,1>(0,0)=pt;
    const static Eigen::Matrix<zsw::Scalar,1,4> tmp_v=Eigen::Matrix<zsw::Scalar,1,4>::Ones()*(1-normal_cond_scale_);
    for(VertexTriple vt : bound_tris) {
      if(!isConstructTolCell(point_type, vt.first->info().pt_type_, vt.second->info().pt_type_, vt.third->info().pt_type_, val)) { continue; }
      tri_pts(0,1)=vt.first->point()[0]; tri_pts(1,1)=vt.first->point()[1]; tri_pts(2,1)=vt.first->point()[2];
      tri_pts(0,2)=vt.second->point()[0]; tri_pts(1,2)=vt.second->point()[1]; tri_pts(2,2)=vt.second->point()[2];
      tri_pts(0,3)=vt.third->point()[0]; tri_pts(1,3)=vt.third->point()[1]; tri_pts(2,3)=vt.third->point()[2];
      Eigen::Matrix<zsw::Scalar,3,1> bc=0.25*(
                                              tri_pts.block<3,1>(0,0)+tri_pts.block<3,1>(0,1)+
                                              tri_pts.block<3,1>(0,2)+tri_pts.block<3,1>(0,3));
      Eigen::Matrix<zsw::Scalar,3,4> scaled_tri_pts=normal_cond_scale_*tri_pts+bc*tmp_v;
      if(!normalCondition(val, scaled_tri_pts, tri_pts, inner_jpts_, outer_jpts_, inner_kdtree_ptr_, outer_kdtree_ptr_)) { return false;}
    }
    return true;
  }

  bool Approximation::isSatisfyErrorBound(const std::vector<VertexTriple> &bound_tris,
                                          const std::vector<const JudgePoint*> &jpts_in_bbox,
                                          const Eigen::Matrix<zsw::Scalar,3,1> &merge_pt,
                                          const zsw::Scalar v_pt,
                                          std::vector<VertexUpdateData> &vup,
                                          std::vector<JudgePointUpdateData> * jup) const
  {
    std::vector<bool> is_updated(jpts_in_bbox.size(), false);
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
      for(size_t jpt_i=0; jpt_i<jpts_in_bbox.size(); ++jpt_i) {
        if(is_updated[jpt_i]) { continue; }
        const JudgePoint *jpt=jpts_in_bbox[jpt_i];
        Eigen::Matrix<zsw::Scalar,4,1> b; b.block<3,1>(0,0)=jpt->pt_; b[3]=1;
        Eigen::Matrix<zsw::Scalar,4,1> x=pplu.solve(b);
        if(x[0]<-zsw::const_val::eps || x[1]<-zsw::const_val::eps
           || x[2]<-zsw::const_val::eps || x[3]<-zsw::const_val::eps) { continue; }
        is_updated[jpt_i]=true;
        zsw::Scalar tmp_val = val.dot(x);
        if(fabs(tmp_val-jpt->val_exp_)>1) { return false; }
        if(jup!=nullptr) { jup->push_back({const_cast<JudgePoint*>(jpt), tmp_val}); }
      }
    }
    std::cerr << "TODO calc vup" << std::endl;
    return true;
  }

  bool Approximation::tryCollapseBoundaryEdge(TTds::Edge &e,
                                              std::unordered_map<std::string,TTds::Edge> &edge_map)
  {
    const TTds &tds = tw_->getTds();
    if(!tw_->isSatisfyLinkCondition(e)) {
      // std::cout << "[INFO] link condition failed!" << std::endl;
      // std::cout << e.first->vertex(e.second)->point()[0] << " "
      //           << e.first->vertex(e.second)->point()[1] << " "
      //           << e.first->vertex(e.second)->point()[2] << std::endl;
      // std::cout << e.first->vertex(e.third)->point()[0] << " "
      //           << e.first->vertex(e.third)->point()[1] << " "
      //           << e.first->vertex(e.third)->point()[2] << std::endl;
      // writeAdjcentCells("/home/wegatron/tmp/link_cond_adj_cells.vtk", e);
      // abort();
      return false;
    }
    std::vector<VertexTriple> bound_tris;
    std::vector<Vhd> opposite_vs;
    tw_->calcBoundTris(e, bound_tris, opposite_vs);
    std::vector<const JudgePoint*> jpts_in_bbox;
    calcJptsInBbox(&bound_tris[0].first, 3*bound_tris.size(), jpts_in_bbox);
    // candicate merge points in kernel region
    KernelRegionJudger krj;
    constructKernelRegionJudger(bound_tris, opposite_vs, krj);
    std::vector<const JudgePoint*> candicate_points;
    std::function<bool(zsw::Scalar v)> judge_func;
    if(e.first->vertex(e.second)->info().pt_type_==zsw::INNER_POINT) {
      judge_func=[](zsw::Scalar v){ return v<0; };
    } else { judge_func=[](zsw::Scalar v){ return v>0; }; }
    for(const JudgePoint * jpt_ptr : jpts_in_bbox) {
      if(!judge_func(jpt_ptr->val_exp_)) { continue; }
      if(krj.judge(jpt_ptr->pt_)) { candicate_points.push_back(jpt_ptr); }
    }

#if 0
    // test adj info and jpts selection is right
    writeAdjcentCells("/home/wegatron/tmp/adj_cell.vtk", e);
    writeJudgePoints("/home/wegatron/tmp/jpts_in_bbox.vtk", jpts_in_bbox);
    writeJudgePoints("/home/wegatron/tmp/candicate_points.vtk", candicate_points);
    abort();
#endif

    // sort jpt by error
    sort(candicate_points.begin(), candicate_points.end(), [](const JudgePoint *a, const JudgePoint *b){
        return fabs(a->val_cur_-a->val_exp_) > fabs(b->val_cur_-b->val_exp_);
      });
    const JudgePoint *merge_pt=nullptr;
    std::vector<VertexUpdateData> vup;
    std::vector<JudgePointUpdateData> jup;
    const zsw::Scalar v_pt=(e.first->vertex(e.second)->info().pt_type_==zsw::INNER_POINT) ? -1 : 1;
    std::cout << "candicate pt size:" << candicate_points.size() << std::endl;
    for(const JudgePoint * jpt_ptr : candicate_points) {
      vup.clear(); jup.clear();
      if(isSatisfyErrorBound(bound_tris, jpts_in_bbox, jpt_ptr->pt_, v_pt, vup, &jup)
         && isTolTetsSatisfyNormalCondition(bound_tris, jpt_ptr->pt_, e.first->vertex(e.second)->info().pt_type_))
        { merge_pt=jpt_ptr; break; }
    }
    if(merge_pt==nullptr) { return false; }
    Vhd vhd=e.first->vertex(e.second);
    tw_->collapseEdge(e, vhd, merge_pt->pt_);
    std::for_each(jup.begin(), jup.end(), [](const JudgePointUpdateData &dt){dt.jpt->val_cur_=dt.val_cur_;});
    updateVertex(vup);
    boundaryEdgeBack(vhd, edge_map);
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

  bool Approximation::tryCollapseZeroEdge(TTds::Edge &e,
                                          std::unordered_map<std::string,TTds::Edge> &z_map,
                                          std::unordered_map<std::string,TTds::Edge> *bz_map)
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
    Vhd vhd=e.first->vertex(e.second);
    tw_->collapseEdge(e, vhd, *merge_pt);
    updateVertex(vup);
    zeroEdgeBack(vhd, z_map);
    if(bz_map!=nullptr) { bzEdgeBack(vhd, *bz_map); }
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

  void Approximation::bzEdgeBack(Vhd vhd, std::unordered_map<std::string, TTds::Edge> &edge_map) const
  {
    const TTds &tds=tw_->getTds();
    std::vector<TTds::Edge> edges;
    tds.incident_edges(vhd, std::back_inserter(edges));
    for(const TTds::Edge &e : edges) {
      if(!tw_->isBZEdge(e)) { continue; }
      std::string key_str=edge2key(e);
      edge_map.insert(std::make_pair(key_str, e));
    }
  }

  void Approximation::sampleAdjCells(const TTds::Edge &e, std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &sample_points) const
  {
    Vhd vhds[2]={e.first->vertex(e.second), e.first->vertex(e.third)};
    std::vector<Chd> cells;
    const TTds &tds=tw_->getTds();
    tds.incident_cells(e.first->vertex(e.second), std::back_inserter(cells));
    std::unordered_map<std::string,Chd> cell_map;
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
      sampleTet(v0,v1,v2,v3,tet_sample_r_, sample_points);
    }
  }

  void Approximation::mutuallTessellation()
  {
    const TTds &tds = tw_->getTds();
    std::vector<Vhd> vhds;
    for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
      if(vit->info().pt_type_==zsw::INNER_POINT) { vhds.push_back(vit); }
    }

    for(auto tmp_vhd : vhds) {
      assert(tmp_vhd->info().pt_type_==zsw::INNER_POINT);
      do {
        std::vector<TTds::Edge> edges;
        tds.incident_edges(tmp_vhd, std::back_inserter(edges));
        auto it=find_if(edges.begin(), edges.end(), [](const TTds::Edge &e){
            return e.first->vertex(e.second)->info().pt_type_==zsw::OUTER_POINT ||
            e.first->vertex(e.third)->info().pt_type_==zsw::OUTER_POINT;});
        if(it==edges.end()) { break; }
        TTds::Edge &e = *it;
        Point &pa = e.first->vertex(e.second)->point();
        Point &pb = e.first->vertex(e.third)->point();
        Point pt((pa[0]+pb[0])/2.0, (pa[1]+pb[1])/2.0, (pa[2]+pb[2])/2.0);
        tw_->insertInEdge(e, pt, zsw::ZERO_POINT);
      } while(1);
    }

    // std::cerr << "tds cell size:" << tds.number_of_cells() << std::endl;
    // for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) {
    //   if(isValidCell(cit)) {
    //     updateJptsInCell(cit, nullptr);
    //   }
    // }
    // for(const JudgePoint &jpt : jpts_) {
    //   assert(fabs(jpt.val_cur_-jpt.val_exp_)<1);
    // }
  }

  void Approximation::simpZeroSurface(std::unordered_map<std::string,TTds::Edge> *z_map,
                                      std::unordered_map<std::string,TTds::Edge> *bz_map)
  {
    const bool z_map_create=(z_map==nullptr);
    if(z_map_create) {
      z_map=new std::unordered_map<std::string, TTds::Edge>();
      const TTds &tds = tw_->getTds();
      for(TTds::Edge_iterator eit=tds.edges_begin();
          eit!=tds.edges_end(); ++eit) {
        if(!tw_->isZeroEdge(*eit)) { continue; }
        std::string key_str=edge2key(*eit);
        z_map->insert(std::make_pair(key_str, *eit));
      }
    }
    size_t z_c_step=0;
    while(!z_map->empty()) {
      TTds::Edge e=z_map->begin()->second; z_map->erase(z_map->begin());
      if(!tw_->isZeroEdge(e)) { continue; }
      // std::cout << "[INFO] try collapse zc edge:" << e.first->vertex(e.second)->info().index_
      //           << " " << e.first->vertex(e.third)->info().index_ << std::endl;
      if(tryCollapseZeroEdge(e, *z_map, bz_map)) {
        if(z_c_step++%50==0) { std::cout << "[INFO] zero edge collapsed " << z_c_step << std::endl; }
      }
    }
    if(z_map_create) { delete z_map; }
    std::cout << "[INFO] zero edge collapsed total:" << z_c_step << std::endl;
  }

  void Approximation::writeZeroSurface(const std::string &filepath) const
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  }

  bool Approximation::simpBZEdges(std::unordered_map<std::string,TTds::Edge> *bz_map,
                                  std::unordered_map<std::string,TTds::Edge> *z_map)
  {
    bool bz_map_create=(bz_map==nullptr);
    if(bz_map_create) {
      bz_map=new std::unordered_map<std::string, TTds::Edge>();
      const TTds &tds=tw_->getTds();
      for(auto eit=tds.edges_begin(); eit!=tds.edges_end(); ++eit) {
        if(!tw_->isBZEdge(*eit)) { continue; }
        std::string key_str=edge2key(*eit);
        bz_map->insert(std::make_pair(key_str, *eit));
      }
    }

    size_t zb_c_step=0;
    while(!bz_map->empty()) {
      TTds::Edge e=bz_map->begin()->second; bz_map->erase(bz_map->begin());
      if(!tw_->isBZEdge(e)) { continue; }
      if(tryCollapseBZEdge(e, *bz_map, z_map)) {
        if(zb_c_step++%50==0) { std::cout << "[INFO] zb collapsed " << zb_c_step << std::endl; }
      }
    }
    if(bz_map_create) { delete bz_map; }
    std::cout << "[INFO] zb collapsed total:" << zb_c_step << std::endl;
    return zb_c_step==0;
  }

  bool Approximation::tryCollapseBZEdge(TTds::Edge &e, std::unordered_map<std::string,TTds::Edge> &bz_map,
                                        std::unordered_map<std::string,TTds::Edge> *z_map)
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
    tw_->collapseEdge(e, vhd, *merge_pt);
    updateVertex(vup);
    bzEdgeBack(vhd, bz_map);
    if(z_map!=nullptr) { zeroEdgeBack(vhd, *z_map); }
    return true;
  }

  void Approximation::writeTetMesh(const std::string &filepath,
                                   std::vector<std::function<bool(const TTds::Cell_handle)>> ignore_tet_funcs) const
  {
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    const TTds &tds=tw_->getTds();
    ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
        << tds.number_of_vertices()-1 << " float" << std::endl;

    CGAL::Unique_hash_map<TTds::Vertex_handle,size_t> v_map;
    size_t v_index=0;
    for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
      if(vit->info().pt_type_==zsw::INVALID_POINT) { continue; }
      ofs << *vit << std::endl;
      v_map[vit] = v_index++;
    }

    stringstream ss;
    size_t valid_cells_number=0;
    for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) {
      bool ignore=false;
      for(auto igf : ignore_tet_funcs) {        if(igf(cit)) { ignore=true; break; }      }
      if(ignore || ignore_invalid(cit)) { continue; }
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

  void Approximation::calcJptsInBbox(Vhd *vhd_ptr, const size_t n, std::vector<const JudgePoint*> &jpts_in_bbox) const
  {
    Eigen::Matrix<zsw::Scalar,3,2> bbox;
    calcVerticesBbox(vhd_ptr, n, bbox);
    for(const JudgePoint &jpt : jpts_) {
      if(jpt.pt_[0]<bbox(0,0) || jpt.pt_[1]<bbox(1,0) || jpt.pt_[2]<bbox(2,0) ||
         jpt.pt_[0]>bbox(0,1) || jpt.pt_[1]>bbox(1,1) || jpt.pt_[2]>bbox(2,1)) { continue; }
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

  void Approximation::writeAdjcentCells(const std::string &filepath, const TTds::Edge &e) const
  {
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    const TTds &tds=tw_->getTds();
    ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
        << tds.number_of_vertices() << " float" << std::endl;

    CGAL::Unique_hash_map<TTds::Vertex_handle,size_t> v_map;
    size_t v_index=0;
    for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
      ofs << *vit << std::endl;
      v_map[vit] = v_index++;
    }

    std::unordered_set<TTds::Cell_handle, CGAL::Handle_hash_function> cell_set;
    {
      std::list<TTds::Cell_handle> cells;
      tds.incident_cells(e.first->vertex(e.second), std::back_inserter(cells));
      for(auto chd : cells) { cell_set.insert(chd); }
    }

    {
      std::list<TTds::Cell_handle> cells;
      tds.incident_cells(e.first->vertex(e.third), std::back_inserter(cells));
      for(auto chd : cells) { cell_set.insert(chd); }
    }

    ofs << "CELLS " << cell_set.size() << " " << cell_set.size()*5 << std::endl;
    for(auto cit=cell_set.begin(); cit!=cell_set.end(); ++cit) {
      ofs << "4 " << v_map[(*cit)->vertex(0)] << " " <<
        v_map[(*cit)->vertex(1)] << " " <<
        v_map[(*cit)->vertex(2)] << " " <<
        v_map[(*cit)->vertex(3)] << std::endl;
    }
    ofs << "CELL_TYPES " << cell_set.size() << std::endl;
    for(size_t ci=0; ci<cell_set.size(); ++ci) { ofs << "10" << std::endl; }
    ofs.close();
  }

  void Approximation::writeJudgePoints(const std::string &filepath, const vector<const JudgePoint*> &jpts) const
  {
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
        << jpts.size() << " float" << std::endl;
    for(const JudgePoint *jpt : jpts) {
      ofs << jpt->pt_.transpose() << std::endl;
    }
    ofs << "CELLS " << jpts.size() << " " << jpts.size()*2 << std::endl;
    for(size_t i=0; i<jpts.size(); ++i) {      ofs << "1 " << i << std::endl;    }
    ofs << "CELL_TYPES " << jpts.size() << std::endl;
    for(size_t i=0; i<jpts.size(); ++i) { ofs << "1" << std::endl; }
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

  void writePoints(const std::string &filepath, const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &pts)
  {
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
        << pts.size() << " float" << std::endl;
    for(const Eigen::Matrix<zsw::Scalar,3,1> &pt : pts) {
      ofs << pt.transpose() << std::endl;
    }
    ofs << "CELLS " << pts.size() << " " << pts.size()*2 << std::endl;
    for(size_t i=0; i<pts.size(); ++i) {      ofs << "1 " << i << std::endl;    }
    ofs << "CELL_TYPES " << pts.size() << std::endl;
    for(size_t i=0; i<pts.size(); ++i) { ofs << "1" << std::endl; }
  }

  void Approximation::writeAdjcentCells(const std::string &filepath, const std::vector<Chd> &chds) const
  {
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    const TTds &tds=tw_->getTds();
    ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
        << tds.number_of_vertices() << " float" << std::endl;

    CGAL::Unique_hash_map<TTds::Vertex_handle,size_t> v_map;
    size_t v_index=0;
    for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
      ofs << *vit << std::endl;
      v_map[vit] = v_index++;
    }

    ofs << "CELLS " << chds.size() << " " << chds.size()*5 << std::endl;
    for(auto cit=chds.begin(); cit!=chds.end(); ++cit) {
      ofs << "4 " << v_map[(*cit)->vertex(0)] << " " <<
        v_map[(*cit)->vertex(1)] << " " <<
        v_map[(*cit)->vertex(2)] << " " <<
        v_map[(*cit)->vertex(3)] << std::endl;
    }
    ofs << "CELL_TYPES " << chds.size() << std::endl;
    for(size_t ci=0; ci<chds.size(); ++ci) { ofs << "10" << std::endl; }
    ofs.close();
  }

  void calcCircumcenter(const Eigen::Matrix<zsw::Scalar,3,4> &tri_pts, Eigen::Matrix<zsw::Scalar,3,1> &center)
  {
    Eigen::Matrix<zsw::Scalar,3,3> A;
    A.block<3,1>(0,0)=tri_pts.block<3,1>(0,1)-tri_pts.block<3,1>(0,0);
    A.block<3,1>(0,1)=tri_pts.block<3,1>(0,2)-tri_pts.block<3,1>(0,0);
    A.block<3,1>(0,2)=tri_pts.block<3,1>(0,3)-tri_pts.block<3,1>(0,0);
    Eigen::PartialPivLU<Eigen::Matrix<zsw::Scalar,3,3>> pplu; pplu.compute(A);
    Eigen::Matrix<zsw::Scalar,3,1> b;
    b[0] = 0.5*(A.block<3,1>(0,0)).dot(tri_pts.block<3,1>(0,1)+tri_pts.block<3,1>(0,0));
    b[1] = 0.5*(A.block<3,1>(0,1)).dot(tri_pts.block<3,1>(0,2)+tri_pts.block<3,1>(0,0));
    b[2] = 0.5*(A.block<3,1>(0,2)).dot(tri_pts.block<3,1>(0,3)+tri_pts.block<3,1>(0,0));
    center = pplu.solve(b);
  }
}

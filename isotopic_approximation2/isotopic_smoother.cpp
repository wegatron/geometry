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

namespace zsw {

  zsw::Scalar calcTetHeightType0(const Eigen::Matrix<zsw::Scalar,3,1> &pt0, const Eigen::Matrix<zsw::Scalar,3,1> &pt1,
                                 const Eigen::Matrix<zsw::Scalar,3,1> &pt2, const Eigen::Matrix<zsw::Scalar,3,1> &pt3)
  {
    Eigen::Matrix<zsw::Scalar,3,1> vn=(pt2-pt1).cross(pt3-pt1); vn.normalized();
    return fabs(vn.dot(pt0-pt1));
  }

  zsw::Scalar calcTetHeightType1(const Eigen::Matrix<zsw::Scalar,3,1> &pt0, const Eigen::Matrix<zsw::Scalar,3,1> &pt1,
                                 const Eigen::Matrix<zsw::Scalar,3,1> &pt2, const Eigen::Matrix<zsw::Scalar,3,1> &pt3)
  {
    Eigen::Matrix<zsw::Scalar,3,1> vn=(pt1-pt0).cross(pt3-pt2); vn.normalized();
    return fabs(vn.dot(pt2-pt0));
  }

  // void Approximation::updateAllBoundaryVerticesMaxDis()
  // {
  //   TTds &tds=tw_->getTds();
  //   // set each vertex's max dis to -1
  //   for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
  //     vit->info().pos_ori_<<vit->point()[0], vit->point()[1],vit->point()[2];
  //     vit->info().max_dis_=-1;
  //   }
  //   // for each tet in tds, calc tet's height h
  //   // calc the max jpt val in tet
  //   // max dis = 2*val/h
  //   for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) {
  //     if(!tw_->isTolCell(cit)) { continue; }
  //     zsw::Scalar val = updateJptsInCell(cit, nullptr);
  //     zsw::Scalar h = calcBoundaryTetHeight(cit);
  //     zsw::Scalar max_dis=val*h;
  //     for(size_t i=0; i<4; ++i) {
  //       if(cit->vertex(i)->info().max_dis_<-zsw::const_val::eps ||
  //          cit->vertex(i)->info().max_dis_>max_dis) { cit->vertex(i)->info().max_dis_=max_dis; }
  //     }
  //   }
  // }

  zsw::Scalar calcBoundaryTetHeight(Chd chd)
  {
    size_t inner_cnt=0;
    size_t outer_cnt=0;
    Eigen::Matrix<zsw::Scalar,3,1> inner_pts[4];
    Eigen::Matrix<zsw::Scalar,3,1> outer_pts[4];
    for(size_t i=0; i<4; ++i) {
      if(chd->vertex(i)->info().pt_type_==zsw::INNER_POINT) {
        inner_pts[inner_cnt] << chd->vertex(i)->point()[0], chd->vertex(i)->point()[1], chd->vertex(i)->point()[2];
        ++inner_cnt;
      } else if(chd->vertex(i)->info().pt_type_==zsw::OUTER_POINT) {
        outer_pts[outer_cnt] <<chd->vertex(i)->point()[0], chd->vertex(i)->point()[1], chd->vertex(i)->point()[2];
        ++outer_cnt;
      }
    }
    assert(inner_cnt!=0 && outer_cnt!=0 && inner_cnt+outer_cnt==4);
    zsw::Scalar ret=0;
    if(inner_cnt==1) {  ret=calcTetHeightType0(inner_pts[0], outer_pts[0], outer_pts[1], outer_pts[2]); }
    else if(inner_cnt==2) { ret=calcTetHeightType1(inner_pts[0], inner_pts[1], outer_pts[0], outer_pts[1]); }
    else { ret=calcTetHeightType0(outer_pts[0], inner_pts[0], inner_pts[1], inner_pts[2]); }
    return ret;
  }

  void Approximation::smoothBoundary()
  {
    TTds &tds=tw_->getTds();
    PointType pt_type[2] = {zsw::INNER_POINT, zsw::OUTER_POINT};
    for(size_t up_i=0; up_i<2; ++up_i) {
      updateAllBoundaryVerticesMaxDis();
      // boundary smooth
      for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
        if(vit->info().pt_type_!=pt_type[up_i]) { continue; }
        // find one ring pts
        std::vector<Eigen::Matrix<zsw::Scalar,3,1>> ring_pts;
        calcBoundaryOneRing(vit, ring_pts);
        Eigen::Matrix<zsw::Scalar,3,1> pos;
        pos<<vit->point()[0], vit->point()[1], vit->point()[2];
        // calc smooth n limit it in the sphere r of pos_ori
        laplaceSmooth(pos, ring_pts, vit->info().pos_ori_, 1, vit->info().max_dis_);
        Point npt(pos[0],pos[1],pos[2]);
        vit->set_point(npt);
      }
    }
  }

  void Approximation::calcBoundaryOneRing(Vhd vhd, std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &ring_pts) const
  {
    const TTds &tds=tw_->getTds();
    std::vector<TTds::Edge> edges;
    tds.incident_edges(vhd, std::back_inserter(edges));
    for(TTds::Edge &e : edges) {
      if(!tw_->isBoundaryEdge(e)) { continue; }
      Eigen::Matrix<zsw::Scalar,3,1> tmp_pt;
      if(e.first->vertex(e.second)!=vhd) {
        tmp_pt << e.first->vertex(e.second)->point()[0],
          e.first->vertex(e.second)->point()[1],
          e.first->vertex(e.second)->point()[2];
      } else {
        tmp_pt << e.first->vertex(e.third)->point()[0],
          e.first->vertex(e.third)->point()[1],
          e.first->vertex(e.third)->point()[2];
      }
      ring_pts.push_back(tmp_pt);
    }
  }

  // void Approximation::updateAllZeroVerticesMaxDis()
  // {
  //   TTds &tds=tw_->getTds();
  //   for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
  //     if(vit->info().pt_type_!=zsw::ZERO_POINT) { continue; }
  //     vit->info().pos_ori_<< vit->point()[0], vit->point()[1], vit->point()[2];
  //     vit->info().max_dis_=-1;
  //   }
  //   // calc incident cell of zero point
  //   // calc tet's height
  //   // ids = val*height
  //   for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) {
  //     if(cit->vertex(0)->info().pt_type_!=zsw::ZERO_POINT && cit->vertex(1)->info().pt_type_!=zsw::ZERO_POINT &&
  //        cit->vertex(2)->info().pt_type_!=zsw::ZERO_POINT && cit->vertex(3)->info().pt_type_!=zsw::ZERO_POINT) { continue; }
  //     zsw::Scalar h=calcZeroTetHeight(cit);
  //     zsw::Scalar val=updateJptsInCell(cit, nullptr);
  //     zsw::Scalar max_dis=val*h;
  //     // std::cerr << "h:" << h << std::endl;
  //     // std::cerr << "val:" << val << std::endl;
  //     // std::cerr << "max_dis:" << max_dis << std::endl;
  //     for(size_t i=0; i<4; ++i) {
  //       if(cit->vertex(i)->info().pt_type_==zsw::ZERO_POINT &&
  //          (cit->vertex(i)->info().max_dis_<-zsw::const_val::eps ||
  //           cit->vertex(i)->info().max_dis_>max_dis)) {
  //         cit->vertex(i)->info().max_dis_=max_dis;
  //         // std::cerr << "updated max_dis:" << max_dis << std::endl;
  //       }
  //     }
  //   }
  // }

  zsw::Scalar calcZeroTetHeight(Chd chd)
  {
    size_t zero_cnt=0;
    size_t bound_cnt=0;
    Eigen::Matrix<zsw::Scalar,3,1> zero_pts[4];
    Eigen::Matrix<zsw::Scalar,3,1> bound_pts[4];
    for(size_t i=0; i<4; ++i) {
      if(chd->vertex(i)->info().pt_type_==zsw::ZERO_POINT) {
        zero_pts[zero_cnt] << chd->vertex(i)->point()[0], chd->vertex(i)->point()[1], chd->vertex(i)->point()[2];
        ++zero_cnt;
      } else  {
        bound_pts[bound_cnt] <<chd->vertex(i)->point()[0], chd->vertex(i)->point()[1], chd->vertex(i)->point()[2];
        ++bound_cnt;
      }
    }
    zsw::Scalar ret=0;
    if(zero_cnt==1) {  ret=calcTetHeightType0(zero_pts[0], bound_pts[0], bound_pts[1], bound_pts[2]); }
    else if(zero_cnt==2) { ret=calcTetHeightType1(zero_pts[0], zero_pts[1], bound_pts[0], bound_pts[1]); }
    else { ret=calcTetHeightType0(zero_pts[0], zero_pts[1], zero_pts[2], bound_pts[0]); }
    return ret;
  }

  // void Approximation::laplaceSmoothZeroSurface()
  // {
  //   updateAllZeroVerticesMaxDis();
  //   TTds &tds=tw_->getTds();
  //   const size_t N=20;
  //   std::vector<std::pair<Vhd, Eigen::Matrix<zsw::Scalar,3,1>>> vpts;
  //   for(size_t times=0; times<N; ++times) {
  //     for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
  //       if(vit->info().pt_type_!=zsw::ZERO_POINT) { continue; }
  //       std::vector<Eigen::Matrix<zsw::Scalar,3,1>> ring_pts;
  //       calcZeroSurfaceOneRing(vit, ring_pts);
  //       Eigen::Matrix<zsw::Scalar,3,1> pos;
  //       pos<<vit->point()[0], vit->point()[1], vit->point()[2];
  //       // laplaceSmooth(pos, ring_pts, vit->info().pos_ori_, N, vit->info().max_dis_);
  //       laplaceSmooth(pos, ring_pts, vit->info().pos_ori_, N, 0.1);
  //       vpts.push_back(std::make_pair(vit, pos));
  //     }
  //     for(std::pair<Vhd, Eigen::Matrix<zsw::Scalar,3,1>> &vpt : vpts) {
  //       Point pt(vpt.second[0], vpt.second[1], vpt.second[2]);
  //       vpt.first->set_point(pt);
  //     }
  //     vpts.clear();
  //   }
  // }

  void Approximation::bilateralSmoothZeroSurface()
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  }

  void Approximation::calcZeroSurfaceOneRing(Vhd vhd, std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &ring_pts) const
  {
    std::vector<Vhd> vhds;
    const TTds &tds=tw_->getTds();
    tds.adjacent_vertices(vhd, std::back_inserter(vhds));
    for(Vhd tmp_vhd : vhds) {
      if(tmp_vhd->info().pt_type_!=zsw::ZERO_POINT) { continue; }
      Eigen::Matrix<zsw::Scalar,3,1> pt;
      pt<<tmp_vhd->point()[0], tmp_vhd->point()[1], tmp_vhd->point()[2];
      ring_pts.push_back(pt);
    }
  }
}

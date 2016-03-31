#include <iostream>
#include "isotopic_debug.h"
using namespace std;

namespace zsw{
  bool isZeroTetExist(const TTds &tds)
  {
    for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) {
      if(cit->vertex(0)->info().pt_type_==zsw::ZERO_POINT &&
         cit->vertex(1)->info().pt_type_==zsw::ZERO_POINT &&
         cit->vertex(2)->info().pt_type_==zsw::ZERO_POINT &&
         cit->vertex(3)->info().pt_type_==zsw::ZERO_POINT) {
        return true;
      }
    }
    return false;
  }

  void testKdtreeFunc(const KdTreeWarper &kdtree, const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &jpts)
  {
    Eigen::Matrix<zsw::Scalar,3,1> q_pt=Eigen::Matrix<zsw::Scalar,3,1>::Random();
    std::vector<size_t> indices;
    std::vector<zsw::Scalar> dist;
    kdtree.queryNearest(q_pt, indices, dist);
    size_t real_min_ind=indices[0];
    zsw::Scalar real_min_dis=dist[0];
    for(const Eigen::Matrix<zsw::Scalar,3,1> &pt : jpts) {
      if((pt-q_pt).squaredNorm() < real_min_dis) {
        std::cerr << "kdtree check failed!!!!" << std::endl;
        abort();
      }
    }
  }

  void Approximation::testKdtree() const
  {
    for(size_t i=0; i<1000; ++i) {
      testKdtreeFunc(outer_kdtree_, outer_jpts_);
    }
    for(size_t i=0; i<1000; ++i) {
      testKdtreeFunc(inner_kdtree_, inner_jpts_);
    }
    std::cerr << "kdtree test_pass!!!" << std::endl;
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
                          inner_kdtree_, outer_kdtree_)) {
        if(!normalCondition(val, scaled_tri_pts, tri_pts, inner_jpts_, outer_jpts_,
                            inner_kdtree_, outer_kdtree_,
                            true, &filepath)) {
          std::cerr << "normal cond failed!!!" << std::endl;
          std::cerr << "normal_cond_debug_i=" << normal_cond_debug_i-1 << std::endl;
          return false;
        }
      }
    }
    return true;
  }

  bool checkZeroPlane(const Plane &plane, Chd chd)
  {
    Eigen::Matrix<zsw::Scalar,3,1> pt;
    zsw::Scalar dis[4];
    for(size_t i=0; i<4; ++i) {
      pt[0] = chd->vertex(i)->point()[0];
      pt[1] = chd->vertex(i)->point()[1];
      pt[2] = chd->vertex(i)->point()[2];
      dis[i] = (pt-plane.v0_).dot(plane.normal_);
      if(chd->vertex(i)->info().pt_type_==zsw::OUTER_POINT){
        if(dis[i]<0) { return false; }
      } else if(chd->vertex(i)->info().pt_type_==zsw::INNER_POINT) {
        if(dis[i]>0) { return false; }
        dis[i]=-dis[i];
      }
    }
    return fabs(dis[1]-dis[0])<zsw::const_val::eps && fabs(dis[2]-dis[0])<zsw::const_val::eps && fabs(dis[3]-dis[0])<zsw::const_val::eps;
  }
}

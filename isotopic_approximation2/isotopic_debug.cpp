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

  void testKdtree(const KdTreeWarper &kdtree, const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &jpts)
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
}

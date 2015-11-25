#ifndef TRIANGULATION_FIX_H
#define TRIANGULATION_FIX_H

#include <Eigen/Dense>
#include <zswlib/config.h>
#include <zswlib/const_val.h>
#include "cgal_common.h"

namespace zsw{
  void removeSliverTet(Delaunay &delaunay);

  bool isSliverTet(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> &v1,
                   const Eigen::Matrix<zsw::Scalar,3,1> &v2, const Eigen::Matrix<zsw::Scalar,3,1> &v3,
                   std::pair<size_t,size_t> &e0, std::pair<size_t,size_t> &e1)
  {
    Eigen::Matrix<zsw::Scalar,3,3> A;
    A.block<3,1>(0,0) = v1-v0; A.block<3,1>(0,0) = v2-v0;
    A.block<3,1>(0,0) = v3-v0;
    if(fabs(A.determinant()) > zsw::const_val::eps) { return false; }
    Eigen::Matrix<zsw::Scalar,3,1> vn = (v1-v0).cross(v2-v0); vn.normalize();
    Eigen::Matrix<zsw::Scalar,3,1> vm = (v1-v0).cross(vn); vm.normalize();
    if((vm.cross(v3-v0)).dot(vm.cross(v2-v0)) < 0) {
      e0 = std::pair<size_t,size_t>(0,1);
      e1 = std::pair<size_t,size_t>(2,3);
    } else {
      e0 = std::pair<size_t,size_t>(0,2);
      e1 = std::pair<size_t,size_t>(1,3);
    }
    return true;
  }
}

#endif /* TRIANGULATION_FIX_H */

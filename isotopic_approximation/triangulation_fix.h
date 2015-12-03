#ifndef TRIANGULATION_FIX_H
#define TRIANGULATION_FIX_H

#include <Eigen/Dense>
#include <zswlib/config.h>
#include <zswlib/const_val.h>
#include "cgal_common.h"
#include "basic_data_structure.h"

namespace zsw{
  void removeSliverTet(Delaunay &delaunay, const std::vector<zsw::Vertex> &vertices);

  bool isSliverTet(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> &v1,
                   const Eigen::Matrix<zsw::Scalar,3,1> &v2, const Eigen::Matrix<zsw::Scalar,3,1> &v3,
                   std::pair<size_t,size_t> &e0, std::pair<size_t,size_t> &e1);

  void haveSliverTet(Delaunay &delaunay, const std::vector<zsw::Vertex> &vertices);

  zsw::Scalar calcTetQuality(const Eigen::Matrix<zsw::Scalar,3,1> &v0,
                             const Eigen::Matrix<zsw::Scalar,3,1> &v1,
                             const Eigen::Matrix<zsw::Scalar,3,1> &v2,
                             const Eigen::Matrix<zsw::Scalar,3,1> &v3);

}

#endif /* TRIANGULATION_FIX_H */

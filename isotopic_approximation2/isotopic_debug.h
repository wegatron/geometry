#ifndef ISOTOPIC_DEBUG_H
#define ISOTOPIC_DEBUG_H

#include "isotopic_approximation.h"

namespace zsw
{
  bool isZeroTetExist(const TTds &tds);
  void testKdtreeFunc(const KdTreeWarper &kdtree, const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &jpts);
  bool checkZeroPlane(const Plane &plane, Chd chd);
}

#endif /* ISOTOPIC_DEBUG_H */

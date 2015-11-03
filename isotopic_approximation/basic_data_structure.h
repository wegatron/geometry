#ifndef BASIC_DATA_STRUCTURE_H
#define BASIC_DATA_STRUCTURE_H

#include <list>
#include <vector>
#include <Eigen/Dense>

namespace zsw{
  enum PointType
  {
    BBOX_POINT=1,
    ZERO_POINT=2,
    OUTER_POINT=4,
    INNER_POINT=8
  };

  struct JudgePoint
  {
    const Eigen::Matrix<zsw::Scalar,3,1> pt_;
    const zsw::Scalar val_exp_;
    zsw::Scalar val_cur_;
  };

  struct Vertex
  {
    const PointType pt_type_;
    Eigen::Matrix<zsw::Scalar,3,1> pt_;
    std::vector<size_t> tet_ids_;
    std::vector<size_t> edge_ids_;
  };

  struct Edge
  {
    size_t vid_[2];
  };

  struct Tet
  {
    bool valid_;
    size_t vid_[4];
    std::list<JudgePoint> jpts_;
  };
}

#endif /* BASIC_DATA_STRUCTURE_H */

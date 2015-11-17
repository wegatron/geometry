#ifndef BASIC_DATA_STRUCTURE_H
#define BASIC_DATA_STRUCTURE_H

#include <list>
#include <vector>
#include <Eigen/Dense>
#include <zswlib/config.h>

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
    Eigen::Matrix<zsw::Scalar,3,1> pt_;
    zsw::Scalar val_exp_;
    zsw::Scalar val_cur_;
  };

  struct Vertex
  {
    bool valid_;
    PointType pt_type_;
    Eigen::Matrix<zsw::Scalar,3,1> pt_;
    std::list<size_t> tet_ids_;
    std::list<size_t> edge_ids_;
  };

  struct Edge
  {
    bool valid_;
    size_t vid_[2];
  };

  struct Tet
  {
    bool valid_;
    size_t vid_[4];
  };
}

#endif /* BASIC_DATA_STRUCTURE_H */

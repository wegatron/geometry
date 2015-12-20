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

  enum EdgeType
  {
    ZERO_EDGE=1,
    OUTER_EDGE=2,
    INNER_EDGE=4,
    OTHETR_EDGE=8
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
    // std::list<JudgePoint> jpts_;
  };

  class BoundSphere
  {
  public:
    BoundSphere(const std::string &filepath, const zsw::Scalar scale, const Eigen::Matrix<zsw::Scalar,3,1> &transform);
    const std::vector<Eigen::Matrix<zsw::Scalar,3,1>>& getVertices() { return vertices_; }
  private:
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> vertices_;
  };
}

#endif /* BASIC_DATA_STRUCTURE_H */

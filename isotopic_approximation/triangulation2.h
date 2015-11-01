/// \file triangulation2.h
/// \brief 3d triangulation structure.
///
/// we build this 3d triangulation from bi and bo points,
/// simplify it then we get the result.
///
/// \author wegatron
/// \email wegatron@hotmail.com
/// \version 1.0
/// \date 2015-11-01

#ifndef TRIANGULATION2_H
#define TRIANGULATION2_H

#include <vector>
#include <Eigen/Dense>
#include <zswlib/config.h>
#include "cgal_common.h"

namespace zsw
{

  /// \brief kernel region judger.
  ///
  ///  judge wether the point is in kernel region.
  ///
  class KernelRegionJudger final
  {
  public:
    KernelRegionJudger() {}
    void addConstraint(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> &v1,
                       const Eigen::Matrix<zsw::Scalar,3,1> &v2, const Eigen::Matrix<zsw::Scalar,3,1> &vr);
    bool judge(const Eigen::Matrix<zsw::Scalar,3,1> &pt);
  private:
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> vec_v0;
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> vec_vn;
  };

  class Triangulation final
  {
  public:
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
      size_t vid_[4];
      std::list<JudgePoint> jpts_;
    };

    Triangulation(const zsw::Scalar r,
                  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bo_points,
                  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bi_points);
    void simpTolerance();
    void mutualTessellation();
    void writeTetMesh(const std::string &filepath, size_t mask) const;
    void writeSurface(const std::string &filepath, PointType pt_tyte) const;
  private:
    void initTets(Delaunay &delaunay);
    bool testCollapse(Edge &e, const Eigen::Matrix<zsw::Scalar,3,1> &pt, std::list<JudgePoint> jpts) const;
    void edgeCollapse(Edge &e, const Eigen::Matrix<zsw::Scalar,3,1> &pt, std::list<JudgePoint> jpts);
    std::vector<Edge> edges_;
    std::vector<Vertex> vertices_;
    std::vector<Tet> tets_;
  };
}

#endif /* TRIANGULATION2_H */

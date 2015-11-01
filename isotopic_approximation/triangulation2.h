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

namespace zsw
{

  /// \brief kernel region judger.
  ///
  ///  judge wether the point is in kernel region.
  ///
  class KernelRegionJudger
  {
  public:
    KernelRegionJudger() {}
    void addConstraint(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> &v1,
                       const Eigen::Matrix<zsw::Scalar,3,1> &v2, const Eigen::Matrix<zsw::Scalar,3,1> &vr);
    bool judge(const Point &pt);
  private:
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> vec_v0;
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> vec_vn;
  };

  class Triangulation final
  {
  public:
    enum PointType
    {
      BBOX=1,
      ZERO_POINT=2,
      OUTER_POINT=4,
      INNER_POINT=8
    };
    struct JudgePoint
    {
      const Eigen::Matrix<zsw::Scalar,3,1> pt_;
      const zsw::Scalar val_exp_;
      const zsw::Scalar val_cur_;
    };
    struct Vertex
    {
      const PointType pt_type_;
      Point pt_;
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
    Triangulation(const zsw::Scalar r, std::vector<Point> &bo_points, std::vector<Point> &bi_points);
    void simpTolerance();
    void mutualTessellation();
    void writeTetMesh(const string &filepath, size_t mask);
    void writeSurface(const string &filepath, PointType pt_tyte);
  private:
    void edgeCollapse(size_t eid, const Point &pt);
    std::vector<Edge> edges_;
    std::vector<Vertex> vertices_;
    std::vector<Tet> tets_;
  };
}

#endif /* TRIANGULATION2_H */

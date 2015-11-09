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
#include <map>
#include <Eigen/Dense>
#include <zswlib/config.h>
#include "cgal_common.h"
#include "basic_data_structure.h"

#define ZSW_DEBUG

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

  typedef bool(*PairCompFunc)(const std::pair<size_t,size_t>&a, const std::pair<size_t,size_t>&b);

  class Triangulation final
  {
  public:
#ifdef ZSW_DEBUG
    Triangulation() {}
    std::vector<Vertex>& getVertices() { return vertices_; }
    std::vector<Tet>& getTets()  { return tets_; }
    std::vector<Edge>& getEdges() { return edges_; }
    bool testLinkCondition(const Edge &e) const { return linkCondition(e); }
#endif
    Triangulation(const zsw::Scalar r,
                  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bo_points,
                  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bi_points);
    void simpTolerance();
    void mutualTessellation();
    void writeTetMesh(const std::string &filepath, size_t mask) const;
    void writeSurface(const std::string &filepath, PointType pt_tyte) const;


    /// \brief a funtion for debug or check.
    ///
    /// write the single tet in the triangulation.
    ///
    ///
    /// \param
    void writeTet(const std::string &filepath, const size_t tet_id) const;

    const std::vector<Vertex>& getVertices() const { return vertices_; }
    const std::vector<Tet>& getTets() const { return tets_; }
    const std::vector<Edge>& getEdges() const { return edges_; }
  private:

    /// \brief init the triangulation's tets from 3d delaunay triangulation.
    ///
    /// including fill the basic tets data, and bi, bo's surface sampling
    /// and point sampling in the tets which is for zero point set edge collapse.
    ///
    /// \param r the sample radius in bi and bo surface
    /// \param Delaunay cgal's delaunay triangulation
    void init(const zsw::Scalar r, Delaunay &delaunay);


    bool linkCondition(const Edge &e) const;

    /// \brief test if the edge can collapse to this point
    ///
    /// check if the judge points' error is satisfied,
    /// in other words check whether the classifcation of S is keeped.
    ///
    /// \param e input edge
    /// \param pt the point this edge collapse to
    /// \param jpts the judge points
    /// \return true if S is keeped or false otherwise
    bool testCollapse(const Edge &e, const PointType pt_type, const Eigen::Matrix<zsw::Scalar,3,1> &pt,
                      const std::list<Eigen::Matrix<size_t,3,1>> &bound_tris,
                      const std::list<JudgePoint> &jpts) const;

    void edgeCollapse(Edge &e, const PointType pt_type,
                      const std::list<Eigen::Matrix<size_t,3,1>> &bound_tris,
                      const Eigen::Matrix<zsw::Scalar,3,1> &pt, std::list<JudgePoint> jpts);

    void addZeroPoints(std::map<std::pair<size_t,size_t>, size_t, PairCompFunc> &ev_map);

    void tessellation3v1(const size_t vo_0, const size_t vo_1, const size_t vo_2, const size_t vi_0,
                         Tet &tet, std::map<std::pair<size_t,size_t>,
                         size_t, PairCompFunc> &ev_map);

    void tessellation2v2(const size_t vo_0, const size_t vo_1, const size_t vi_0, const size_t vi_1,
                         Tet &tet, std::map<std::pair<size_t,size_t>, size_t,
                         PairCompFunc> &ev_map );

    void tessellation1v3(const size_t vo_0, const size_t vi_0, const size_t vi_1, const size_t vi_2,
                         Tet &tet, std::map<std::pair<size_t,size_t>, size_t,
                         PairCompFunc> &ev_map);

    std::vector<Edge> edges_;
    std::vector<Vertex> vertices_;
    std::vector<Tet> tets_;
  };
}

#endif /* TRIANGULATION2_H */

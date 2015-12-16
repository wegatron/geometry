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
#include <queue>
#include <unordered_set>
#include <Eigen/Dense>
#include <zswlib/config.h>
#include <zswlib/mesh/zsw_flann2.h>
#include "cgal_common.h"
#include "basic_data_structure.h"
#include "constraint.h"

//#define ZSW_DEBUG

#define NORMAL_CONT_TOL 0.9

namespace zsw
{
  typedef bool(*PairCompFunc)(const std::pair<size_t,size_t>&a, const std::pair<size_t,size_t>&b);

  class Triangulation final
  {
  public:

    // functions for debug start{
    std::vector<Vertex>& getVertices() { return vertices_; }
    std::vector<Tet>& getTets()  { return tets_; }
    std::vector<Edge>& getEdges() { return edges_; }
    bool testLinkCondition(const Edge &e) const { return linkCondition(e); }
    void testCollapseDebug(const size_t vid0, const size_t vid1);
    void checkTetEdgeExist(const size_t n0, const size_t n1, const size_t n2, const size_t n3);
    // functions for debug end }

    Triangulation() {}

    /// \brief construct
    ///
    /// A construct a triangulation with bo and bi points, and sample
    /// with r radius to get the judge points, The S in paper.
    ///
    /// \param r, sample triangle with r, as judge points. the S in paper.
    /// \param bo_points points in outer boundary
    /// \param bi_points points in inner boundary
    size_t construct(const zsw::Scalar r,
                  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bo_points,
                  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bi_points);

    /// \brief check triangulation is valid(we can use the triangulation to get a result)
    ///
    /// check there is no edge link inner vertex and bbox vertex, if there exist,
    /// we cant generate no zero surface. Also check no valid tet include invalid vertex
    ///
    /// \note check error using error code, not use bool as return value
    /// \warning
    /// \return 0 if the triangulation is good, error line number otherwise.
    size_t isGood() const;

    void simpTolerance();
    void simpZeroSurface();
    void mutualTessellation();
    void writeAllJpts(const std::string &filepath) const;
    void writeTetMeshAdjVs(const std::string &filepath, const std::vector<size_t> &vids) const;
    void writeTetMesh(const std::string &filepath, std::vector<std::function<bool(const Tet &tet)>> ignore_tet_funcs) const;
    void writeSurface(const std::string &filepath, PointType pt_type) const;
    void writeSurface2(const std::string &filepath, PointType pt_type) const;

    /// \brief a funtion for debug or check.
    ///
    /// write the single tet in the triangulation.
    ///
    /// \param
    void writeTet(const std::string &filepath, const size_t tet_id) const;

    void writeBoundTris(const std::string &filepath, const size_t vid0, const size_t vid1);

    bool ignoreWithPtType(const Tet &tet, PointType pt_type);

    bool ignoreOnlyWithPtType(const Tet &tet, PointType pt_type);

    bool ignoreNotWithPtType(const Tet &tet, PointType pt_type);

    void tryCollapseBoundaryEdge(const size_t e_id,
                                 std::set<size_t> &eids_set);

    void tryCollapseZeroEdge(const size_t e_id,
                                 std::set<size_t> &eids_set);

    void reAssignJpts();

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

    void initQem(const size_t vid0, const size_t vid1, const size_t vid2);

    bool linkCondition(const Edge &e) const;

    void initNormalCond(NormalConditionJudger &ncj, const Edge &e) const;

    bool isKeepJpts(const zsw::Scalar pt_val, const Eigen::Matrix<zsw::Scalar,3,1> &pt,
                    const std::vector<Eigen::Matrix<size_t,3,1>> &bound_tris,
                    std::map<size_t,zsw::Scalar> &jpts_update) const;

    void edgeCollapse(const std::vector<size_t> &tet_ids,
                      const std::vector<Eigen::Matrix<size_t,3,1>> &bound_tris,
                      const Eigen::Matrix<zsw::Scalar,3,1> &pt,
                      const zsw::PointType pt_type,
                      const std::map<size_t, zsw::Scalar> &jpts_update,
                      Edge &e,
                      std::function<void(const size_t e_id)> eb_func);

    void addZeroPoints(std::map<std::pair<size_t,size_t>, size_t, PairCompFunc> &ev_map);

    void invalidEdge(const size_t e_id);

    void invalidTet(Tet &tet);

    zsw::Scalar tet_sample_r_;
    std::vector<Edge> edges_;
    std::vector<Vertex> vertices_;
    std::vector<Tet> tets_;
    std::vector<JudgePoint> jpts_;
    std::vector<JudgePoint> bi_jpts_;
    std::vector<JudgePoint> bo_jpts_;
    std::shared_ptr<zsw::Flann2<zsw::Scalar,2>> jpts_ptr_bi_;
    std::shared_ptr<zsw::Flann2<zsw::Scalar,2>> jpts_ptr_bo_;
  };
}

#endif /* TRIANGULATION2_H */

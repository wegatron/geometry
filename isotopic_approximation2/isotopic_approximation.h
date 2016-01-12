#ifndef ISOTOPIC_APPROXIMATION_H
#define ISOTOPIC_APPROXIMATION_H

#include <unordered_map>
#include <zswlib/config.h>
#include <zswlib/mesh/zsw_flann.h>
#include "basic_data_structure.h"
#include "constraint.h"

namespace zsw {

  struct JudgePointUpdateData
  {
    size_t index_;
    zsw::Scalar val_cur_;
  };

  struct VertexUpdateData
  {
    Vhd vhd_;
    zsw::Scalar max_ids_;
    Eigen::Matrix<zsw::Scalar,3,1> pos_ori_;
  };

  class Approximation final
  {
  public:
    Approximation() {}
    /* void init(const zsw::Scalar &surf_sample_r, const zsw::Scalar &tet_sample_r, */
    /*           const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bi_vertices, */
    /*           const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bo_vertices); */
    void init(const zsw::Scalar err_epsilon,
              const zsw::Scalar tri_sample_r,
              const zsw::Scalar tet_sample_r,
              std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &inner_jpts,
              std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &outer_jpts,
              std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bs_jpts);
    void simpTolerance();
    void mutuallTessellation();
    void simpZeroSurface();

    void writeZeroSurface(const std::string &filepath) const;
    void writeTetMesh(const std::string &filepath,
                      std::vector<std::function<bool(const TTds::Cell_handle)>> ignore_tet_funcs) const;
    void writeAdjcentCells(const std::string &filepath, const TTds::Edge &e) const;
    void writeJudgePoints(const std::string &filepath) const;
    void writeJudgePoints(const std::string &filepath, const std::vector<const JudgePoint*> &jpts) const;

    bool isSatisfyErrorBound(std::vector<VertexTriple> &bound_tris,
                             const std::vector<const JudgePoint*> &jpts_in_bbox,
                             const Eigen::Matrix<zsw::Scalar,3,1> &merge_pt,
                             std::vector<VertexUpdateData> &vup,
                             std::vector<JudgePointUpdateData> * jup=nullptr) const
    {
      std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
    }

    void constructKernelRegionJudger(const std::vector<VertexTriple> &bound_tris,
                                     std::vector<Vhd> &opposite_vs, KernelRegionJudger &krj) const;

    void tryCollapseBoundaryEdge(TTds::Edge &e,
                                 std::unordered_map<std::string,TTds::Edge> &edge_map);

    void tryCollapseZeroEdge(TTds::Edge &e,
                             std::unordered_map<std::string,TTds::Edge> &edge_map);

    //void createJudgePoints();

    void updateVertex(const std::vector<VertexUpdateData> &vup)
    {
      std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
    }

    void updateJudgePoint(const std::vector<JudgePointUpdateData> &jup)
    {
      std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
    }

    void boundaryEdgeBack(Vhd vhd, std::unordered_map<std::string, TTds::Edge> &edge_map) const
    {
      std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
    }

    void calcJptsInBbox(std::vector<VertexTriple> &bound_tris, std::vector<const JudgePoint*> &jpts_in_bbox) const;

    void sampleIncidentCells(const TTds::Edge &e, std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &sample_points)
    {
      std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
    }

    void zeroEdgeBack(Vhd vhd, std::unordered_map<std::string,TTds::Edge> &edge_map)
    {
      std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
    }

    // data access

    // for debug
    void testTdsValid();
    void testCollapse();
  private:
    std::vector<JudgePoint> jpts_;
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bi_jpts_;
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bo_jpts_;
    std::shared_ptr<zsw::Flann<zsw::Scalar>> jpts_ptr_bi_;
    std::shared_ptr<zsw::Flann<zsw::Scalar>> jpts_ptr_bo_;
    std::shared_ptr<TriangulationWapper> tw_;
    zsw::Scalar err_epsilon_;
    zsw::Scalar tri_sample_r_;
    zsw::Scalar tet_sample_r_;
  };

  void calcVertexTripleBBox(const std::vector<VertexTriple> &bound_tris, Eigen::Matrix<zsw::Scalar,3,2> &bbox);
  void writePoints(const std::string &filepath, const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &pts);
}
#endif /* ISOTOPIC_APPROXIMATION_H */

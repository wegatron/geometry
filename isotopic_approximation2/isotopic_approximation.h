#ifndef ISOTOPIC_APPROXIMATION_H
#define ISOTOPIC_APPROXIMATION_H

#include <unordered_map>
#include <queue>
#include <zswlib/config.h>
#include <zswlib/mesh/zsw_flann.h>
#include "basic_data_structure.h"
#include "constraint.h"

namespace zsw {

  struct JudgePointUpdateData
  {
    JudgePoint *jpt;
    zsw::Scalar val_cur_;
  };

  struct VertexUpdateData
  {
    Vhd vhd_;
    zsw::Scalar max_ids_;
    Eigen::Matrix<zsw::Scalar,3,1> pos_ori_;
  };

  class ErrorMaxComparison
  {
  public:
    bool operator()(const std::pair<zsw::Scalar,JudgePoint*> &lv,
                    const std::pair<zsw::Scalar,JudgePoint*> &rv){
      return lv.first>rv.first;
    }
  };

  class Approximation final
  {
  public:
    Approximation() { normal_cond_scale_=0.7; }
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
    void simpZeroSurface(std::unordered_map<std::string,TTds::Edge> *z_map=nullptr,
                         std::unordered_map<std::string,TTds::Edge> *bz_map=nullptr);
    bool simpBZEdges(std::unordered_map<std::string,TTds::Edge> *bz_map=nullptr,
                     std::unordered_map<std::string,TTds::Edge> *z_map=nullptr);
    void simpAllEdges() {
      std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
    }

    void writeZeroSurface(const std::string &filepath) const;
    void writeTetMesh(const std::string &filepath,
                      std::vector<std::function<bool(const TTds::Cell_handle)>> ignore_tet_funcs) const;
    void writeAdjcentCells(const std::string &filepath, const TTds::Edge &e) const;
    void writeAdjcentCells(const std::string &filepath, const std::vector<Chd> &chds) const;
    void writeJudgePoints(const std::string &filepath) const;
    void writeJudgePoints(const std::string &filepath, const std::vector<const JudgePoint*> &jpts) const;

    void refine();
    void updateJptsInCell(Chd chd,
                          std::priority_queue<std::pair<zsw::Scalar,JudgePoint*>,std::vector<std::pair<zsw::Scalar,JudgePoint*>>,
                          ErrorMaxComparison> *err_queue);
    void checkUpNormalCondition(Chd chd, std::queue<Chd> &chds_queue,
                                std::unordered_set<std::string> *cell_key_set_pre,
                                std::unordered_set<std::string> *cell_key_set_cur);

    bool isSatisfyErrorBound(const std::vector<VertexTriple> &bound_tris,
                             const std::vector<const JudgePoint*> &jpts_in_bbox,
                             const Eigen::Matrix<zsw::Scalar,3,1> &merge_pt,
                             const zsw::Scalar v_pt,
                             std::vector<VertexUpdateData> &vup,
                             std::vector<JudgePointUpdateData> * jup=nullptr) const;

    bool isTolTetsSatisfyNormalCondition(const std::vector<VertexTriple> &bound_tris,
                                         const Eigen::Matrix<zsw::Scalar,3,1> &pt,
                                         const PointType point_type) const;

    void constructKernelRegionJudger(const std::vector<VertexTriple> &bound_tris,
                                     std::vector<Vhd> &opposite_vs, KernelRegionJudger &krj) const;

    bool tryCollapseBoundaryEdge(TTds::Edge &e,
                                 std::unordered_map<std::string,TTds::Edge> &edge_map);

    bool tryCollapseZeroEdge(TTds::Edge &e,
                             std::unordered_map<std::string,TTds::Edge> &z_map,
                             std::unordered_map<std::string,TTds::Edge> *bz_map);

    bool tryCollapseBZEdge(TTds::Edge &e,
                           std::unordered_map<std::string, TTds::Edge> &bz_map,
                           std::unordered_map<std::string,TTds::Edge> *z_map);

    //void createJudgePoints();

    void updateVertex(const std::vector<VertexUpdateData> &vup)
    {
      std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
    }

    void boundaryEdgeBack(Vhd vhd, std::unordered_map<std::string, TTds::Edge> &edge_map) const;
    void zeroEdgeBack(Vhd vhd, std::unordered_map<std::string,TTds::Edge> &edge_map) const;
    void bzEdgeBack(Vhd vhd, std::unordered_map<std::string, TTds::Edge> &edge_map) const;

    void calcJptsInBbox(Vhd *vhd, const size_t n, std::vector<const JudgePoint*> &jpts_in_bbox) const;

    void sampleAdjCells(const TTds::Edge &e, std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &sample_points) const;



    // data access

    // for debug
    void testTdsValid();
    void testCollapse();
  private:
    std::vector<JudgePoint> jpts_;
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> inner_jpts_;
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> outer_jpts_;
    std::shared_ptr<zsw::Flann<zsw::Scalar>> inner_kdtree_ptr_;
    std::shared_ptr<zsw::Flann<zsw::Scalar>> outer_kdtree_ptr_;
    std::shared_ptr<TriangulationWapper> tw_;
    zsw::Scalar err_epsilon_;
    zsw::Scalar tri_sample_r_;
    zsw::Scalar tet_sample_r_;
    zsw::Scalar normal_cond_scale_;
  };

  void calcVerticesBbox(Vhd *vhd_ptr, const size_t n, Eigen::Matrix<zsw::Scalar,3,2> &bbox);
  void writePoints(const std::string &filepath, const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &pts);
  void calcCircumcenter(const Eigen::Matrix<zsw::Scalar,3,4> &tri_pts, Eigen::Matrix<zsw::Scalar,3,1> &center);
}
#endif /* ISOTOPIC_APPROXIMATION_H */

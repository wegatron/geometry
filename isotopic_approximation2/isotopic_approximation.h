#ifndef ISOTOPIC_APPROXIMATION_H
#define ISOTOPIC_APPROXIMATION_H

#include <unordered_map>
#include <queue>
#include <zswlib/config.h>
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
    Approximation() { normal_cond_scale_=0.3; tmp_outdir_="/home/wegatron/tmp/"; need_smooth_=false; }
    void setTmpOutDir(const std::string &tmp_outdir) { tmp_outdir_=tmp_outdir; }
    void setNeedSmooth(bool need_smooth) { need_smooth_=need_smooth; }
    void init(const zsw::Scalar err_epsilon,
              const zsw::Scalar tri_sample_r,
              const zsw::Scalar tet_sample_r,
              std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &inner_jpts,
              std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &outer_jpts,
              std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bs_jpts);
    void simp(const std::string &tmp_output_dir);

    void writeZeroSurface(const std::string &filepath) const;
    void writeTetMesh(const std::string &filepath,
                      std::vector<std::function<bool(const TTds::Cell_handle)>> ignore_tet_funcs,
                      const TTds *tds_ptr=nullptr) const;

    void writeAdjcentCells(const std::string &filepath, const TTds::Edge &e) const;
    void writeAdjcentCells(const std::string &filepath, const std::vector<Chd> &chds) const;
    void writeJudgePoints(const std::string &filepath) const;
    void writeJudgePoints(const std::string &filepath, const std::vector<const JudgePoint*> &jpts) const;

    void refine();
    void simpTolerance();
    void mutuallTessellation(TTds *tds_ptr=nullptr);
  private:
    void simpZeroSurface(std::unordered_map<std::string,TTds::Edge> *z_map=nullptr,
                         std::unordered_map<std::string,TTds::Edge> *bz_map=nullptr);
    bool simpBZEdges(std::unordered_map<std::string,TTds::Edge> *bz_map=nullptr,
                     std::unordered_map<std::string,TTds::Edge> *z_map=nullptr);

    zsw::Scalar updateJptsInCell(Chd chd,
                                 /*std::priority_queue<std::pair<zsw::Scalar,JudgePoint*>,std::vector<std::pair<zsw::Scalar,JudgePoint*>>,
                                   ErrorMaxComparison> *err_queue*/
                                 std::vector<JudgePoint*> * updated_jpts);

    void updateJptsInCells(const std::vector<Chd> &chds,     std::priority_queue<std::pair<zsw::Scalar,JudgePoint*>,
                           std::vector<std::pair<zsw::Scalar,JudgePoint*>>,
                           ErrorMaxComparison> &err_queue);

    bool checkUpNormalCondition(Chd chd, std::vector<Chd> &chds_queue);
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
    void updateVertex(const std::vector<VertexUpdateData> &vup)
    {
      std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
    }
    void boundaryEdgeBack(Vhd vhd, std::unordered_map<std::string, TTds::Edge> &edge_map) const;
    void zeroEdgeBack(Vhd vhd, std::unordered_map<std::string,TTds::Edge> &edge_map) const;
    void bzEdgeBack(Vhd vhd, std::unordered_map<std::string, TTds::Edge> &edge_map) const;
    void calcJptsInBbox(Vhd *vhd, const size_t n, std::vector<const JudgePoint*> &jpts_in_bbox) const;
    void sampleAdjCells(const TTds::Edge &e, std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &sample_points) const;

    void smoothBoundary();
    void updateAllBoundaryVerticesMaxDis();
    void laplaceSmoothZeroSurface();
    void bilateralSmoothZeroSurface();
    void updateAllZeroVerticesMaxDis();
    void calcBoundaryOneRing(Vhd vhd, std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &ring_pts) const;
    void calcZeroSurfaceOneRing(Vhd vhd, std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &ring_pts) const;
  public:
    bool checkNormalCondition() const;
    void testTdsValid();
    void testCollapse();
    void testKdtree() const;
  private:
    std::vector<JudgePoint> jpts_;
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> inner_jpts_;
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> outer_jpts_;
    KdTreeWarper inner_kdtree_;
    KdTreeWarper outer_kdtree_;
    std::shared_ptr<TriangulationWapper> tw_;
    zsw::Scalar err_epsilon_;
    zsw::Scalar tri_sample_r_;
    zsw::Scalar tet_sample_r_;
    zsw::Scalar normal_cond_scale_;
    std::string tmp_outdir_;
    bool need_smooth_;
  };

  zsw::Scalar calcZeroTetHeight(Chd chd);
  zsw::Scalar calcBoundaryTetHeight(Chd chd);
  zsw::Scalar calcTetHeightType0(const Eigen::Matrix<zsw::Scalar,3,1> &pt0, const Eigen::Matrix<zsw::Scalar,3,1> &pt1,
                                 const Eigen::Matrix<zsw::Scalar,3,1> &pt2, const Eigen::Matrix<zsw::Scalar,3,1> &pt3);
  zsw::Scalar calcTetHeightType1(const Eigen::Matrix<zsw::Scalar,3,1> &pt0, const Eigen::Matrix<zsw::Scalar,3,1> &pt1,
                                 const Eigen::Matrix<zsw::Scalar,3,1> &pt2, const Eigen::Matrix<zsw::Scalar,3,1> &pt3);
  void calcVerticesBbox(Vhd *vhd_ptr, const size_t n, Eigen::Matrix<zsw::Scalar,3,2> &bbox);
  void writePoints(const std::string &filepath, const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &pts);
  void calcCircumcenter(const Eigen::Matrix<zsw::Scalar,3,4> &tri_pts, Eigen::Matrix<zsw::Scalar,3,1> &center);
}
#endif /* ISOTOPIC_APPROXIMATION_H */

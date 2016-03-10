#ifndef ISOTOPIC_APPROXIMATION_H
#define ISOTOPIC_APPROXIMATION_H

#include <unordered_map>
#include <queue>
#include <zswlib/config.h>
#include <zswlib/zsw_clock_c11.h>
#include "basic_data_structure.h"
#include "constraint.h"

namespace zsw {

  struct JudgePointUpdateData
  {
    JudgePoint *jpt;
    zsw::Scalar val_cur_;
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
    Approximation() {
      alpha_=0.1;
      normal_cond_scale_=0.3;
      g_scale_ = 1.0;
      sample_tet_by_ref_ = false;

      print_sp_size_ = false;
      tmp_outdir_="/home/wegatron/tmp/";
      bz_judge_pt_cnt_=0;
      bz_krj_need_judge_cnt_=0;
      bz_normal_cond_judge_cnt_=0;
      bz_error_bound_judge_cnt_=0;
    }
    void setTetByRef(bool sample_tet_by_ref) { sample_tet_by_ref_ = sample_tet_by_ref; }
    void setTmpOutDir(const std::string &tmp_outdir) { tmp_outdir_=tmp_outdir; }
    void init(const zsw::Scalar err_epsilon,
              const zsw::Scalar tri_sample_r,
              const zsw::Scalar tet_sample_r,
              std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &inner_jpts,
              std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &outer_jpts,
              std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bs_jpts);

    void initD(const zsw::Scalar err_epsilon,
               const zsw::Scalar tri_sample_r,
               const zsw::Scalar tet_sample_r,
               std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &inner_jpts,
               std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &outer_jpts,
               std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bs_jpts,
               std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_inner_jpts,
               std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_outer_jpts,
               std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_bs_jpts
               );

    void simp(const std::string &tmp_output_dir);

    void writeZeroSurface(const std::string &filepath) const;
    void writeTetMesh(const std::string &filepath,
                      std::vector<std::function<bool(const TTds::Cell_handle)>> ignore_tet_funcs,
                      const TTds *tds_ptr = nullptr,
                      bool cur_pts = true) const;

    void writeAdjcentCells(const std::string &filepath, const TTds::Edge &e) const;
    void writeAdjcentCells(const std::string &filepath, const std::vector<Chd> &chds) const;
    void writeJudgePoints(const std::string &filepath) const;
    void writeJudgePoints(const std::string &filepath, const std::vector<const JudgePoint*> &jpts) const;
    void setGscale(zsw::Scalar g_scale) { g_scale_ = g_scale; }
    void refine(const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bs_jpts);
    void refineD(std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &ori_bs_jpts,
                 std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_inner_jpts,
                 std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_outer_jpts,
                 std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_bs_jpts);
    void simpTolerance();
    void mutuallTessellation(TTds *tds_ptr=nullptr);
  private:
    void simpZeroSurface();
    void simpBZEdges();

    void updateJptsInCell(Chd chd, std::vector<JudgePoint*> * updated_jpts, bool using_cur_pts);
    void updateJptsInCells(const std::vector<Chd> &chds,     std::priority_queue<std::pair<zsw::Scalar,JudgePoint*>,
                           std::vector<std::pair<zsw::Scalar,JudgePoint*>>,
                           ErrorMaxComparison> &err_queue);

    void updateJptsInCellsD(const std::vector<Chd> &chds,     std::priority_queue<std::pair<zsw::Scalar,JudgePoint*>,
                            std::vector<std::pair<zsw::Scalar,JudgePoint*>>,
                            ErrorMaxComparison> &err_queue);

    bool checkUpNormalCondition(Chd chd, std::vector<Chd> &chds_queue, bool using_cur_pts);
    bool isSatisfyErrorBound(const std::vector<VertexTriple> &bound_tris,
                             const std::vector<const JudgePoint*> &jpts_in_bbox,
                             const Eigen::Matrix<zsw::Scalar,3,1> &merge_pt,
                             const zsw::Scalar v_pt,
                             std::vector<JudgePointUpdateData> * jup=nullptr);
    bool isTetsSatisfyNormalCondition(const std::vector<VertexTriple> &bound_tris,
                                      const Eigen::Matrix<zsw::Scalar,3,1> &pt,
                                      const PointType point_type) const;
    void constructKernelRegionJudger(const std::vector<VertexTriple> &bound_tris,
                                     std::vector<Vhd> &opposite_vs, KernelRegionJudger &krj) const;
    bool tryCollapseBoundaryEdge(TTds::Edge &e,
                                 std::unordered_map<std::string,TTds::Edge> &edge_map);

    bool tryCollapseBZEdge(TTds::Edge &e,
                           std::queue<std::pair<TTds::Edge, size_t>> &z_q,
                           size_t cur_update,
                           bool is_bz_back = false);

    void boundaryEdgeBack(Vhd vhd, std::unordered_map<std::string, TTds::Edge> &edge_map) const;
    void zeroEdgeBack(Vhd vhd, std::queue<std::pair<TTds::Edge, size_t>> &z_q, size_t cur_update) const;
    void bzEdgeBack(Vhd vhd, std::queue<std::pair<TTds::Edge, size_t>> &zb_q, size_t cur_update) const;
    void calcJptsInBbox(Vhd *vhd, const size_t n, std::vector<const JudgePoint*> &jpts_in_bbox, bool using_cur_pts) const;
    void sampleAdjCells(const TTds::Edge &e, std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &sample_points) const;
    //void sampleKrj(const zsw::KernelRegionJudger &krj, std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &sample_points) const;
    size_t countZeroPoints() const;
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
    zsw::Scalar alpha_;
    zsw::Scalar g_scale_;
    std::string tmp_outdir_;
    zsw::common::ClockC11 clock_;
    bool sample_tet_by_ref_;

    bool print_sp_size_;
    size_t bz_judge_pt_cnt_;
    size_t bz_krj_need_judge_cnt_;
    size_t bz_normal_cond_judge_cnt_;
    size_t bz_error_bound_judge_cnt_;
  };

  zsw::Scalar calcZeroTetHeight(Chd chd);
  zsw::Scalar calcBoundaryTetHeight(Chd chd);
  zsw::Scalar calcTetHeightType0(const Eigen::Matrix<zsw::Scalar,3,1> &pt0, const Eigen::Matrix<zsw::Scalar,3,1> &pt1,
                                 const Eigen::Matrix<zsw::Scalar,3,1> &pt2, const Eigen::Matrix<zsw::Scalar,3,1> &pt3);
  zsw::Scalar calcTetHeightType1(const Eigen::Matrix<zsw::Scalar,3,1> &pt0, const Eigen::Matrix<zsw::Scalar,3,1> &pt1,
                                 const Eigen::Matrix<zsw::Scalar,3,1> &pt2, const Eigen::Matrix<zsw::Scalar,3,1> &pt3);
  void calcVerticesBbox(Vhd *vhd_ptr, const size_t n, Eigen::Matrix<zsw::Scalar,3,2> &bbox);
  void calcVerticesBboxD(Vhd *vhd_ptr, const size_t n, Eigen::Matrix<zsw::Scalar,3,2> &bbox);
  void writePoints(const std::string &filepath, const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &pts);
  void calcCircumcenter(const Eigen::Matrix<zsw::Scalar,3,4> &tri_pts, Eigen::Matrix<zsw::Scalar,3,1> &center);
}
#endif /* ISOTOPIC_APPROXIMATION_H */

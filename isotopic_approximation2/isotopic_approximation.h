#ifndef ISOTOPIC_APPROXIMATION_H
#define ISOTOPIC_APPROXIMATION_H

#include <unordered_map>
#include <zswlib/config.h>
#include <zswlib/mesh/zsw_flann.h>
#include "basic_data_structure.h"
#include "constraint.h"

namespace zsw {

  class Approximation final
  {
  public:
    Approximation() {}
    void init(const zsw::Scalar &surf_sample_r, const zsw::Scalar &tet_sample_r,
              const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bi_vertices,
              const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bo_vertices);
    void simpTolerance();
    void mutuallTessellation();
    void simpZeroSurface();
    void writeZeroSurface(const std::string &filepath) const;
    void writeTetMesh(const std::string &filepath,
                      std::vector<std::function<bool(const TTds::Cell_handle)>> ignore_tet_funcs) const;
    void writeJudgePoints(const std::string &filepath) const;
  private:
    void constructKernelRegionJudger(const std::vector<Fhd> &bound_tris, KernelRegionJudger &krj) const;
    void tryCollapseBoundaryEdge(TTds::Edge &e,
                                 std::unordered_map<std::string,TTds::Edge> &edge_map);
    void createJudgePoints();
    std::vector<JudgePoint> jpts_;
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bi_jpts_;
    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> bo_jpts_;
    std::shared_ptr<zsw::Flann<zsw::Scalar>> jpts_ptr_bi_;
    std::shared_ptr<zsw::Flann<zsw::Scalar>> jpts_ptr_bo_;
    std::shared_ptr<TriangulationWapper> tw_;
    zsw::Scalar surf_sample_r_;
    zsw::Scalar tet_sample_r_;
  };

}
#endif /* ISOTOPIC_APPROXIMATION_H */

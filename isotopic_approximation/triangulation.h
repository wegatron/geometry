#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <zswlib/mesh/mesh_type.h>

#include "cgal_common.h"

#define ZSW_DEBUG
#define FAKE_KERNEL_REGION_POINT

namespace zsw {

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

  class TetMesh
  {
  public:
    struct Vertex
    {
      char pt_type_; // outer_boundary 1, zero_set 0, inner_boundary -1
      Point pt_;
      size_t father_;
      std::vector<size_t> tet_ids_;
      std::vector<size_t> edge_ids_;
    };

    struct Tet
    {
      size_t vind0_;
      size_t vind1_;
      size_t vind2_;
      size_t vind3_;
      bool valid_;
    };

    struct Edge
    {
      size_t vind0_;
      size_t vind1_;
      size_t fv_cnt_;
      bool valid_;
      std::vector<size_t> fv_; // the vertex id, that construct a valid face with this edge. the vertex construct a face with this edge must be sorted
    };

    TetMesh(const std::vector<Point> bz_points, const std::vector<Point> &bo_points, const std::vector<Point> &bi_points,
            const zsw::Scalar sample_dense);

    void simplify();
    void cleanVertices();
    void writeVtk(const std::string &filepath) const;
    void writeZeroSetSurface(const std::string &filepath);
#ifdef ZSW_DEBUG
    bool testCollapseEdge(size_t vind0, size_t vind1);
#endif
    const std::vector<Edge>& getEdges() const { return edges_; }
    const std::vector<Vertex>& getVertices() const { return vertices_; }
    const std::vector<Tet>& getTets() const { return tets_; }
  private:
    void addTets(const Delaunay &td, size_t &tet_id);
    void addEdges(const Delaunay &ti, const Delaunay &to);
    void updateFv(Edge &e);

    bool collapseEdge(Edge &e);
    void collapseEdge(Edge &e, const Point &pt);
    bool findKernelRegionPoint(const Edge &e, Point &pt) const;

    zsw::Scalar sample_dense_;
    std::vector<Vertex> vertices_;
    std::vector<Tet> tets_;
    std::vector<Edge> edges_;
    std::vector<Point> sample_points_;
  };
}


#endif /* TRIANGULATION_H */

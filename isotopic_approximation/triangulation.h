#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <zswlib/mesh/mesh_type.h>

#include "cgal_common.h"

namespace zsw {

  /* class KernelRegionJudger */
  /* { */
  /* public: */
  /*   KernelRegionJudger(); */
  /*   void addConstraint(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> &v1, */
  /*                      const Eigen::Matrix<zsw::Scalar,3,1> &v2, const Eigen::Matrix<zsw::Scalar,3,1> &vr); */
  /*   bool judge(const Point &pt); */
  /* private: */
  /*   std::vector<Eigen::Matrix<zsw::Scalar,3,1>> vec_v0; */
  /*   std::vector<Eigen::Matrix<zsw::Scalar,3,1>> vec_v1; */
  /*   std::vector<Eigen::Matrix<zsw::Scalar,3,1>> vec_v2; */
  /*   std::vector<Eigen::Matrix<zsw::Scalar,3,1>> vec_vr; */
  /* }; */

  class TetMesh
  {
  public:
    struct Vertex
    {
      char pt_type_; // outer_boundary 1, zero_set 0, inner_boundary -1
      Point pt_;
      size_t father_;
      std::vector<size_t> tet_ids_;
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
      std::vector<size_t> fv_; // the vertex construct a face with this edge
    };

    TetMesh(const std::vector<Point> bz_points, const std::vector<Point> &bo_points, const std::vector<Point> &bi_points,
            const zsw::Scalar sample_dense);

    void simplify();
    void writeVtk(const std::string &filepath);
    void writeZeroSetSurface(const std::string &filepath);
  private:
    void cleanVertices();
    void addTets(const Delaunay &td, size_t &tet_id);
    void addEdges(const Delaunay &ti, const Delaunay &to);
    bool collapseZEdge(Edge &edge);
    bool isValidEdge(Edge &edge)
    {
      for(; vertices_[edge.vind0_].father_!=-1; edge.vind0_=vertices_[edge.vind0_].father_);
      for(; vertices_[edge.vind1_].father_!=-1; edge.vind1_=vertices_[edge.vind1_].father_);
      return edge.vind0_ != edge.vind1_;
    }

    zsw::Scalar sample_dense_;
    std::vector<Vertex> vertices_;
    std::vector<Tet> tets_;
    std::vector<Edge> edges_;
    std::vector<Point> sample_points_;
  };
}


#endif /* TRIANGULATION_H */

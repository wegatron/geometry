#ifndef BILATER_NORMAL_FILTER_H
#define BILATER_NORMAL_FILTER_H

#include <vector>
#include <set>
#include <jtflib/mesh/trimesh.h>

#include <Eigen/Sparse>

//#define ONE_RING_I 1

#define DEBUG 1

namespace zsw
{
  class BilateralNormalFilter final
  {
  public:
    enum RING_TYPE {
      ONE_EDGE_RING,
      ONE_VERTEX_RING
    };
    BilateralNormalFilter();
    void setSt(size_t st) { st_ = st; }
    void setRingType(RING_TYPE ring_type) { ring_type_ = ring_type; }
    void filter(jtf::mesh::tri_mesh &trimesh);
  private:

    void preProcess(const jtf::mesh::tri_mesh &trimesh);
    double calBc(const size_t fid, const zjucad::matrix::matrix<size_t> &mesh, const zjucad::matrix::matrix<double> &node);

    void filterNormal(jtf::mesh::tri_mesh &trimesh);
    void updateVertex(jtf::mesh::tri_mesh &trimesh);

    zjucad::matrix::matrix<double> fc_; // face_center_
    std::vector<std::set<size_t>> one_ring_;
    Eigen::SparseMatrix<double> weight_;

    RING_TYPE ring_type_;
    size_t st_;  // number of iterations of smooth normal
    size_t ut_; // number of iterations of update vertexes
    double b_c_; // 2*sigma_c^2
    double b_s_; // 2*sigma_s^2
  };

  void processEdgeOneRing(const jtf::mesh::tri_mesh &trimesh, std::vector<std::set<size_t>> &one_ring);

  void processVertexOneRing(const jtf::mesh::tri_mesh &trimesh, std::vector<std::set<size_t>> &one_ring);

  void writeTriMesh(const std::string &filename, const zjucad::matrix::matrix<size_t> &mesh,
                    const zjucad::matrix::matrix<double> &node,
                    const zjucad::matrix::matrix<double> &normal);

  void writeVtk(const std::string &filename, const zjucad::matrix::matrix<size_t> &mesh,
                const zjucad::matrix::matrix<double> &node);
}

#endif /* BILATER_NORMAL_FILTER_H */

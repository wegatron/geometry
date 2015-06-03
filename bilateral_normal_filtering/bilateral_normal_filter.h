#ifndef BILATER_NORMAL_FILTER_H
#define BILATER_NORMAL_FILTER_H

#include <vector>
#include <jtflib/mesh/trimesh.h>

#define ONE_RING_I 1

#define DEBUG 1

namespace zsw
{
  class BilateralNormalFilter final
  {
  public:
    BilateralNormalFilter();
    void filter(jtf::mesh::tri_mesh &trimesh);
    void setSt(size_t st) { st_ = st; }
  private:
    void filterNormal(jtf::mesh::tri_mesh &trimesh);
    void updateVertex(jtf::mesh::tri_mesh &trimesh);
    // return false if face is boundary face
    bool queryFidOneRingI(const size_t fid, const jtf::mesh::tri_mesh &trimesh, std::vector<size_t> &fid_one_ring);
    void queryFidOneRingII(const size_t fid, const zjucad::matrix::matrix<size_t> &mesh, std::vector<size_t> &fid_one_ring);
    void preProcess(const jtf::mesh::tri_mesh &trimesh);
    double calBc(const size_t fid, const std::vector<size_t> &fid_one_ring, const zjucad::matrix::matrix<size_t> &mesh,
                 const zjucad::matrix::matrix<double> &node);

    std::multimap<size_t, size_t> v2f_;
    size_t st_;  // number of iterations of smooth normal
    size_t ut_; // number of iterations of update vertexes
    double b_c_; // 2*sigma_c^2
    double b_s_; // 2*sigma_s^2
  };

  void writeTriMesh(const std::string &filename, const zjucad::matrix::matrix<size_t> &mesh,
                    const zjucad::matrix::matrix<double> &node,
                    const zjucad::matrix::matrix<double> &normal);

  void writeVtk(const std::string &filename, const zjucad::matrix::matrix<size_t> &mesh,
                const zjucad::matrix::matrix<double> &node);
}

#endif /* BILATER_NORMAL_FILTER_H */

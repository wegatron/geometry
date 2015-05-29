#ifndef BILATER_NORMAL_FILTER_H
#define BILATER_NORMAL_FILTER_H

#include <vector>
#include <jtflib/mesh/trimesh.h>

namespace zsw
{
  class BilateralNormalFilter final
  {
  public:
    BilateralNormalFilter();
    void filter(jtf::mesh::tri_mesh &trimesh);
  private:
    void filterNormal(jtf::mesh::tri_mesh &trimesh);
    void updateVertex(jtf::mesh::tri_mesh &trimesh);
    // return false if face is boundary face
    bool queryFidOneRing(const size_t fid, const jtf::mesh::tri_mesh &trimesh, std::vector<size_t> &fid_one_ring);

    size_t st_;  // number of iterations of smooth normal
    size_t ut_; // number of iterations of update vertexes
    double b_c_; // 2*sigma_c^2
    double b_s_; // 2*sigma_s^2
  };

  void writeTriMesh(const std::string &filename, const zjucad::matrix::matrix<size_t> &mesh,
                    const zjucad::matrix::matrix<double> &node,
                    const zjucad::matrix::matrix<double> &normal);
}

#endif /* BILATER_NORMAL_FILTER_H */

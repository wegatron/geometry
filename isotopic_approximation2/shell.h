#ifndef SHELL_H
#define SHELL_H

#include <vector>
#include <zswlib/mesh/mesh_type.h>

namespace zsw{
  void genAndSampleShell(zsw::mesh::TriMesh &input_mesh,
                         const zsw::Scalar err_epsilon,
                         const zsw::Scalar tri_sample_r,
                         std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &inner_jpts,
                         std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &outer_jpts,
                         std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bs_jpts);
}

#endif /* SHELL_H */

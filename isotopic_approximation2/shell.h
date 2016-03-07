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

  void genAndSampleShellD(zsw::mesh::TriMesh &input_mesh,
                            zsw::mesh::TriMesh &deformed_mesh,
                            const zsw::Scalar err_epsilon,
                            const zsw::Scalar tri_sample_r,
                            std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &inner_jpts,
                            std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &outer_jpts,
                            std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bs_jpts,
                            std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_inner_jpts,
                            std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_outer_jpts,
                            std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_bs_jpts);

  void calcBoundSphere(const zsw::mesh::TriMesh &mesh,
                       Eigen::Matrix<zsw::Scalar,3,1> &center,
                       zsw::Scalar &radius);

  zsw::Scalar calcDeformScale(zsw::mesh::TriMesh &ori_mesh, zsw::mesh::TriMesh &deformed_mesh);

  void calcTriMeshBBox(const zsw::mesh::TriMesh &mesh, Eigen::Matrix<zsw::Scalar,3,2> &bbox);
}

#endif /* SHELL_H */

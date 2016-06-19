#ifndef RES_ANALYSIS_H
#define RES_ANALYSIS_H

#include <jtflib/mesh2/trimesh.h>
#include "config.h"

namespace zsw
{
  namespace mesh{
    void ZSW_API diffMeshVertex(jtf::mesh2::meshes &mesh_a, jtf::mesh2::meshes &mesh_b,
                                const std::string &output_vtk);
    void ZSW_API diffMeshNormal(const jtf::mesh2::meshes &mesh_a, const jtf::mesh2::meshes &mesh_b,
                                const std::string &output_vtk);
    void ZSW_API sortMeshVertex(Eigen::Matrix<size_t, 3, Eigen::Dynamic> &mesh, Eigen::Matrix<double, 3, Eigen::Dynamic> &node);
  }
}

#endif /* RES_ANALYSIS_H */

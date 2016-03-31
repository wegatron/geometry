#ifndef SURFACE_GENERATOR_H
#define SURFACE_GENERATOR_H

#include <vector>
#include <zswlib/mesh/mesh_type.h>

namespace zsw
{
    void genPoints(const zsw::Scalar dis, zsw::mesh::TriMesh &tm,
                   std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bi_points,
                   std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bo_points);
}


#endif /* SURFACE_GENERATOR_H */

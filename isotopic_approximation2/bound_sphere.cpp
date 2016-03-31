#include "bound_sphere.h"
#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>

using namespace std;

zsw::BoundSphere::BoundSphere(const std::string &filepath,
                              const zsw::Scalar scale,
                              const Eigen::Matrix<zsw::Scalar,3,1> &transform)
{
  zsw::mesh::TriMesh mesh;
  if(!OpenMesh::IO::read_mesh(mesh, filepath)) {
    std::cerr << "[ERROR] can't read mesh:" << filepath << std::endl;
  }
  vertices_.resize(mesh.n_vertices());
  auto bs_vit = vertices_.begin();
  auto m_vit = mesh.vertices_begin();
  while(bs_vit!=vertices_.end()) {
    *bs_vit = (scale * mesh.point(*m_vit))+transform;
    ++bs_vit; ++m_vit;
  }
}

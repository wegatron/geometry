#ifndef MESH_OP_H
#define MESH_OP_H

#include <zswlib/config.h>
#include <zswlib/mesh/mesh_type.h>
#include <zswlib/data_type.h>

namespace zsw
{
  namespace mesh {
    void ZSW_API rRingVertex(const zsw::mesh::TriMesh &tm, const size_t r, std::vector<zsw::FakeSet<size_t>> &ring);
  }
}


#endif /* MESH_OP_H */

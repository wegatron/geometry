#ifndef MESH_TYPE_H
#define MESH_TYPE_H

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <zswlib/config.h>

namespace zsw
{
  namespace mesh {
#ifdef USING_DOUBLE_PRECISION
    typedef OpenMesh::TriMesh_ArrayKernelT<OpenMesh::EigenDoubleTraits> TriMesh;
#else
    typedef OpenMesh::TriMesh_ArrayKernelT<OpenMesh::EigenFloatTraits> TriMesh;
#endif
  }
}

#endif /* MESH_TYPE_H */

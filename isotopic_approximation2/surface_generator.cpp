#include <OpenMesh/Core/IO/MeshIO.hh>
#include "surface_generator.h"

#include <iostream>
//#define ZSW_DEBUG

void zsw::genPoints(const zsw::Scalar dis, zsw::mesh::TriMesh &tm,
                    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bi_points,
                    std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bo_points)
{
  if(!tm.has_vertex_normals()) {
    tm.request_face_normals();
    tm.request_vertex_normals();
    tm.update_normals();
  }
#ifdef ZSW_DEBUG
  zsw::mesh::TriMesh tmo(tm);
  zsw::mesh::TriMesh tmi(tm);
#endif
  for(zsw::mesh::TriMesh::ConstVertexIter vit=tm.vertices_begin(); vit!=tm.vertices_end(); ++vit) {
    Eigen::Matrix<zsw::Scalar, 3, 1> ep = tm.point(*vit);
    Eigen::Matrix<zsw::Scalar, 3, 1> offset = tm.normal(*vit) * dis;
#ifdef ZSW_DEBUG
    tmi.set_point(zsw::mesh::TriMesh::VertexHandle(vit->idx()), ep-offset);
    tmo.set_point(zsw::mesh::TriMesh::VertexHandle(vit->idx()), ep+offset);
#endif
    bo_points.push_back(ep+offset);
    bi_points.push_back(ep-offset);
  }

#ifdef ZSW_DEBUG
  OpenMesh::IO::write_mesh(tmi, "/home/wegatron/tmp/gen_point_tmi.obj");
  OpenMesh::IO::write_mesh(tmo, "/home/wegatron/tmp/gen_point_tmo.obj");
#endif
}

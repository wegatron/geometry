#include <OpenMesh/Core/IO/MeshIO.hh>
#include "shell.h"

#include <iostream>
#include <fstream>

#include "sampling.h"

namespace zsw{
  void boundSphere(const std::string &filepath,
                   const zsw::Scalar scale,
                   const Eigen::Matrix<zsw::Scalar,3,1> &transform,
                   std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bs_jpts)
  {
    zsw::mesh::TriMesh mesh;
    if(!OpenMesh::IO::read_mesh(mesh, filepath)) {
      std::cerr << "[ERROR] can't read mesh:" << filepath << std::endl;
    }
    bs_jpts.resize(mesh.n_vertices());
    auto bs_vit = bs_jpts.begin();
    auto m_vit = mesh.vertices_begin();
    while(bs_vit!=bs_jpts.end()) {
      *bs_vit = (scale * mesh.point(*m_vit))+transform;
      ++bs_vit; ++m_vit;
    }
  }

  void genAndSampleShell(zsw::mesh::TriMesh &input_mesh,
                         const zsw::Scalar err_epsilon,
                         const zsw::Scalar tri_sample_r,
                         std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &inner_jpts,
                         std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &outer_jpts,
                         std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bs_jpts)
  {
    if(!input_mesh.has_vertex_normals()) {
      input_mesh.request_face_normals();
      input_mesh.request_vertex_normals();
      input_mesh.update_normals();
    }
    Eigen::Matrix<zsw::Scalar,3,2> bbox;
    bbox.block<3,1>(0,0) = input_mesh.point(*input_mesh.vertices_begin());
    bbox.block<3,1>(0,1) = bbox.block<3,1>(0,0);
    for(auto fit=input_mesh.faces_begin(); fit!=input_mesh.faces_end(); ++fit) {
      Eigen::Matrix<zsw::Scalar,3,3> in_tri, out_tri;
      size_t i=0;
      for(zsw::mesh::TriMesh::FaceVertexIter fvit=input_mesh.fv_iter(*fit); fvit.is_valid(); ++fvit) {
        Eigen::Matrix<zsw::Scalar,3,1> offset=input_mesh.normal(*fvit)*err_epsilon;
        in_tri.block<3,1>(0,i)=input_mesh.point(*fvit)-offset;
        out_tri.block<3,1>(0,i)=input_mesh.point(*fvit)+offset;
        for(size_t di=0; di<3; ++di) {
          if(out_tri(di,i)<bbox(di,0)) { bbox(di,0)=out_tri(di,i); }
          else if(out_tri(di,i)>bbox(di,1)){ bbox(di,1)=out_tri(di,i); }
        }
        ++i;
      }

      assert(i==3);

      sampleTriangle(in_tri, tri_sample_r, inner_jpts);
      sampleTriangle(out_tri, tri_sample_r, outer_jpts);
    }
    zsw::Scalar scale=0.5*(bbox.block<3,1>(0,1)-bbox.block<3,1>(0,0)).norm();
    Eigen::Matrix<zsw::Scalar,3,1> transform =0.5*(bbox.block<3,1>(0,0)+bbox.block<3,1>(0,1));
    boundSphere("bound_sphere.obj", scale, transform, bs_jpts);
  }
}

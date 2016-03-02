#include <OpenMesh/Core/IO/MeshIO.hh>
#include "shell.h"

#include <iostream>
#include <fstream>

#include "sampling.h"

//#define TEST_OUR_IDEA

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
#ifdef TEST_OUR_IDEA
      std::cout << "Testing our idea x scaled to 0.25!!!" << std::endl;
      in_tri.block<1,3>(0,0)*=0.25;
      out_tri.block<1,3>(0,0)*=0.25;
#endif
      sampleTriangle(in_tri, tri_sample_r, inner_jpts);
      sampleTriangle(out_tri, tri_sample_r, outer_jpts);
    }
    zsw::Scalar scale=0.5*(bbox.block<3,1>(0,1)-bbox.block<3,1>(0,0)).norm();
    Eigen::Matrix<zsw::Scalar,3,1> transform =0.5*(bbox.block<3,1>(0,0)+bbox.block<3,1>(0,1));
    boundSphere("bound_sphere.obj", scale, transform, bs_jpts);
  }

  zsw::Scalar genAndSampleDeformedShell(zsw::mesh::TriMesh &ori_mesh,
                                        zsw::mesh::TriMesh &deformed_mesh,
                                        const zsw::Scalar err_epsilon,
                                        const zsw::Scalar tri_sample_r,
                                        std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &inner_jpts,
                                        std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &outer_jpts,
                                        std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bs_jpts)
  {
    zsw::Scalar deform_scale=calcDeformScale(ori_mesh, deformed_mesh);
    std::cout << "!!!!!!!!!!!!!!!!scale=" << deform_scale << std::endl;
    if(!deformed_mesh.has_vertex_normals()) {
      deformed_mesh.request_face_normals();
      deformed_mesh.request_vertex_normals();
      deformed_mesh.update_normals();
    }
    // if(!ori_mesh.has_vertex_normals()) {
    //   ori_mesh.request_face_normals();
    //   ori_mesh.request_vertex_normals();
    //   ori_mesh.update_normals();
    // }
    // Eigen::Matrix<zsw::Scalar,3,2> ori_bbox;
    // ori_bbox.block<3,1>(0,0) = ori_mesh.point(*ori_mesh.vertices_begin());
    // ori_bbox.block<3,1>(0,1) = ori_bbox.block<3,1>(0,0);

    Eigen::Matrix<zsw::Scalar,3,2> bbox;
    bbox.block<3,1>(0,0) = deformed_mesh.point(*deformed_mesh.vertices_begin());
    bbox.block<3,1>(0,1) = bbox.block<3,1>(0,0);
    //Eigen::Matrix<zsw::Scalar,3,3> ori_in_tri, ori_out_tri;
    Eigen::Matrix<zsw::Scalar,3,3> deformed_in_tri, deformed_out_tri;
    const size_t nf=ori_mesh.n_faces();
    for(size_t fi=0; fi<nf; ++fi) {
      zsw::mesh::TriMesh::FaceHandle fh=zsw::mesh::TriMesh::FaceHandle(int(fi));
      size_t vi=0;
      // for(zsw::mesh::TriMesh::CFVIter fv_it=ori_mesh.cfv_iter(fh); fv_it.is_valid(); ++fv_it) {
      //   Eigen::Matrix<zsw::Scalar,3,1> ori_offset=ori_mesh.normal(*fv_it)*err_epsilon;
      //   ori_in_tri.block<3,1>(0,vi)=ori_mesh.point(*fv_it)-ori_offset;
      //   ori_out_tri.block<3,1>(0,vi)=ori_mesh.point(*fv_it)+ori_offset;
      //   ++vi;
      //   for(size_t di=0; di<3; ++di) {
      //     if(ori_out_tri(di,vi)<ori_bbox(di,0)) { ori_bbox(di,0)=ori_out_tri(di,vi); }
      //     else if(ori_out_tri(di,vi)>ori_bbox(di,1)){ ori_bbox(di,1)=ori_out_tri(di,vi); }
      //   }
      // }
      vi=0;
      for(zsw::mesh::TriMesh::CFVIter fv_it=deformed_mesh.cfv_iter(fh); fv_it.is_valid(); ++fv_it) {
        Eigen::Matrix<zsw::Scalar,3,1> deformed_offset=deformed_mesh.normal(*fv_it)*err_epsilon*deform_scale;
        deformed_in_tri.block<3,1>(0,vi)=deformed_mesh.point(*fv_it)-deformed_offset;
        deformed_out_tri.block<3,1>(0,vi)=deformed_mesh.point(*fv_it)+deformed_offset;
        for(size_t di=0; di<3; ++di) {
          if(deformed_out_tri(di,vi)<bbox(di,0)) { bbox(di,0)=deformed_out_tri(di,vi); }
          else if(deformed_out_tri(di,vi)>bbox(di,1)){ bbox(di,1)=deformed_out_tri(di,vi); }
        }
        ++vi;
      }
      // sampleTriangleRefTriangle(deformed_in_tri, ori_in_tri, tri_sample_r, inner_jpts);
      // sampleTriangleRefTriangle(deformed_out_tri, ori_out_tri, tri_sample_r, outer_jpts);
      sampleTriangle(deformed_in_tri, tri_sample_r*deform_scale, inner_jpts);
      sampleTriangle(deformed_out_tri, tri_sample_r*deform_scale, outer_jpts) ;
    }
    zsw::Scalar scale=0.5*(bbox.block<3,1>(0,1)-bbox.block<3,1>(0,0)).norm();
    Eigen::Matrix<zsw::Scalar,3,1> transform =0.5*(bbox.block<3,1>(0,0)+bbox.block<3,1>(0,1));
    boundSphere("bound_sphere.obj", scale, transform, bs_jpts);
    return deform_scale;
  }

  void genAndSampleAllShell(zsw::mesh::TriMesh &ori_mesh,
                            zsw::mesh::TriMesh &deformed_mesh,
                            const zsw::Scalar err_epsilon,
                            const zsw::Scalar tri_sample_r,
                            std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &inner_jpts,
                            std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &outer_jpts,
                            std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bs_jpts,
                            std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_inner_jpts,
                            std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_outer_jpts,
                            std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &deformed_bs_jpts)
  {
    if(!ori_mesh.has_vertex_normals()) {
      ori_mesh.request_face_normals();
      ori_mesh.request_vertex_normals();
      ori_mesh.update_normals();
    }
    if(!deformed_mesh.has_vertex_normals()) {
      deformed_mesh.request_face_normals();
      deformed_mesh.request_vertex_normals();
      deformed_mesh.update_normals();
    }
    Eigen::Matrix<zsw::Scalar,3,1> ori_center;
    zsw::Scalar ori_radius;
    Eigen::Matrix<zsw::Scalar,3,1> deformed_center;
    zsw::Scalar deformed_radius;
    calcBoundSphere(ori_mesh, ori_center, ori_radius);
    calcBoundSphere(deformed_mesh, deformed_center, deformed_radius);
    zsw::Scalar deformed_err_epsilon = err_epsilon * deformed_radius / ori_radius;
    boundSphere("bound_sphere.obj", ori_radius+err_epsilon, ori_center, bs_jpts);
    boundSphere("bound_sphere.obj", deformed_radius+deformed_err_epsilon, deformed_center, deformed_bs_jpts);

    Eigen::Matrix<zsw::Scalar,3,3> ori_in_tri, ori_out_tri;
    Eigen::Matrix<zsw::Scalar,3,3> deformed_in_tri, deformed_out_tri;
    const size_t nf=ori_mesh.n_faces();
    for(size_t fi=0; fi<nf; ++fi) {
      zsw::mesh::TriMesh::FaceHandle fh=zsw::mesh::TriMesh::FaceHandle(int(fi));
      size_t vi=0;
      for(zsw::mesh::TriMesh::CFVIter fv_it=ori_mesh.cfv_iter(fh); fv_it.is_valid(); ++fv_it) {
        Eigen::Matrix<zsw::Scalar,3,1> ori_offset=ori_mesh.normal(*fv_it)*err_epsilon;
        ori_in_tri.block<3,1>(0,vi)=ori_mesh.point(*fv_it)-ori_offset;
        ori_out_tri.block<3,1>(0,vi)=ori_mesh.point(*fv_it)+ori_offset;
        ++vi;
      }
      vi=0;
      for(zsw::mesh::TriMesh::CFVIter fv_it=deformed_mesh.cfv_iter(fh); fv_it.is_valid(); ++fv_it) {
        Eigen::Matrix<zsw::Scalar,3,1> deformed_offset=deformed_mesh.normal(*fv_it)*deformed_err_epsilon;
        deformed_in_tri.block<3,1>(0,vi)=deformed_mesh.point(*fv_it)-deformed_offset;
        deformed_out_tri.block<3,1>(0,vi)=deformed_mesh.point(*fv_it)+deformed_offset;
        ++vi;
      }
      sampleTriangleAndDeformedTriangle(deformed_in_tri, ori_in_tri, tri_sample_r,deformed_inner_jpts, inner_jpts);
      sampleTriangleAndDeformedTriangle(deformed_out_tri, ori_out_tri, tri_sample_r, deformed_outer_jpts, outer_jpts);
    }
  }

  void calcBoundSphere(const zsw::mesh::TriMesh &mesh,
                       Eigen::Matrix<zsw::Scalar,3,1> &center,
                       zsw::Scalar &radius)
  {
    Eigen::Matrix<zsw::Scalar,3,2> bbox;
    calcTriMeshBBox(mesh, bbox);
    center = ( bbox.block<3,1>(0,0) + bbox.block<3,1>(0,1) ) * 0.5;
    radius = ( bbox.block<3,1>(0,0) - bbox.block<3,1>(0,1) ).norm() * 0.5;
  }

  zsw::Scalar calcDeformScale(zsw::mesh::TriMesh &ori_mesh, zsw::mesh::TriMesh &deformed_mesh)
  {
    Eigen::Matrix<zsw::Scalar,3,2> ori_bbox, deformed_bbox;
    calcTriMeshBBox(ori_mesh, ori_bbox);
    calcTriMeshBBox(deformed_mesh, deformed_bbox);
    zsw::Scalar n=(ori_bbox.block<3,1>(0,0)-ori_bbox.block<3,1>(0,1)).norm();
    zsw::Scalar m=(deformed_bbox.block<3,1>(0,0)-deformed_bbox.block<3,1>(0,1)).norm();
    return m/n;
  }

  void calcTriMeshBBox(const zsw::mesh::TriMesh &mesh, Eigen::Matrix<zsw::Scalar,3,2> &bbox)
  {
    bbox.block<3,1>(0,0) = mesh.point( *mesh.vertices_begin() );
    bbox.block<3,1>(0,1) = bbox.block<3,1>(0,0);
    for(auto vit=mesh.vertices_begin(); vit!=mesh.vertices_end(); ++vit) {
      Eigen::Matrix<zsw::Scalar,3,1> tmp_pt=mesh.point(*vit);
      for(size_t di=0; di<3; ++di) {
        if(tmp_pt[di] < bbox(di,0)) { bbox(di,0) = tmp_pt[di]; }
        else if(tmp_pt[di] > bbox(di,1)) { bbox(di,1) = tmp_pt[di]; }
      }
    }
  }
}

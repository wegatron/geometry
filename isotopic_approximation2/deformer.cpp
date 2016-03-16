#include "deformer.h"

#include <iostream>
#include <Eigen/Dense>

#include "sampling.h"

namespace zsw
{
  void GaussDisWeightFunc::calcWeight(
                                      const Eigen::Matrix<zsw::Scalar,3,1> &vt,
                                      const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &vs,
                                      std::vector<zsw::Scalar> &weight)
  {
    const zsw::Scalar e = 2.718281828459;
    zsw::Scalar total = 0;
#pragma omp parallel for
    for(size_t i=0; i<vs.size(); ++i) {
      weight[i] = pow(e, c_*(vt-vs[i]).squaredNorm());
#pragma omp atomic
      total += weight[i];
    }
#pragma omp parallel for
    for(size_t i=0; i<vs.size(); ++i) {      weight[i] /= total;    }
  }

  void LocalTranslateDeformFunc::calcDeformedPos(
                                                 const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &vs,
                                                 const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &dvs,
                                                 const Eigen::Matrix<zsw::Scalar,3,1> &vt,
                                                 Eigen::Matrix<zsw::Scalar,3,1> &dvt)
  {
    assert(dvs.size() == vs.size());
    std::vector<zsw::Scalar> weight(vs.size(), 0);
    dis_weight_func_->calcWeight(vt, vs, weight);
    dvt = vt;
    for(size_t i=0; i<vs.size(); ++i) { dvt += weight[i] * (dvs[i] - vs[i]); }
  }

  Deformer::Deformer(const zsw::mesh::TriMesh &ori_mesh,
                     const zsw::mesh::TriMesh &deformed_mesh,
                     const zsw::Scalar sample_r)
  {
    // sample two mesh
    Eigen::Matrix<zsw::Scalar,3,3> ori_tri, deformed_tri;
    const size_t nf=ori_mesh.n_faces();
    for(size_t fi=0; fi<nf; ++fi) {
      zsw::mesh::TriMesh::FaceHandle fh=zsw::mesh::TriMesh::FaceHandle(int(fi));
      size_t vi=0;
      for(zsw::mesh::TriMesh::CFVIter fv_it=ori_mesh.cfv_iter(fh); fv_it.is_valid(); ++fv_it) {
        ori_tri.block<3,1>(0,vi) = ori_mesh.point(*fv_it);
      }
      vi=0;
      for(zsw::mesh::TriMesh::CFVIter fv_it=deformed_mesh.cfv_iter(fh); fv_it.is_valid(); ++fv_it) {
        deformed_tri.block<3,1>(0,vi++) = deformed_mesh.point(*fv_it);
      }
      zsw::sampleTriangleAndDeformedTriangle(deformed_tri, ori_tri, sample_r, ref_dvs_, ref_vs_);
    }

    // build two kdtree
    vs_kdt_.buildTree(ref_vs_[0].data(), ref_vs_.size());
    dvs_kdt_.buildTree(ref_dvs_[0].data(), ref_dvs_.size());
  }

  LocalTranslateDeformer::LocalTranslateDeformer(const zsw::mesh::TriMesh &ori_mesh,
                                                 const zsw::mesh::TriMesh &deformed_mesh,
                                                 const zsw::Scalar sample_r)
    : Deformer(ori_mesh, deformed_mesh, sample_r)  { ref_r_ = sample_r; }

  void LocalTranslateDeformer::deformTo(const std::vector<zsw::Vector3s> &vs, std::vector<zsw::Vector3s> &dvs)
  {
    std::vector<zsw::KdTreeNode> ref_nodes;
    std::vector<zsw::Vector3s> ref_vs_cur, ref_dvs_cur;
    for(size_t i=0; i<vs.size(); ++i) {
      ref_nodes.clear();
      vs_kdt_.findWithinR(vs[i], ref_r_, std::back_inserter(ref_nodes));
      ref_vs_cur.clear(); ref_dvs_cur.clear();
      for(const zsw::KdTreeNode &node : ref_nodes) {
        ref_vs_cur.push_back(ref_vs_[node.index_]); ref_dvs_cur.push_back(ref_dvs_[node.index_]);
      }
      df_.calcDeformedPos(ref_vs_cur, ref_dvs_cur, vs[i], dvs[i]);
    }
  }

  void LocalTranslateDeformer::deformBack(const std::vector<zsw::Vector3s> &dvs, std::vector<zsw::Vector3s> &vs)
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  }

}

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

  void LocalTranslateDeformer::calcDeformedPos(
                                               const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &vs,
                                               const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &dvs,
                                               const Eigen::Matrix<zsw::Scalar,3,1> &vt,
                                               Eigen::Matrix<zsw::Scalar,3,1> &dvt) const
  {
    assert(dvs.size() == vs.size());
    std::vector<zsw::Scalar> weight(vs.size(), 0);
    dis_weight_func_->calcWeight(vt, vs, weight);
    dvt = vt;
    for(size_t i=0; i<vs.size(); ++i) { dvt += weight[i] * (dvs[i] - vs[i]); }
  }

  Deformer::Deformer(const zsw::mesh::TriMesh &ori_mesh,
                     const zsw::mesh::TriMesh &deformed_mesh,
                     const zsw::Scalar sample_r,
                     std::shared_ptr<zsw::DisWeightFunc> dis_weight_func) : dis_weight_func_(dis_weight_func)
  {
    // sample two mesh
    Eigen::Matrix<zsw::Scalar,3,3> ori_tri, deformed_tri;
    const size_t nf=ori_mesh.n_faces();
    for(size_t fi=0; fi<nf; ++fi) {
      zsw::mesh::TriMesh::FaceHandle fh=zsw::mesh::TriMesh::FaceHandle(int(fi));
      size_t vi=0;
      for(zsw::mesh::TriMesh::CFVIter fv_it=ori_mesh.cfv_iter(fh); fv_it.is_valid(); ++fv_it) {
        ori_tri.block<3,1>(0,vi++) = ori_mesh.point(*fv_it);
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
                                                 const zsw::Scalar sample_r,
                                                 std::shared_ptr<DisWeightFunc> dis_weight_func)
    : Deformer(ori_mesh, deformed_mesh, sample_r, dis_weight_func)  { ref_r_ = 5 * sample_r; }

  void LocalTranslateDeformer::deformTo(const std::vector<zsw::Vector3s> &vs, std::vector<zsw::Vector3s> &dvs) const
  {
    std::vector<zsw::KdTreeNode> ref_nodes;
    std::vector<zsw::Vector3s> ref_vs_cur, ref_dvs_cur;
    dvs.resize(vs.size());
    for(size_t i=0; i<vs.size(); ++i) {
      ref_nodes.clear();
      vs_kdt_.findWithinR(vs[i], ref_r_*ref_r_, std::back_inserter(ref_nodes));
      ref_vs_cur.clear(); ref_dvs_cur.clear();
      for(const zsw::KdTreeNode &node : ref_nodes) {
        ref_vs_cur.push_back(ref_vs_[node.index_]); ref_dvs_cur.push_back(ref_dvs_[node.index_]);
      }
      calcDeformedPos(ref_vs_cur, ref_dvs_cur, vs[i], dvs[i]);
    }
  }

  void LocalTranslateDeformer::deformBack(const std::vector<zsw::Vector3s> &dvs, std::vector<zsw::Vector3s> &vs) const
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  }

  LocalVectorFieldDeformer::LocalVectorFieldDeformer(const zsw::mesh::TriMesh &ori_mesh,
                                                     const zsw::mesh::TriMesh &deformed_mesh,
                                                     const zsw::Scalar sample_r, std::shared_ptr<zsw::DisWeightFunc> dis_weight_func)
    : Deformer(ori_mesh, deformed_mesh, sample_r, dis_weight_func)
  {
    ref_r_ = 2 * sample_r;
    calcJac();
  }

  void LocalVectorFieldDeformer::calcJac()
  {
    jac_.resize(ref_vs_.size());
    // for each vertex, search ref_r_
    for(size_t i=0; i<ref_vs_.size(); ++i) {
      std::vector<KdTreeNode> nodes;
      vs_kdt_.findWithinR(ref_vs_[i], ref_r_, std::back_inserter(nodes));
      std::vector<zsw::Vector3s> ref_vs_cur(nodes.size());
      for(size_t j=0; j<nodes.size(); ++j) { ref_vs_cur[j] = ref_vs_[nodes[j].index_]; }

      std::vector<zsw::Scalar> weight(nodes.size(), 0);
      dis_weight_func_->calcWeight(ref_vs_[i], ref_vs_cur, weight);
      // for each vertex in ref_r_,  minmize  \sum w|| e' - j*e||, calc j_{k,0} j_{k,1} j_{k,2}
      Eigen::Matrix<zsw::Scalar,3,3> A = Eigen::Matrix<zsw::Scalar,3,3>::Zero();
      Eigen::Matrix<zsw::Scalar,3,1> B[3];
      B[0].setZero(); B[1].setZero(); B[2].setZero();
      for(size_t j=0; j<nodes.size(); ++j) {
        Eigen::Matrix<zsw::Scalar,3,1> e = ref_vs_[i] - ref_vs_cur[j];
        Eigen::Matrix<zsw::Scalar,3,1> de = ref_dvs_[i] - ref_dvs_[nodes[j].index_];
        A += weight[j] * e * e.transpose();
        B[0] += weight[j] * de[0] * e;
        B[1] += weight[j] * de[1] * e;
        B[2] += weight[j] * de[2] * e;
      }
#if 0
      zsw::Vector3s exp_b[3];
      exp_b[0] << 2, 0, 0;      exp_b[1] << 0, 2, 0;      exp_b[2] << 0, 0, 2;
      if((A*exp_b[0] - B[0]).norm() > 1e-3) {
        std::cout << "err0!!! : " << (A*exp_b[0] - B[0]).norm() << std::endl;
        abort();
      }
      if((A*exp_b[1] - B[1]).norm() > 1e-3) {
        std::cout << "err1!!! : " << (A*exp_b[1] - B[1]).norm() << std::endl;
        abort();
      }
      if((A*exp_b[2] - B[2]).norm() > 1e-3) {
        std::cout << "err2!!! : " << (A*exp_b[2] - B[2]).norm() << std::endl;
        abort();
      }
#endif
      Eigen::PartialPivLU<Eigen::Matrix<zsw::Scalar,3,3>> pplu;
      pplu.compute(A);
      jac_[i].block<1,3>(0,0) = pplu.solve(B[0]).transpose();
      jac_[i].block<1,3>(1,0) = pplu.solve(B[1]).transpose();
      jac_[i].block<1,3>(2,0) = pplu.solve(B[2]).transpose();

#if 0
      if((A*jac_[i].block<1,3>(0,0).transpose() - B[0]).norm() > 1e-3) {
        std::cout << "ferr0!!! : " << (A*jac_[i].block<1,3>(0,0).transpose() - B[0]).norm() << std::endl;
        abort();
      }
      if((A*jac_[i].block<1,3>(1,0).transpose() - B[1]).norm() > 1e-3) {
        std::cout << "ferr1!!! : " << (A*jac_[i].block<1,3>(1,0).transpose() - B[1]).norm() << std::endl;
        abort();
      }
      if((A*jac_[i].block<1,3>(2,0).transpose() - B[2]).norm() > 1e-3) {
        std::cout << "ferr2!!! : " << (A*jac_[i].block<1,3>(2,0).transpose() - B[2]).norm() << std::endl;
        abort();
      }
#endif
    }
  }

  void LocalVectorFieldDeformer::deformBack(const std::vector<zsw::Vector3s> &dvs, std::vector<zsw::Vector3s> &vs) const
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  }

  void LocalVectorFieldDeformer::deformTo(const std::vector<zsw::Vector3s> &vs, std::vector<zsw::Vector3s> &dvs) const
  {
    // for each vertex v_i of vs
    // find all v_j in ref_r, and calc v'_i = w_{ij}(J_j(v_i - v_j) + v'_j')
    for(size_t i=0; i<vs.size(); ++i) {
      std::vector<zsw::KdTreeNode> nodes;
      vs_kdt_.findWithinR(vs[i], ref_r_, std::back_inserter(nodes));
      std::vector<zsw::Vector3s> ref_vs_cur(nodes.size());
      for(size_t j=0; j<nodes.size(); ++j) { ref_vs_cur[i] = ref_vs_[nodes[i].index_]; }
      std::cout << __FILE__ << __LINE__ << std::endl;
      std::cout << "size:" << ref_vs_cur.size() << std::endl;
      std::vector<zsw::Scalar> weight(ref_vs_cur.size(), 0);
      dis_weight_func_->calcWeight(vs[i], ref_vs_cur, weight);
      dvs[i].setZero();
      std::cout << __FILE__ << __LINE__ << std::endl;
      for(size_t j=0; j<nodes.size(); ++j) {
        const size_t ref_index = nodes[j].index_;
        dvs[i] += weight[j] * (jac_[ref_index] * (vs[i] - ref_vs_[ref_index]) + ref_dvs_[ref_index]);
      }
    }
  }

}

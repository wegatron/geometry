#include "deformer.h"

#include <iostream>
#include <queue>
#include <Eigen/Dense>

#include "sampling.h"
#include "isotopic_approximation.h"

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

  // void LocalTranslateDeformer::calcDeformedPos(
  //                                              const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &vs,
  //                                              const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &dvs,
  //                                              const Eigen::Matrix<zsw::Scalar,3,1> &vt,
  //                                              Eigen::Matrix<zsw::Scalar,3,1> &dvt) const
  // {
  //   assert(dvs.size() == vs.size());
  //   std::vector<zsw::Scalar> weight(vs.size(), 0);
  //   dis_weight_func_->calcWeight(vt, vs, weight);
  //   dvt = vt;
  //   for(size_t i=0; i<vs.size(); ++i) { dvt += weight[i] * (dvs[i] - vs[i]); }
  // }

  Deformer::Deformer(const zsw::mesh::TriMesh &ori_mesh,
                     const zsw::mesh::TriMesh &deformed_mesh,
                     const zsw::Scalar sample_r,
                     const size_t near_count,
                     std::shared_ptr<zsw::DisWeightFunc> dis_weight_func) : dis_weight_func_(dis_weight_func)
  {
    near_count_ = near_count;
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
      zsw::Vector3s ori_normal =
        (ori_tri.block<3,1>(0,1) - ori_tri.block<3,1>(0,0)).cross(ori_tri.block<3,1>(0,2) - ori_tri.block<3,1>(0,0));
      zsw::Vector3s deformed_normal =
        (deformed_tri.block<3,1>(0,1) - deformed_tri.block<3,1>(0,0)).cross(deformed_tri.block<3,1>(0,2) - deformed_tri.block<3,1>(0,0));
      ori_normal.normalize(); deformed_normal.normalize();
      size_t cnt = ref_dvs_.size();
      zsw::sampleTriangleAndDeformedTriangle(deformed_tri, ori_tri, sample_r, ref_dvs_, ref_vs_);
      cnt = ref_dvs_.size() - cnt;
      ref_vs_normal_.insert(ref_vs_normal_.end(), cnt, ori_normal);
      ref_dvs_normal_.insert(ref_dvs_normal_.end(), cnt, deformed_normal);
    }

    // build two ann
    vs_ann_.reset(new Flann<zsw::Scalar>(ref_vs_[0].data(), ref_vs_.size()));
    dvs_ann_.reset(new Flann<zsw::Scalar>(ref_dvs_[0].data(), ref_dvs_.size()));
  }

  // LocalTranslateDeformer::LocalTranslateDeformer(const zsw::mesh::TriMesh &ori_mesh,
  //                                                const zsw::mesh::TriMesh &deformed_mesh,
  //                                                const zsw::Scalar sample_r,
  //                                                std::shared_ptr<DisWeightFunc> dis_weight_func)
  //   : Deformer(ori_mesh, deformed_mesh, sample_r, dis_weight_func)  { ref_r_ = 5 * sample_r; }

  // void LocalTranslateDeformer::deformTo(const std::vector<zsw::Vector3s> &vs, std::vector<zsw::Vector3s> &dvs) const
  // {
  //   std::vector<zsw::KdTreeNode> ref_nodes;
  //   std::vector<zsw::Vector3s> ref_vs_cur, ref_dvs_cur;
  //   dvs.resize(vs.size());
  //   for(size_t i=0; i<vs.size(); ++i) {
  //     ref_nodes.clear();
  //     vs_kdt_.findWithinR(vs[i], ref_r_*ref_r_, std::back_inserter(ref_nodes));
  //     ref_vs_cur.clear(); ref_dvs_cur.clear();
  //     for(const zsw::KdTreeNode &node : ref_nodes) {
  //       ref_vs_cur.push_back(ref_vs_[node.index_]); ref_dvs_cur.push_back(ref_dvs_[node.index_]);
  //     }
  //     calcDeformedPos(ref_vs_cur, ref_dvs_cur, vs[i], dvs[i]);
  //   }
  // }

  // void LocalTranslateDeformer::deformBack(const std::vector<zsw::Vector3s> &dvs, std::vector<zsw::Vector3s> &vs) const
  // {
  //   std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  // }

  LocalVectorFieldDeformer::LocalVectorFieldDeformer(const zsw::mesh::TriMesh &ori_mesh,
                                                     const zsw::mesh::TriMesh &deformed_mesh,
                                                     const zsw::Scalar sample_r,
                                                     const size_t near_count,
                                                     std::shared_ptr<zsw::DisWeightFunc> dis_weight_func)
    : Deformer(ori_mesh, deformed_mesh, sample_r, near_count, dis_weight_func)
  {
    calcJac();
  }

  void LocalVectorFieldDeformer::calcJac()
  {
    jac_.resize(ref_vs_.size());
    // for each vertex, search near_count neighbour
    std::vector<std::vector<size_t>> indices;
    std::vector<std::vector<zsw::Scalar>> dists;
    vs_ann_->queryKnn(ref_vs_, indices, dists, near_count_);
    valid_.resize(ref_vs_.size(), false);
    for(size_t i=0; i<ref_vs_.size(); ++i) {
      std::vector<zsw::Vector3s> ref_vs_cur(near_count_);
      for(size_t j=0; j<near_count_; ++j) { ref_vs_cur[j] = ref_vs_[indices[i][j]]; }
      std::vector<zsw::Scalar> weight(near_count_, 0);
      dis_weight_func_->calcWeight(ref_vs_[i], ref_vs_cur, weight);
      // for each vertex in ref_r_,  minmize  \sum w|| e' - j*e||, calc j_{k,0} j_{k,1} j_{k,2}
      Eigen::Matrix<zsw::Scalar,3,3> A = Eigen::Matrix<zsw::Scalar,3,3>::Zero();
      Eigen::Matrix<zsw::Scalar,3,1> B[3];
      B[0].setZero(); B[1].setZero(); B[2].setZero();
      for(size_t j=0; j<near_count_; ++j) {
        Eigen::Matrix<zsw::Scalar,3,1> e = ref_vs_[i] - ref_vs_cur[j];
        Eigen::Matrix<zsw::Scalar,3,1> de = ref_dvs_[i] - ref_dvs_[indices[i][j]];
        A += weight[j] * e * e.transpose();
        B[0] += weight[j] * de[0] * e;
        B[1] += weight[j] * de[1] * e;
        B[2] += weight[j] * de[2] * e;
      }
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix<zsw::Scalar,3,3>> eigensolver(A);
      if (eigensolver.info() != Eigen::Success) abort();
      zsw::Vector3s ev = eigensolver.eigenvalues(); ev.normalize();
      valid_[i] = fabs(ev[0])>6e-4 && fabs(ev[1])>6e-4 && fabs(ev[2])>6e-4;

      if(valid_[i]) {
        Eigen::LDLT<Eigen::Matrix<zsw::Scalar,3,3>> ldlt;
        ldlt.compute(A);
        jac_[i].block<1,3>(0,0) = ldlt.solve(B[0]).transpose();
        jac_[i].block<1,3>(1,0) = ldlt.solve(B[1]).transpose();
        jac_[i].block<1,3>(2,0) = ldlt.solve(B[2]).transpose();
#if 0
        zsw::Vector3s exp_b[3];
        exp_b[0] << 2, 0, 0;      exp_b[1] << 0, 2, 0;      exp_b[2] << 0, 0, 2;
        Eigen::Matrix<zsw::Scalar,3,3> jac_exp;
        jac_exp.block<1,3>(0,0) = exp_b[0].transpose();
        jac_exp.block<1,3>(1,0) = exp_b[1].transpose();
        jac_exp.block<1,3>(2,0) = exp_b[2].transpose();
        std::cout << "err:" << (jac_exp-jac_[i]).norm() << std::endl;
        // std::cout << "ldt info():" << (ldlt.info()==Eigen::Success) << std::endl;
        // std::cout << "det(A)/A.norm():" << A.determinant()/A.norm() << std::endl;
        // Eigen::SelfAdjointEigenSolver<Eigen::Matrix<zsw::Scalar,3,3>> eigensolver(A);
        // if (eigensolver.info() != Eigen::Success) abort();
        // zsw::Vector3s ev = eigensolver.eigenvalues(); ev.normalize();
        // std::cout << "The eigenvalues of A:" << ev.transpose() << std::endl;
#endif
      }
#if 0
      {
        if(!valid_[i]) {
          zsw::Vector3s exp_b[3];
          exp_b[0] << 1, 0, 0;      exp_b[1] << 0, 1, 0;      exp_b[2] << 0, 0, 1;
          jac_[i].block<1,3>(0,0) = exp_b[0].transpose();
          jac_[i].block<1,3>(1,0) = exp_b[1].transpose();
          jac_[i].block<1,3>(2,0) = exp_b[2].transpose();
        }
        continue;
        if((A*jac_[i].block<1,3>(0,0).transpose() - B[0]).norm() > 1e-5) {
          std::cout << "ferr0!!! : " << (A*jac_[i].block<1,3>(0,0).transpose() - B[0]).norm() << std::endl;
          abort();
        }
        if((A*jac_[i].block<1,3>(1,0).transpose() - B[1]).norm() > 1e-5) {
          std::cout << "ferr1!!! : " << (A*jac_[i].block<1,3>(1,0).transpose() - B[1]).norm() << std::endl;
          abort();
        }
        if((A*jac_[i].block<1,3>(2,0).transpose() - B[2]).norm() > 1e-5) {
          std::cout << "ferr2!!! : " << (A*jac_[i].block<1,3>(2,0).transpose() - B[2]).norm() << std::endl;
          abort();
        }
      }
      {
        zsw::Vector3s exp_b[3];
        exp_b[0] << 2, 0, 0;      exp_b[1] << 0, 2, 0;      exp_b[2] << 0, 0, 2;
        if((A*exp_b[0] - B[0]).norm() > 1e-5) {
          std::cout << "err0!!! : " << (A*exp_b[0] - B[0]).norm() << std::endl;
          abort();
        }
        if((A*exp_b[1] - B[1]).norm() > 1e-5) {
          std::cout << "err1!!! : " << (A*exp_b[1] - B[1]).norm() << std::endl;
          abort();
        }
        if((A*exp_b[2] - B[2]).norm() > 1e-5) {
          std::cout << "err2!!! : " << (A*exp_b[2] - B[2]).norm() << std::endl;
          abort();
        }
        // Eigen::Matrix<zsw::Scalar,3,3> jac_exp;
        // jac_exp.block<1,3>(0,0) = exp_b[0].transpose();
        // jac_exp.block<1,3>(1,0) = exp_b[1].transpose();
        // jac_exp.block<1,3>(2,0) = exp_b[2].transpose();
        // if((jac_exp-jac_[i]).norm() > 1e-3) {          std::cout << "jac:\n" << jac_[i] << std::endl;       abort();  }
      }
#endif
    }
    writeValidRefVs();
    // resolve invalid jacobian
    resolveInvalidJacobian(indices);
  }

  class CntMaxComp
  {
  public:
    // lv is after rv ?
    bool operator()(const std::pair<size_t, size_t> &lv,
                    const std::pair<size_t,size_t> &rv){
      return lv.first < rv.first;
    }
  };

  void LocalVectorFieldDeformer::resolveInvalidJacobian(const std::vector<std::vector<size_t>> &indices)
  {
    std::vector<size_t> valid_ne_cnt(ref_vs_.size(), 0);
    std::priority_queue<std::pair<size_t, size_t>, std::vector<std::pair<size_t, size_t>>, CntMaxComp> max_n_q;
    for(size_t i=0; i<ref_vs_.size(); ++i) {
      if(valid_[i]) { continue; }
      size_t valid_neighbour = 0;
      for(size_t j=0; j<near_count_; ++j) {        if(valid_[indices[i][j]]) { ++valid_neighbour; }      }
      max_n_q.push(std::make_pair(valid_neighbour, i));
      valid_ne_cnt[i] = valid_neighbour;
    }
    std::cout << "invalid:" << max_n_q.size() << std::endl;
    // size_t resolved_cnt = 0;
    while(!max_n_q.empty()) {
      std::pair<size_t, size_t> n_ind = max_n_q.top(); max_n_q.pop();
      const size_t cur_ind = n_ind.second;
      if(valid_[cur_ind]) { continue; }
      /// resolve n_ind->second
      std::vector<zsw::Vector3s> ref_vs_cur;
      for(auto ind : indices[cur_ind]) { ref_vs_cur.push_back(ref_vs_[ind]); }
      std::vector<zsw::Scalar> weight(near_count_, 0);
      dis_weight_func_->calcWeight(ref_vs_[cur_ind], ref_vs_cur, weight);
      // calc cur_normal scale from neighbour normal scale
      zsw::Scalar scale = 0;
      zsw::Scalar total_weight = 0;
      size_t true_valid_cnt = 0;
      for(size_t i=0; i<near_count_; ++i) {
        const size_t ind = indices[cur_ind][i];
        if(!valid_[ind]) { continue; }
        ++true_valid_cnt;
        scale += weight[i] * (jac_[ind] * ref_vs_normal_[ind]).norm();
        total_weight += weight[i];
      }
      if(true_valid_cnt != n_ind.first) { std::cout << __FILE__ << __LINE__ << std::endl; abort(); }
      scale = scale / total_weight;
      // add cur normal scale pt into A and calc jac
      zsw::Vector3s vn_ori = ref_vs_normal_[cur_ind];
      zsw::Vector3s vn_deformed = ref_dvs_normal_[cur_ind] * scale;
      //Eigen::Matrix<zsw::Scalar,3,3> A_ori = Eigen::Matrix<zsw::Scalar,3,3>::Zero();
      Eigen::Matrix<zsw::Scalar,3,3> A = 0.5 * vn_ori * vn_ori.transpose();
      zsw::Vector3s B[3];
      B[0] = 0.5 * vn_deformed[0] * vn_ori;
      B[1] = 0.5 * vn_deformed[1] * vn_ori;
      B[2] = 0.5 * vn_deformed[2] * vn_ori;
      for(size_t j=0; j<near_count_; ++j) {
        Eigen::Matrix<zsw::Scalar,3,1> e = ref_vs_[cur_ind] - ref_vs_cur[j];
        Eigen::Matrix<zsw::Scalar,3,1> de = ref_dvs_[cur_ind] - ref_dvs_[indices[cur_ind][j]];
        A += weight[j] * e * e.transpose();
        //A_ori += weight[j] * e * e.transpose();
        B[0] += weight[j] * de[0] * e;
        B[1] += weight[j] * de[1] * e;
        B[2] += weight[j] * de[2] * e;
      }

      Eigen::LDLT<Eigen::Matrix<zsw::Scalar,3,3>> ldlt;
      ldlt.compute(A);
      jac_[cur_ind].block<1,3>(0,0) = ldlt.solve(B[0]).transpose();
      jac_[cur_ind].block<1,3>(1,0) = ldlt.solve(B[1]).transpose();
      jac_[cur_ind].block<1,3>(2,0) = ldlt.solve(B[2]).transpose();

      if(!jac_[cur_ind].allFinite()) {
        std::cout << "valid_cnt:" << n_ind.first << std::endl;
        std::cout << "true_valid_cnt:" << true_valid_cnt << std::endl;
        std::cout << "total weight:" << total_weight << std::endl;
        std::cout << "scale" << scale << std::endl;
        std::cout << "failed again!!!" << std::endl;
        std::cout << "ind=" << cur_ind << std::endl;
        std::cout << "pos:" << ref_vs_[cur_ind].transpose() << std::endl;
        std::cout << "vn_ori:" << vn_ori.transpose() << std::endl;
        std::cout << "vn_dformed:" <<  vn_deformed.transpose() << std::endl;
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<zsw::Scalar,3,3>> eigensolver(A);
        if (eigensolver.info() != Eigen::Success) abort();
        zsw::Vector3s ev = eigensolver.eigenvalues(); ev.normalize();
        std::cout << "ev:" << ev.transpose() << std::endl;
        // eigensolver.compute(A_ori);
        // ev = eigensolver.eigenvalues(); ev.normalize();
        // std::cout << "ori_ev=" << ev.transpose() << std::endl;
        {
          zsw::Vector3s exp_b[3];
          exp_b[0] << 2, 0, 0;      exp_b[1] << 0, 2, 0;      exp_b[2] << 0, 0, 2;
          if((A*exp_b[0] - B[0]).norm() > 1e-5) {
            std::cout << "err0!!! : " << (A*exp_b[0] - B[0]).norm() << std::endl;
            abort();
          }
          if((A*exp_b[1] - B[1]).norm() > 1e-5) {
            std::cout << "err1!!! : " << (A*exp_b[1] - B[1]).norm() << std::endl;
            abort();
          }
          if((A*exp_b[2] - B[2]).norm() > 1e-5) {
            std::cout << "err2!!! : " << (A*exp_b[2] - B[2]).norm() << std::endl;
            abort();
          }
          // Eigen::Matrix<zsw::Scalar,3,3> jac_exp;
          // jac_exp.block<1,3>(0,0) = exp_b[0].transpose();
          // jac_exp.block<1,3>(1,0) = exp_b[1].transpose();
          // jac_exp.block<1,3>(2,0) = exp_b[2].transpose();
          // if((jac_exp-jac_[i]).norm() > 1e-3) {          std::cout << "jac:\n" << jac_[i] << std::endl;       abort();  }
        }
      }
      valid_[cur_ind] = true;
      // update the neighbour
      for(size_t ind : indices[cur_ind]) {
        if(!valid_[ind]) { max_n_q.push(std::make_pair(++valid_ne_cnt[ind], ind)); }
      }
      // if(++resolved_cnt % 50 == 0) {
      //   std::cout << "resolved cnt=" << resolved_cnt << std::endl;
      //   std::cout << "jac:\n" << jac_[cur_ind] << std::endl;
      // }
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
    std::vector<std::vector<size_t>> indices;
    std::vector<std::vector<zsw::Scalar>> dists;
    vs_ann_->queryKnn(vs, indices, dists, near_count_);
    std::vector<zsw::Vector3s> ref_vs_cur(near_count_);
    dvs.resize(vs.size(), zsw::Vector3s::Zero());
    for(size_t i=0; i<vs.size(); ++i) {
      for(size_t j=0; j<near_count_; ++j) { ref_vs_cur[j] = ref_vs_[indices[i][j]]; }
      std::vector<zsw::Scalar> weight(near_count_, 0);
      dis_weight_func_->calcWeight(vs[i], ref_vs_cur, weight);
      for(size_t j=0; j<near_count_; ++j) {
        const size_t ref_index = indices[i][j];
        dvs[i] += weight[j] * (jac_[ref_index] * (vs[i] - ref_vs_[ref_index]) + ref_dvs_[ref_index]);
        //dvs[i] += 1.0/near_count_ * (jac_[ref_index] * (vs[i] - ref_vs_[ref_index]) + ref_dvs_[ref_index]);
      }
    }
  }

  void LocalVectorFieldDeformer::writeValidRefVs() const
  {
    std::vector<zsw::Vector3s> tmp_vs;
    for(size_t i=0; i<ref_vs_.size(); ++i) {
      if(valid_[i]) {        tmp_vs.push_back(ref_vs_[i]);      }
    }
    zsw::writePoints("/home/wegatron/tmp/deform_test/valid_ref_vs.vtk", tmp_vs);
  }
}

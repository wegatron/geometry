#define _EXPORTING

#include "lowpass_filter.h"

#include <boost/foreach.hpp>
#include <zswlib/zsw_clock.h>
#include <zswlib/mesh/mesh_op.h>
#include "../basic_op.h"


zsw::mesh::LowPassFilter::LowPassFilter(std::shared_ptr<zsw::data::FilterDataInterFace> df)
{
  df_ = df;
}

void zsw::mesh::LowPassFilter::run()
{
  zsw::mesh::rRingVertex(df_->trimesh_, 1, df_->adj_ring_->ring_dat_);
//  zsw::mesh::rRingVertex(df_->trimesh_, 1, df_->adj_ring_->ring_dat_);
  zsw::mesh::trimesh2p(df_->trimesh_, p_);
  calcOpMat(df_->adj_ring_->ring_dat_, mat_l_);
  calcOpMat(df_->adj_ring_->ring_dat_, smat_l_);
  filter();
}

void zsw::mesh::LowPassFilter::calcOpMat(std::vector<zsw::FakeSet<size_t>> &ring, Eigen::Matrix<zsw::Scalar, Eigen::Dynamic, Eigen::Dynamic> &mat)
{
  size_t n = ring.size();
  mat.resize(n,n); mat.setZero(n,n);
  for(size_t i=0; i<n; ++i) {
    mat(i,i) = 1.0;
    double val = -1.0/ring[i].size();
    BOOST_FOREACH(size_t j, ring[i]) {
      mat(i,j) = val;
    }
  }
}

void zsw::mesh::LowPassFilter::calcOpMat(std::vector<zsw::FakeSet<size_t>> &ring, Eigen::SparseMatrix<zsw::Scalar> &mat)
{
  size_t n = ring.size();
  mat.resize(n,n);
  std::list<Eigen::Triplet<double>> triplet_list;
  for(size_t i=0; i<n; ++i) {
    triplet_list.push_back(Eigen::Triplet<double>(i,i, 1.0));
    double val = -1.0/ring[i].size();
    BOOST_FOREACH(size_t j, ring[i]) {
      triplet_list.push_back(Eigen::Triplet<double>(i,j, val));
    }
  }
  mat.setFromTriplets(triplet_list.begin(), triplet_list.end());
}

#if 1
void zsw::mesh::LowPassFilter::filter()
{
  std::cerr << "Eigen Decomposition ..." << std::endl;
  zsw::common::Clock clock;
  Eigen::EigenSolver<Eigen::Matrix<zsw::Scalar, Eigen::Dynamic, Eigen::Dynamic>> es(mat_l_);
  std::cout << "time cost: " << clock.time() << std::endl;
  Eigen::VectorXcd eivals = es.eigenvalues();
  std::map<zsw::Scalar, size_t> v2ind;
  for(size_t i=0; i<eivals.size(); ++i) {
    v2ind.insert(std::pair<zsw::Scalar, size_t>(eivals(i).real(),i));
  } // sort by val
//  std::cerr << __FILE__ << " " << __LINE__ << std::endl;
  Eigen::Matrix<zsw::Scalar, Eigen::Dynamic, Eigen::Dynamic> mat_e = es.eigenvectors().real();
  Eigen::Matrix<zsw::Scalar, Eigen::Dynamic, Eigen::Dynamic> mat_inv_e = mat_e.inverse();
  x_h_.resize(p_.rows(), 3);
  x_h_.col(0) = mat_inv_e*p_.col(0);
  x_h_.col(1) = mat_inv_e*p_.col(1);
  x_h_.col(2) = mat_inv_e*p_.col(2);
//  std::cerr << __FILE__ << " " << __LINE__ << std::endl;
  std::map<zsw::Scalar, size_t>::const_reverse_iterator it = v2ind.crbegin();
  for(size_t i=0; i<v2ind.size()*3/4; ++i) {
    size_t id = it->second; ++it;
//    std::cerr << __FILE__ << ":" << it->first << std::endl;
    x_h_(id,0) = x_h_(id,1) = x_h_(id,2) = 0.0;
  }
//  std::cerr << __FILE__ << " " << __LINE__ << std::endl;
  std::cerr << mat_e.rows()  << "," << mat_e.cols() << " # " << x_h_.rows() << "," << x_h_.cols() << std::endl;
  res_p_ = mat_e*x_h_;
  p2trimesh(res_p_, df_->trimesh_);
}

#else

void zsw::mesh::LowPassFilter::filter()
{
  size_t n = smat_l_.rows();
  arma::SpMat<zsw::Scalar> a_spm(n,n);
  zsw::common::Clock clock;
  {
    std::cerr << "Eigen spmat to arma spmat ..." << std::endl;
    eigen2arma(smat_l_, a_spm);
    std::cout << "time cost: " << clock.time() << std::endl;
  }
  Eigen::Matrix<zsw::Scalar, Eigen::Dynamic, 1> val(n);
  Eigen::Matrix<zsw::Scalar, Eigen::Dynamic, Eigen::Dynamic> mat_e(n, n);
  {
    std::cerr << "Eigen Decomposition ..." << std::endl;
    arma::vec eigval;
    arma::mat eigvec;
    eigs_sym(eigval, eigvec, B, n*0.5);  // find 50% eigenvalues
    arma2eigtn(eigval, val);
    arma2eigen(eigvec, mat_e);
    std::cout << "time cost: " << clock.time() << std::endl;
  }

  std::cerr << __FILE__ << " " << __LINE__ << std::endl;
  Eigen::Matrix<zsw::Scalar, Eigen::Dynamic, Eigen::Dynamic> mat_e = es.eigenvectors().real();
  Eigen::Matrix<zsw::Scalar, Eigen::Dynamic, Eigen::Dynamic> mat_inv_e = mat_e.inverse();
  x_h_.resize(p_.rows(), 3);
  x_h_.col(0) = mat_inv_e*p_.col(0);
  x_h_.col(1) = mat_inv_e*p_.col(1);
  x_h_.col(2) = mat_inv_e*p_.col(2);
  std::cerr << __FILE__ << " " << __LINE__ << std::endl;
  std::map<zsw::Scalar, size_t>::const_reverse_iterator it = v2ind.crbegin();
  for(size_t i=0; i<v2ind.size()*3/4; ++i) {
    size_t id = it->second; ++it;
    std::cerr << __FILE__ << ":" << it->first << std::endl;
    x_h_(id,0) = x_h_(id,1) = x_h_(id,2) = 0.0;
  }
  std::cerr << __FILE__ << " " << __LINE__ << std::endl;
  std::cerr << mat_e.rows()  << "," << mat_e.cols() << " # " << x_h_.rows() << "," << x_h_.cols() << std::endl;
  res_p_ = mat_e*x_h_;
  p2trimesh(res_p_, df_->trimesh_);
}

#endif

#define _EXPORTING
#include "lh_filter.h"

#include <queue>
#include <set>
#include <boost/foreach.hpp>

//#include <Eigen/CholmodSupport>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/Dense>

#include <zswlib/mesh/mesh_op.h>

#include "../basic_op.h"

zsw::mesh::LHFilter::LHFilter(std::shared_ptr<zsw::data::LHFDataInterFace> df)
{
  df_ = df;
}

int zsw::mesh::LHFilter::run()
{
  // process ring_l_ and ring_h_
  zsw::mesh::rRingVertex(df_->trimesh_, df_->l_size_, *(df_->ring_l_));
  zsw::mesh::rRingVertex(df_->trimesh_, df_->h_size_, *(df_->ring_h_));

  // construct L H
  std::cerr << "construct L H" << std::endl;
  calcOpMat(*(df_->ring_l_), mat_l_);
  calcOpMat(*(df_->ring_h_), mat_h_);

  // fill p_
  std::cerr << "fill p_" << std::endl;
  trimesh2p(df_->trimesh_, p_);

  // calc Matrix A and vect b
  std::cerr << "calc Matrix A and vect b" << std::endl;
  buildEquation();

  // solve A* res_p_ = b
  std::cerr << "solve A* res_p_ = b" << std::endl;
  if(solveEquation()) {
    return __LINE__;
  }

  // p2trimesh
  p2trimesh(res_p_, df_->trimesh_);
  return 0;
}

void zsw::mesh::LHFilter::calcOpMat(std::vector<zsw::FakeSet<size_t>> &ring, Eigen::SparseMatrix<zsw::Scalar> &mat)
{
  size_t n = ring.size();
  mat.resize(3*n, 3*n);
  std::list<Eigen::Triplet<zsw::Scalar>> triplet_list;
  for(size_t i=0; i<n; ++i) {
    triplet_list.push_back(Eigen::Triplet<zsw::Scalar>(3*i,3*i, 1));
    triplet_list.push_back(Eigen::Triplet<zsw::Scalar>(3*i+1,3*i+1, 1));
    triplet_list.push_back(Eigen::Triplet<zsw::Scalar>(3*i+2,3*i+2, 1));
    zsw::Scalar val = -1.0/ring[i].size();
    BOOST_FOREACH(size_t j, ring[i]) {
      triplet_list.push_back(Eigen::Triplet<zsw::Scalar>(3*i, 3*j, val));
      triplet_list.push_back(Eigen::Triplet<zsw::Scalar>(3*i+1, 3*j+1, val));
      triplet_list.push_back(Eigen::Triplet<zsw::Scalar>(3*i+2, 3*j+2, val));
    }
  }
  mat.setFromTriplets(triplet_list.begin(), triplet_list.end());
}

void zsw::mesh::LHFilter::buildEquation()
{
  A_ = (1-df_->alpha_)*mat_h_.transpose()*mat_h_ + df_->alpha_*mat_l_.transpose()*mat_l_;
  b_ = (1-df_->alpha_)*mat_h_.transpose()*mat_h_*p_;
}

int zsw::mesh::LHFilter::solveEquation()
{
  std::cout <<"!!!!!!!!!!" << std::endl;
  //Eigen::ConjugateGradient<Eigen::SparseMatrix<zsw::Scalar>> solver;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<zsw::Scalar>> solver;
  solver.compute(A_);
  if(solver.info()!=Eigen::Success) {
    std::cerr << "Decomposition Failed!!!" << std::endl;
    return __LINE__;
  }
  res_p_ = solver.solve(b_);
  if(solver.info()!=Eigen::Success) {
    std::cerr << "Solving Failed!!!" << std::endl;
    return __LINE__;
  }
  Eigen::Matrix<zsw::Scalar, Eigen::Dynamic, 1> diff1 = A_ * res_p_ - b_;
//  Eigen::Matrix<zsw::Scalar, Eigen::Dynamic, 1> diff2 = A_ * p_ - b_;
  std::cout << "diff1:" << diff1.norm() << std::endl;
//  std::cout << "diff2:" << diff2.norm() << std::endl;

  return 0;
}

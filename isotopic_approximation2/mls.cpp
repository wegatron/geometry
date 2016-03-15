#include "mls.h"

#include <Eigen/Dense>

namespace zsw
{
  void DeformFunc::calcDeformedPos(
                                   const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &vs,
                                   const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &dvs,
                                   const Eigen::Matrix<zsw::Scalar,3,1> &vt,
                                   Eigen::Matrix<zsw::Scalar,3,1> &dvt)
  {
    std::vector<zsw::Scalar> weight(vs.size(), 0);
    dis_weight_func_->calcWeight(vt, vs, weight);

    assert(vs.size() == dvs.size());
    A_.setZero(); B_.setZero();
    for(size_t i=0; i<vs.size(); ++i) { updateAB(vt, vs[i], dvs[i], weight[i]); }
    Eigen::Matrix<zsw::Scalar,12,1> x;
    Eigen::PartialPivLU<Eigen::Matrix<zsw::Scalar,12,12>> pplu; pplu.compute(A_);
    x = pplu.solve(B_);
    Eigen::Matrix<zsw::Scalar,3,3> a;
    Eigen::Matrix<zsw::Scalar,3,1> b;
    std::copy(x.data(), x.data()+9, a.data());
    std::copy(x.data()+9, x.data()+12, b.data());
    dvt = a * vt + b;
  }

  void DeformFunc::updateAB(const Eigen::Matrix<zsw::Scalar,3,1> &vt,
                            const Eigen::Matrix<zsw::Scalar,3,1> &v,
                            const Eigen::Matrix<zsw::Scalar,3,1> &dv,
                            const zsw::Scalar weight)
  {
    Eigen::Matrix<zsw::Scalar,12,12> tmp_A;
    Eigen::Matrix<zsw::Scalar,12,1> tmp_B;
    for(size_t i=0; i<9; i+=3) { tmp_B.block<3,1>(i,0) = v[i] * dv; }
    tmp_B.block<3,1>(9,0) = dv;
    Eigen::Matrix<zsw::Scalar,3,12> tmp_c = Eigen::Matrix<zsw::Scalar,3,12>::Zero();
    for(size_t i=0; i<3; ++i) {
      for(size_t j=0; j<3; ++j) {
        tmp_c(i, i+j*3) = v[j];
      }
      tmp_c(i, i+9) = 1;
    }
    for(size_t i=0; i<3; ++i) {      tmp_A.block<3,12>(i*3, 0) = v[i] * tmp_c;    }
    tmp_A.block<3,12>(9,0) = tmp_c;
    A_ = weight*tmp_A + A_;
    B_ = weight*tmp_B + B_;
  }

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
}

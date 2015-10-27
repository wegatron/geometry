#include "optimizer.h"

using namespace std;
using namespace Ipopt;

namespace zsw
{
  Optimizer::Optimizer(Eigen::Matrix<Number,3,1> cx) : cx_(cx) {
    cn_=0;
  }

  // ax+by+cz>=d
  void Optimizer::addConstraint(Number a, Number b, Number c, Number d)
  {
    vec_a_.push_back(a);
    vec_b_.push_back(b);
    vec_c_.push_back(c);
    vec_d_.push_back(d);
    ++cn_;
  }

  bool Optimizer::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                               Index& nnz_h_lag, IndexStyleEnum& index_style)
  {
    n=3;
    m=cn_;
    nnz_jac_g = 3*cn_;
    nnz_h_lag=0;
    index_style=TNLP::C_STYLE;
    return true;
  }

  bool Optimizer::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                  Index m, Number* g_l, Number* g_u)
  {
    assert(n==3);    assert(m==cn_); assert(cn_==vec_d_.size());
    copy(vec_d_.begin(), vec_d_.end(), g_l);
    fill(x_l, x_l+3,-10);     fill(x_u, x_u+3,10);
    fill(g_u, g_u+cn_, 1e19);
    return true;
  }

  bool Optimizer::get_starting_point(Index n, bool init_x, Number* x,
                                     bool init_z, Number* z_L, Number* z_U,
                                     Index m, bool init_lambda,
                                     Number* lambda)
  {
    assert(n==3);
    assert(init_x==true);
    x[0]=cx_[0]; x[1]=cx_[1]; x[2]=cx_[2];
    assert(init_z==false); assert(init_lambda==false);
    return true;
  }

  bool Optimizer::eval_f(Index n, const Number* x, bool new_x,
                         Number& obj_value)
  {
    assert(n==3);
    obj_value = (Eigen::Map<const Eigen::Matrix<Number,3,1>>(x)- cx_).squaredNorm();
    return true;
  }


  bool Optimizer::eval_grad_f(Index n, const Number* x, bool new_x,
                              Number* grad_f)
  {
    Eigen::Map<Eigen::Matrix<Number,3,1>> grad_f_e(grad_f);
    grad_f_e = (Eigen::Map<const Eigen::Matrix<Number,3,1>>(x)-cx_)*2;
    return true;
  }

  bool Optimizer::eval_g(Index n, const Number* x, bool new_x,
                         Index m, Number* g)
  {
    assert(m==cn_);
    Eigen::Map<Eigen::Matrix<Number, Eigen::Dynamic,1>> col_a(vec_a_.data(), vec_a_.size());
    Eigen::Map<Eigen::Matrix<Number, Eigen::Dynamic,1>> col_b(vec_b_.data(), vec_b_.size());
    Eigen::Map<Eigen::Matrix<Number, Eigen::Dynamic,1>> col_c(vec_c_.data(), vec_c_.size());
    Eigen::Map<Eigen::Matrix<Number, Eigen::Dynamic,1>>ge(g, cn_);
    ge=col_a*x[0] + col_b*x[1] + col_c*x[2];
    return true;
  }

  bool Optimizer::eval_jac_g(Index n, const Number* x, bool new_x,
                             Index m, Index nele_jac, Index* iRow,
                             Index *jCol, Number* values)
  {
    if(values==NULL) {
      int nnz=0;
      for(int i=0; i<cn_; ++i) {
        iRow[nnz]=i; jCol[nnz]=0; ++nnz;
        iRow[nnz]=i; jCol[nnz]=1; ++nnz;
        iRow[nnz]=i; jCol[nnz]=2; ++nnz;
      }
    }
    else {
      Number *ptr=values;
      for(int i=0; i<cn_; ++i) {
        *ptr=vec_a_[i]; ++ptr;
        *ptr=vec_b_[i]; ++ptr;
        *ptr=vec_c_[i]; ++ptr;
      }
    }
    return true;
  }

  bool Optimizer::eval_h(Index n, const Number* x, bool new_x,
                         Number obj_factor, Index m, const Number* lambda,
                         bool new_lambda, Index nele_hess,
                         Index* iRow, Index* jCol, Number* values)
  {
    // hessian is all zero
    return true;
  }

  void Optimizer::finalize_solution(SolverReturn status,
                                    Index n, const Number* x, const Number* z_L, const Number* z_U,
                                    Index m, const Number* g, const Number* lambda,
                                    Number obj_value,
                                    const IpoptData* ip_data,
                                    IpoptCalculatedQuantities* ip_cq)
  {
    if(status == SUCCESS) {
      std::cout << "success!\n x:"  << Eigen::Map<const Eigen::Matrix<Number,1,3>>(x) << std::endl;
      std::cout << "g:" << Eigen::Map<const Eigen::Matrix<Number,1,Eigen::Dynamic>>(g, cn_) << std::endl;
      copy(x, x+3, res_x_);
    } else {
      std::cerr << "failed!" << std::endl;
    }
  }

#ifdef ZSW_DEBUG
  bool Optimizer::verify() const
  {
    assert(cn_==vec_a_.size());
    assert(cn_==vec_b_.size());
    assert(cn_==vec_c_.size());
    assert(cn_==vec_d_.size());

    const double eps=1e-5;
    for(int i=0; i<cn_; ++i) {
      if(vec_a_[i]*res_x_[0]+vec_b_[i]*res_x_[1]+vec_c_[i]*res_x_[2] < vec_d_[i]-eps) {
        std::cerr << "ERROR : " << vec_a_[i]*res_x_[0]+vec_b_[i]*res_x_[1]+vec_c_[i]*res_x_[2]-vec_d_[i] << std::endl;
        return false;
      }
    }
    return true;
  }

  bool Optimizer::verify(const Eigen::Matrix<Ipopt::Number,3,1> &res) const
  {
    assert(cn_==vec_a_.size());
    assert(cn_==vec_b_.size());
    assert(cn_==vec_c_.size());
    assert(cn_==vec_d_.size());

    const double eps=1e-5;
    for(int i=0; i<cn_; ++i) {
      if(vec_a_[i]*res[0]+vec_b_[i]*res[1]+vec_c_[i]*res[2] < vec_d_[i]-eps) {
        std::cerr << "ERROR : " << vec_a_[i]*res[0]+vec_b_[i]*res[1]+vec_c_[i]*res[2]-vec_d_[i] << std::endl;
        return false;
      }
    }
    return true;
  }

#endif
}

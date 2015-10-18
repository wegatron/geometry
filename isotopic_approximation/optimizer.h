#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <cassert>
#include <iostream>

#include <Eigen/Dense>

#define HAVE_CSTDDEF
#include <Ipopt/coin/IpTNLP.hpp>
#undef HAVE_CSTDDEF
#include <Ipopt/coin/IpIpoptApplication.hpp>

namespace zsw
{
  class Optimizer: public Ipopt::TNLP
  {
  public:
    Optimizer(Eigen::Matrix<Number,3,1> cx);

    // ax+by+cz>=d
    void addConstraint(Number a, Number b, Number c, Number d);

    bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                      Index& nnz_h_lag, IndexStyleEnum& index_style);

    bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                         Index m, Number* g_l, Number* g_u);

    bool get_starting_point(Index n, bool init_x, Number* x,
                            bool init_z, Number* z_L, Number* z_U,
                            Index m, bool init_lambda,
                            Number* lambda);

    bool eval_f(Index n, const Number* x, bool new_x,
                Number& obj_value);

    bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

    bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

    bool eval_jac_g(Index n, const Number* x, bool new_x,
                    Index m, Index nele_jac, Index* iRow,
                    Index *jCol, Number* values);

    bool eval_h(Index n, const Number* x, bool new_x,
                Number obj_factor, Index m, const Number* lambda,
                bool new_lambda, Index nele_hess,
                Index* iRow, Index* jCol, Number* values);

    void finalize_solution(SolverReturn status,
                           Index n, const Number* x, const Number* z_L, const Number* z_U,
                           Index m, const Number* g, const Number* lambda,
                           Number obj_value,
                           const IpoptData* ip_data,
                           IpoptCalculatedQuantities* ip_cq);
  private:
    Ipopt::Index cn_; // number of constraint
    Eigen::Matrix<Number,3,1> cx_;
    vector<Number> vec_a_;
    vector<Number> vec_b_;
    vector<Number> vec_c_;
    vector<Number> vec_d_;
  };

}

#endif /* OPTIMIZER_H */

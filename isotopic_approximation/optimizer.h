#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <cassert>
#include <iostream>

#include <Eigen/Dense>

#define HAVE_CSTDDEF
#include <Ipopt/coin/IpTNLP.hpp>
#undef HAVE_CSTDDEF
#include <Ipopt/coin/IpIpoptApplication.hpp>

#define ZSW_DEBUG

namespace zsw
{
  class Optimizer: public Ipopt::TNLP
  {
  public:
    Optimizer(Eigen::Matrix<Ipopt::Number,3,1> cx);

    // ax+by+cz>=d
    void addConstraint(Ipopt::Number a, Ipopt::Number b, Ipopt::Number c, Ipopt::Number d);

    bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
                      Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style);

    bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                         Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u);

    bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
                            bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                            Ipopt::Index m, bool init_lambda,
                            Ipopt::Number* lambda);

    bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                Ipopt::Number& obj_value);

    bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f);

    bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g);

    bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                    Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
                    Ipopt::Index *jCol, Ipopt::Number* values);

    bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
                bool new_lambda, Ipopt::Index nele_hess,
                Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values);

    void finalize_solution(Ipopt::SolverReturn status,
                           Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L, const Ipopt::Number* z_U,
                           Ipopt::Index m, const Ipopt::Number* g, const Ipopt::Number* lambda,
                           Ipopt::Number obj_value,
                           const Ipopt::IpoptData* ip_data,
                           Ipopt::IpoptCalculatedQuantities* ip_cq);
    template<typename SCALAR>
      void getResult(Eigen::Matrix<SCALAR,3,1> &res) const { std::copy(&res_x_[0], &res_x_[0]+3, res.data()); }
#ifdef ZSW_DEBUG
    bool verify() const;
    bool verify(const Eigen::Matrix<Ipopt::Number,3,1> &res) const ;
#endif
  private:
    Ipopt::Index cn_; // number of constraint
    Ipopt::Number res_x_[3];
    Eigen::Matrix<Ipopt::Number,3,1> cx_;
    std::vector<Ipopt::Number> vec_a_;
    std::vector<Ipopt::Number> vec_b_;
    std::vector<Ipopt::Number> vec_c_;
    std::vector<Ipopt::Number> vec_d_;
  };

}

#endif /* OPTIMIZER_H */

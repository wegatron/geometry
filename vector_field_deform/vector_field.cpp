#include "vector_field.h"

#include <iostream>

#include <Eigen/Dense>

using namespace zsw;
using namespace std;

zsw::VectorField::VectorField()
{
}

void zsw::VectorField::val(const double *x, double *val)
{
  assert(ex_func_!=nullptr);
  assert(fx_func_!=nullptr);
  assert(br_func_!=nullptr);
  assert(rx_func_!=nullptr);
  // judge region and caculate vector field
  zsw::RegionFunc::REGION_TYPE region_type = rx_func_->judgeRegion(x);
  if( region_type == zsw::RegionFunc::OUTER_REGION) {
    fill(val, val+3, 0.0);
  } else if(region_type == zsw::RegionFunc::INNER_REGION) {
    Eigen::Vector3d jac_e, jac_f, res;
    ex_func_->jac(x, jac_e.data());
    fx_func_->jac(x, jac_f.data());
    res = jac_e.cross(jac_f);
    std::copy(res.data(), res.data()+3, val);
  } else {
    Eigen::Vector3d jac_e, jac_f, jac_b, res;
    ex_func_->jac(x, jac_e.data());
    fx_func_->jac(x, jac_f.data());
    rx_func_->jac(x, jac_b.data());
    double g_b_r = 0.0, r = rx_func_->val(x);
    br_func_->jac(&r, &g_b_r);
    jac_b *= g_b_r;

    Eigen::Vector3d pv=-jac_b*ex_func_->val(x)  +(1-br_func_->val(&r))*jac_e;
    Eigen::Vector3d qv = -jac_b*fx_func_->val(x)+ (1-br_func_->val(&r))*jac_f;
    res = pv.cross(qv);
    std::copy(res.data(), res.data()+3, val);
  }
}

void zsw::VectorField::setExFunc(std::shared_ptr<zsw::Function> ex_func)
{
  ex_func_ = ex_func;
}


void zsw::VectorField::setFxFunc(std::shared_ptr<zsw::Function> fx_func)
{
  fx_func_ = fx_func;
}

void zsw::VectorField::setBrFunc(std::shared_ptr<zsw::BlendFunc> br_func)
{
  br_func_ = br_func;
}

void zsw::VectorField::setRxFunc(std::shared_ptr<zsw::RegionFunc> rx_func)
{
  rx_func_ = rx_func;
}

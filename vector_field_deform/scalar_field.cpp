#include "scalar_field.h"

#include <iostream>
#include <stdexcept>
#include <string>
#include <string.h>
#include <assert.h>
#include <cmath>

extern "C" {
  void quad_scalar_field_(double *val, const double *x, const double *a, const double *c);
  void quad_scalar_field_jac_(double *jac, const double *x, const double *a, const double *c);
  void quad_scalar_field2_(double *val, const double *x, const double *a, const double *c);
  void quad_scalar_field2_jac_(double *jac, const double *x, const double *a, const double *c);
}

zsw::LinearScalarField::LinearScalarField(const double *u, const double *c)
{
  memcpy(u_, u, sizeof(double)*3);
  memcpy(c_, c, sizeof(double)*3);
}

double zsw::LinearScalarField::val(const double *x)
{
  return u_[0]*(x[0]-c_[0])+u_[1]*(x[1]-c_[1])+u_[2]*(x[2]-c_[2]);
}

void zsw::LinearScalarField::jac(const double *x, double *g)
{
  memcpy(g, u_, sizeof(double)*3);
}

zsw::LinearScalarField::~LinearScalarField(){}


zsw::QuadraticScalarField::QuadraticScalarField(const double *a, const double *c)
{
  assert(a!=nullptr && c!=nullptr);
  std::copy(a, a+3, a_);
  std::copy(c, c+3, c_);
}

double zsw::QuadraticScalarField::val(const double *x)
{
  assert(x != nullptr);
  double ret = 0.0;
  quad_scalar_field_(&ret, x, a_, c_);
  return ret;
}

void zsw::QuadraticScalarField::jac(const double *x, double *g)
{
  assert(x!=nullptr && g!=nullptr);
  quad_scalar_field_jac_(g, x, a_, c_);
}

zsw::QuadraticScalarField2::QuadraticScalarField2(const double *a, const double *center)
{
  std::copy(a, a+3, a_);
  std::copy(center, center+3, center_);
}

double zsw::QuadraticScalarField2::val(const double *x)
{
  double ret = 0.0;
  quad_scalar_field2_(&ret, x, a_, center_);
  return ret;
}

void zsw::QuadraticScalarField2::jac(const double *x, double *g)
{
  quad_scalar_field2_jac_(g, x, a_, center_);
}

zsw::BlendFunc::BlendFunc(const double ri, const double ro)
{
  ri_ = ri;
  ro_ = ro;
  if(ri_ > ro_) {
    std::cerr << "ro < ri " << ro_  << " " << ri_ << std::endl;
    throw std::logic_error("ro < ri "+std::to_string(ro_) + " " + std::to_string(ri_));
  }
}

double zsw::BlendFunc::val(const double *r)
{
  if(*r<ri_ || *r>ro_) {
    std::cerr << "[ERROR] r in blend func should in (ri,ro):" << ri_ << " " << ro_ <<  " " << *r << std::endl;
    throw std::logic_error("[ERROR] r in blend func should in (ri,ro):"+std::to_string(ri_)+" "+std::to_string(ro_));
  }
  double v= (*r-ri_)/(ro_-ri_);
  return 4*v*v*v*(1-v) + v*v*v*v;
}

void zsw::BlendFunc::jac(const double *r, double *g)
{
  g[0] = 12*(*r-ri_)*(*r-ri_)*(ro_-*r)/((ro_-ri_)*(ro_-ri_)*(ro_-ri_)*(ro_-ri_));
}

zsw::BlendFunc::~BlendFunc() {}

zsw::SphereRegionFunc::SphereRegionFunc(const double ri, const double ro, const double *center)
{
  ri_ = ri;
  ro_ = ro;
  if(ri_ > ro_) {
    std::cerr << "ro < ri " << ro_  << " " << ri_ << std::endl;
    throw std::logic_error("ro < ri "+std::to_string(ro_) + " " + std::to_string(ri_));
  }
  memcpy(center_, center, sizeof(double)*3);
}

double zsw::SphereRegionFunc::val(const double *x)
{
  return  sqrt((x[0]-center_[0])*(x[0]-center_[0])+(x[1]-center_[1])*(x[1]-center_[1])+(x[2]-center_[2])*(x[2]-center_[2]));
}

void zsw::SphereRegionFunc::jac(const double *x, double *g)
{
  double v_tmp = sqrt((x[0]-center_[0])*(x[0]-center_[0]) +
                          (x[1]-center_[1])*(x[1]-center_[1]) + (x[2]-center_[2])*(x[2]-center_[2]));
  // std::cerr << "v_tmp:" << v_tmp << std::endl;
  // std::cerr << "x[0]-center_[0]:" << x[0]-center_[0]<< std::endl;
  g[0] = (x[0]-center_[0])/v_tmp;   g[1] = (x[1]-center_[1])/v_tmp;   g[2] = (x[2]-center_[2])/v_tmp;
}

zsw::RegionFunc::REGION_TYPE zsw::SphereRegionFunc::judgeRegion(const double *x)
{
  double v = (x[0]-center_[0])*(x[0]-center_[0])+(x[1]-center_[1])*(x[1]-center_[1])+(x[2]-center_[2])*(x[2]-center_[2]);
  // std::cerr <<"[DEBUG]" <<  v << " " << ro_*ro_ << " " << ri_*ri_ << std::endl;
  if(v >= ro_*ro_)   {
    return zsw::RegionFunc::OUTER_REGION;
  }
  if(v < ri_*ri_) {
    return zsw::RegionFunc::INNER_REGION;
  }
  return zsw::RegionFunc::BLENDER_REGION;
}

zsw::IsosurfacesRegionFunc::IsosurfacesRegionFunc(const double *b, const double *center, const double ri, const double ro)
{
  std::copy(b, b+3, b_);
  std::copy(center, center+3, center_);
  ro_ = ro;
  ri_ = ri;
}
double zsw::IsosurfacesRegionFunc::val(const double*x)
{
  return (x[0]-center_[0])*b_[0] + (x[1]-center_[1])*b_[1] + (x[2]-center_[2])*b_[2];
}
void zsw::IsosurfacesRegionFunc::jac(const double *x, double *g)
{
  std::copy(b_, b_+3, g);
}

zsw::RegionFunc::REGION_TYPE zsw::IsosurfacesRegionFunc::judgeRegion(const double *x)
{
  double v = val(x);
  if(v<ri_) {
    return zsw::RegionFunc::INNER_REGION;
  } else if(v>=ro_) {
    return zsw::RegionFunc::OUTER_REGION;
  }
  return zsw::RegionFunc::BLENDER_REGION;
}

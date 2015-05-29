#include <iostream>

#include <Eigen/Dense>
#include <chrono>

#include "../vector_field.h"

using namespace std;
using namespace zsw;

#define EPS 1e-5

void testValSpecific()
{
  bool suc = true;

  // simplest case 0 : inner region
  {
    VectorField vf;
    Eigen::Vector3d u[2], c;
    u[0] << 1,0,0; u[1] << 0,1,0;  c << 0,0,0;
    std::shared_ptr<Function> ex_func(new LinearScalarField(u[0].data(), c.data()));
    std::shared_ptr<Function> fx_func(new LinearScalarField(u[1].data(), c.data()));

    double r[2] = {1.0, 2.0};
    std::shared_ptr<BlendFunc> br_func(new BlendFunc(r[0], r[1]));
    std::shared_ptr<RegionFunc> rx_func(new SphereRegionFunc(r[0], r[1], c.data()));
    vf.setExFunc(ex_func);
    vf.setFxFunc(fx_func);
    vf.setBrFunc(br_func);
    vf.setRxFunc(rx_func);

    // test val
    Eigen::Vector3d val_vf, x, val_expected;
    x<< 0.7, 0.7, 0; val_expected << 0, 0, 1;

    vf.val(x.data(), val_vf.data());
    if((val_vf - val_expected).squaredNorm() > EPS) {
      suc = false;
      std::cerr << "[ERROR]" << __FILE__ << " line: " <<  __LINE__  << std::endl;
      std::cerr << val_vf.transpose() << std::endl;
    }
  }

  // simplest case 1: outer region
  {
    VectorField vf;
    Eigen::Vector3d u[2], c;
    c << -14.13, 43.2, -6.3;
    u[0] << -3.34, -8.3, 6.7;
    u[1] << 1.3, 5.8, -7.833134;
    std::shared_ptr<Function> ex_func(new LinearScalarField(u[0].data(), c.data()));
    std::shared_ptr<Function> fx_func(new LinearScalarField(u[1].data(), c.data()));

    double r[2] = {4.3, 6.9};
    std::shared_ptr<BlendFunc> br_func(new BlendFunc(r[0], r[1]));
    std::shared_ptr<RegionFunc> rx_func(new SphereRegionFunc(r[0], r[1], c.data()));
    vf.setExFunc(ex_func);
    vf.setFxFunc(fx_func);
    vf.setBrFunc(br_func);
    vf.setRxFunc(rx_func);

    // test val
    Eigen::Vector3d val_vf, x;
    x << -10.3, 40.2, -1.3;
    vf.val(x.data(), val_vf.data());
    if(val_vf.squaredNorm() > EPS) {
      suc = false;
      std::cerr << "[ERROR]" << __FILE__ << " line: " << __LINE__  << std::endl;
    }
  }

  // simplest case 2: blend region
  {
    VectorField vf;
    Eigen::Vector3d u[2], c;
    c << -14.13, 43.2, -6.3;
    u[0] << -3.34, -8.3, 6.7;
    u[1] << 1.3, 5.8, -7.833134;
    std::shared_ptr<Function> ex_func(new LinearScalarField(u[0].data(), c.data()));
    std::shared_ptr<Function> fx_func(new LinearScalarField(u[1].data(), c.data()));

    double r[2] = {4.3, 6.9};
    std::shared_ptr<BlendFunc> br_func(new BlendFunc(r[0], r[1]));
    std::shared_ptr<RegionFunc> rx_func(new SphereRegionFunc(r[0], r[1], c.data()));
    vf.setExFunc(ex_func);
    vf.setFxFunc(fx_func);
    vf.setBrFunc(br_func);
    vf.setRxFunc(rx_func);

    // test val
    Eigen::Vector3d val_vf, x, val_expected, jac_p, jac_q;
    x << -10.3, 42.1, -4.1;
    jac_p << -3.69317, -8.16737,  6.46809;
    jac_q << 1.90862,   5.6044, -7.45463;

    val_expected = jac_p.cross(jac_q);
    vf.val(x.data(), val_vf.data());
    if((val_vf-val_expected).squaredNorm() > EPS) {
      suc = false;
      std::cerr << "[ERROR]" << __FILE__ << " line: " << __LINE__  << std::endl;
      std::cerr <<  "val_expected:" << val_expected.transpose() << std::endl;
      std::cerr << "val_get:" << val_vf.transpose() << std::endl;
    }
  }

  if(suc) { std::cout << "[INFO]" << __FUNCTION__ << "passed!" << std::endl; }
  suc = false;
}

static VectorField* genRandomVectorField(double *x)
{
  VectorField *vf = new VectorField();
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> distribution(-100.0,100.0);
  std::uniform_real_distribution<double> distribution_ri(0.0,100.0);
  std::uniform_real_distribution<double> distribution_ro(100.0,300.0);
  Eigen::Vector3d u[2], c;
  double r[2] = {distribution_ri(generator), distribution_ro(generator)};
  for(size_t i=0; i<3; ++i) {
    u[0][i] = distribution(generator);
    c[i] = distribution(generator);
    x[i] = distribution(generator);
  }
  u[1][0] = -u[0][1]; u[1][1] = u[0][0]; u[1][2] = 0;
  // std::cout << "u[0]:" << u[0].transpose() << std::endl;
  // std::cout << "u[1]:" << u[1].transpose() << std::endl;
  // std::cout << "c:" << c.transpose() << std::endl;
  // std::cout << "x:" << x[0] << " " << x[1] << " " << x[2] << std::endl;
  // std::cout << "r:" << r[0] << "  " << r[1] << std::endl;
  std::shared_ptr<Function> ex_func(new LinearScalarField(u[0].data(), c.data()));
  std::shared_ptr<Function> fx_func(new LinearScalarField(u[1].data(), c.data()));

  std::shared_ptr<BlendFunc> br_func(new BlendFunc(r[0], r[1]));
  std::shared_ptr<RegionFunc> rx_func(new SphereRegionFunc(r[0], r[1], c.data()));
  vf->setExFunc(ex_func);
  vf->setFxFunc(fx_func);
  vf->setBrFunc(br_func);
  vf->setRxFunc(rx_func);
  return vf;
}

static void cacuValExpected(std::shared_ptr<VectorField> vf, Eigen::Vector3d &x, Eigen::Vector3d &val_expected)
{
  Eigen::Vector3d pv, qv;
  std::shared_ptr<Function> ex_func = vf->getExFunc();
  std::shared_ptr<Function> fx_func = vf->getFxFunc();
  std::shared_ptr<RegionFunc> rx_func = vf->getRxFunc();
  std::shared_ptr<BlendFunc> br_func = vf->getBrFunc();
  if(rx_func->judgeRegion(x.data()) == zsw::RegionFunc::INNER_REGION) {
    std::cerr << "InnerRegion" << std::endl;
    ex_func->jac(x.data(), pv.data());
    fx_func->jac(x.data(), qv.data());
  } else if(rx_func->judgeRegion(x.data()) == zsw::RegionFunc::OUTER_REGION) {
    std::cerr << "OutRegion" << std::endl;
    pv.setZero();
    qv.setZero();
  } else {
    std::cerr << "BlendRegion" << std::endl;
    const double eps = 1e-6;
    for(size_t i=0; i<3; ++i) {
      double save = x[i];
      double p_v[4];
      double q_v[4];
      x[i] = save+eps;
      double r = rx_func->val(x.data());
      p_v[0] = (1-br_func->val(&r))*ex_func->val(x.data());
      q_v[0] = (1-br_func->val(&r))*fx_func->val(x.data());
      x[i] = save - eps;
      r = rx_func->val(x.data());
      p_v[1] = (1-br_func->val(&r))*ex_func->val(x.data());
      q_v[1] = (1-br_func->val(&r))*fx_func->val(x.data());
      x[i] = save+2*eps;
      r = rx_func->val(x.data());
      p_v[2] = (1-br_func->val(&r))*ex_func->val(x.data());
      q_v[2] = (1-br_func->val(&r))*fx_func->val(x.data());
      x[i] = save - 2*eps;
      r = rx_func->val(x.data());
      p_v[3] = (1-br_func->val(&r))*ex_func->val(x.data());
      q_v[3] = (1-br_func->val(&r))*fx_func->val(x.data());
      pv[i] -= (8*(p_v[0]-p_v[1])-p_v[2]+p_v[3])/(12*eps);
      qv[i] -= (8*(q_v[0]-q_v[1])-q_v[2]+q_v[3])/(12*eps);
    }
  }
  val_expected = pv.cross(qv);
}

void testValRandom(size_t times)
{
  bool suc = true;
  do {
    Eigen::Vector3d x, val_vf, val_expected;
    std::shared_ptr<VectorField> vf(genRandomVectorField(x.data()));
    vf->val(x.data(), val_vf.data());
    cacuValExpected(vf, x, val_expected);
    if((val_vf-val_expected).squaredNorm() > EPS) {
      suc = false;
      std::cerr << "[ERROR]" << __FILE__ << " line: " << __LINE__  << std::endl;
      std::cerr <<  "val_expected:" << val_expected.transpose() << std::endl;
      std::cerr << "val_get:" << val_vf.transpose() << std::endl;
      break;
    }
  } while(--times);

  if(suc) { std::cout << "[INFO]" << __FUNCTION__ << "passed!" << std::endl; }
}

int main(int argc, char *argv[])
{
  testValSpecific();
  testValRandom(100);
  return 0;
}

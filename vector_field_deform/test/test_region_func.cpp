#include <iostream>
#include <chrono>
#include <memory>
#include <Eigen/Dense>

#include "../scalar_field.h"
#include "jac_hes_err.h"

using namespace std;
using namespace zsw;

#define EPS 1e-6

void testValSpecific()
{
  bool suc = true;

  // group 1
  {
    double c[3] = {0,0,0};
    SphereRegionFunc spf(10, 35, c);
    double x[3] = {3.2,-1.4, 2.7};
    if( fabs(spf.val(x) -  4.414748010928823)> EPS) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      std::cerr << "val:" << spf.val(x) << std::endl;
      suc = false;
    }
    if(spf.judgeRegion(x) != zsw::RegionFunc::INNER_REGION) {
      suc =false;
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
    }
  }

  // group 2
  {
    double c[3] = {87.2,45.3,92.4};
    SphereRegionFunc spf(100, 120, c);
    double x[3] = {3.2,-1.4, 2.7};
    if( fabs(spf.val(x) -  131.4647481266365)> EPS) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      std::cerr << "val:" << spf.val(x) << std::endl;
      suc = false;
    }
    if(spf.judgeRegion(x) != zsw::RegionFunc::OUTER_REGION) {
      suc = false;
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
    }
  }

  //boundary ri
  {
    double c[3] = {67.182, 87.12, -123.12};
    SphereRegionFunc spf(92.86, 126, c);
    double x[3] = {1.2, 30.4, -90.7};
    if(fabs(spf.val(x) - 92.85385896127312)>EPS) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      std::cerr << "val:" << spf.val(x) << std::endl;
      suc = false;
    }
    if(spf.judgeRegion(x) != zsw::RegionFunc::INNER_REGION) {
      suc = false;
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
    }
  }

  //boundary ri
  {
    double c[3] = {67.182, 87.12, -123.12};
    SphereRegionFunc spf(92.86, 126, c);
    double x[3] = {0, 30.4, -90.7};
    if(fabs(spf.val(x) - 93.71038322405901)>EPS) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      std::cerr << "val:" << spf.val(x) << std::endl;
      suc = false;
    }
    if(spf.judgeRegion(x) != zsw::RegionFunc::BLENDER_REGION) {
      suc = false;
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
    }
  }

  //boundary ro
  {
    double c[3] = {62.19, -26.2, 57.0};
    SphereRegionFunc spf(92.85, 165, c);
    double x[3] = {11.2, 25.4, -90.7};
    if(fabs(spf.val(x) - 164.553426278519)>EPS) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      std::cerr << "val:" << spf.val(x) << std::endl;
      suc = false;
    }
    if(spf.judgeRegion(x) != zsw::RegionFunc::BLENDER_REGION) {
      suc = false;
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
    }
  }

  //boundary ro
  {
    double c[3] = {62.19, -26.2, 57.0};
    SphereRegionFunc spf(92.85, 165, c);
    double x[3] = {11.2, 25.4, -91.7};
    if(fabs(spf.val(x) - 165.4515944317249)>EPS) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      std::cerr << "val:" << spf.val(x) << std::endl;
      suc = false;
    }
    if(spf.judgeRegion(x) != zsw::RegionFunc::OUTER_REGION) {
      suc = false;
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
    }
  }

  if(suc) { std::cout << "[INFO]" << __FUNCTION__ << "passed!" << std::endl; }
}

static SphereRegionFunc* genRandomFunc(double *x, double *c, double *r)
{
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> distribution(-100.0,100.0);
  std::uniform_real_distribution<double> distribution_ro(150.0,200.0);
  std::uniform_real_distribution<double> distribution_ri(10.0,120.0);
  for(size_t i=0; i<3; ++i) {
    x[i] = distribution(generator);
    c[i] = distribution(generator);
  }
  r[0] = distribution_ri(generator);
  r[1] = distribution_ro(generator);
  return new SphereRegionFunc(r[0], r[1], c);
}

void testValRandom(size_t time)
{
  bool suc = true;
  Eigen::Vector3d x(3), c(3);
  double r[2];
  do {
    std::shared_ptr<SphereRegionFunc> spf(genRandomFunc(x.data(), c.data(), r));
    if(fabs(spf->val(x.data()) - (x-c).norm())  > EPS) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      std::cerr << "c: " << c.transpose() << std::endl;
      std::cerr << "x:" << x.transpose() << std::endl;
      suc = false;
      break;
    }
    zsw::RegionFunc::REGION_TYPE region_type = zsw::RegionFunc::BLENDER_REGION;
    if((x-c).norm() > r[1]) { region_type = zsw::RegionFunc::OUTER_REGION; }
    else if((x-c).norm() < r[0] ) { region_type = zsw::RegionFunc::INNER_REGION; }
    if(spf->judgeRegion(x.data()) != region_type) {
      suc = false;
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
    }
  } while(--time);
  if(suc) { std::cout << "[INFO]"  << __FUNCTION__ << " passed!" << std::endl; }
}

void testGraErr(size_t time)
{
  bool suc = true;
  Eigen::VectorXd x(3), c(3);
  double r[2];
  do{
    std::shared_ptr<SphereRegionFunc> spf(genRandomFunc(x.data(), c.data(), r));
    double max_err = graErr(*spf, x);
    if(max_err > 1e-3) {
      suc = false;
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      std::cerr << "max_err:" << max_err << std::endl;
      std::cerr << "x:" << x.transpose() << std::endl;
      std::cerr << "c:" << c.transpose() << std::endl;
      Eigen::VectorXd g(3); spf->jac(x.data(), g.data());
      std::cerr << "jac:" << g.transpose() << std::endl;
      break;
    }
  }while(--time);

  if(suc) { std::cout << "[INFO]"  << __FUNCTION__ << " passed!" << std::endl; }
}

int main(int argc, char *argv[])
{
  testValSpecific();
  testValRandom(100);
  testGraErr(100);
  return 0;
}

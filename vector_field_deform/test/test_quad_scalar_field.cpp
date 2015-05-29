#include <iostream>
#include <vector>
#include <random>

#include <Eigen/Dense>
#include <chrono>

#include "../scalar_field.h"

#include "jac_hes_err.h"

#define EPS 1e-6

using namespace std;
using namespace zsw;


void testValSpecific()
{
  bool suc = true;
  Eigen::Vector3d a, x, c;
  {
    a << 4.25, 4.52, 3.46;
    x << 6.54,-2.3, 1.2;
    c << 7.8, 8, -2.3;
    QuadraticScalarField quad_scalar_field(a.data(), c.data());
    if(fabs(quad_scalar_field.val(x.data()) - 4467.9667692) > EPS) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      suc = false;
    }
  }
  if(suc) { std::cout << "[INFO]" << __FUNCTION__ << "passed!" << std::endl; }
}

void testValRandom(size_t time)
{
  bool suc = true;
  Eigen::Vector3d a, x, c;
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_real_distribution<double> distribution(-10.0,10.0);
  do {
    for(size_t i=0; i<3; ++i) {
      a[i] = distribution(generator);
      x[i] = distribution(generator);
      c[i] = distribution(generator);
    }
    QuadraticScalarField qfs(a.data(), c.data());
    double val_expected = (a.cross(x-c)).squaredNorm();
    if(fabs(qfs.val(x.data())- val_expected)  > EPS) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      std::cerr <<  "x: " << x.transpose() << std::endl;
      std::cerr << "a: " << a.transpose() << std::endl;
      std::cerr << "c:" << c.transpose() << std::endl;
      suc = false;
      break;
    }
  } while(--time);
  if(suc) { std::cout << "[INFO]"  << __FUNCTION__ << " passed!" << std::endl; }
}


void testGraErr(size_t time)
{
  bool suc = true;
  Eigen::VectorXd x(3), a(3), c(3), g(3);

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_real_distribution<double> distribution(-100.0,100.0);
  do {
    for(size_t i=0; i<3; ++i) {
      x[i] = distribution(generator);
      a[i] = distribution(generator);
      c[i] = distribution(generator);
    }
    QuadraticScalarField qfs(a.data(), c.data());
    double max_err = graErr(qfs, x);
    if(max_err > 0.1) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      std::cerr << "err val:" << graErr(qfs, x) << std::endl;
      std::cerr <<  "x: " << x << std::endl;
      std::cerr << "a: " << a << std::endl;
      std::cerr << "c:" << c << std::endl;
      suc = false;
      break;
    }
  } while(--time);
  if(suc) { std::cout << "[INFO]"  << __FUNCTION__ << " passed!" << std::endl; }
}


int main(int argc, char *argv[])
{
  testValSpecific();
  testValRandom(100);
  testGraErr(100);
  return 0;
}

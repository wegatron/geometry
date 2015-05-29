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
  double a[3] = {1.4, -23.1, -20.1};
  double c[3] = {12.2, 14.5, -87.1};
  double x[3] = {0.23, 4.5, -0.2};
  QuadraticScalarField2 qsf2(a, c);
  if(fabs(qsf2.val(x) - 2348396.872704) >EPS) {
    std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
    suc = false;
  }
  if(suc) { std::cout << "[INFO]" << __FUNCTION__ << "passed!" << std::endl; }
}

void testValRandom(size_t time)
{
  bool suc = true;
  Eigen::Vector3d a, c, x;
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_real_distribution<double> distribution(-10.0,10.0);
  do{
    for(size_t i=0; i<3; ++i) {
      a[i] = distribution(generator);
      c[i] = distribution(generator);
      x[i] = distribution(generator);
    }
    QuadraticScalarField2 qsf2(a.data(), c.data());
    double v_expected = a.dot(x-c);
    v_expected *= v_expected;
    if(fabs(qsf2.val(x.data()) - v_expected) > EPS) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      suc = false;
      break;
    }
  }while(--time);
  if(suc) { std::cout << "[INFO]" << __FUNCTION__ << "passed!" << std::endl; }
}

void testGarErr(size_t time)
{
  bool suc = true;
  Eigen::VectorXd a(3), c(3), x(3);
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_real_distribution<double> distribution(-100.0,100.0);
  do{
    for(size_t i=0; i<3; ++i) {
      a[i] = distribution(generator);
      c[i] = distribution(generator);
      x[i] = distribution(generator);
    }
    QuadraticScalarField2 qsf2(a.data(), c.data());
    double max_err = graErr(qsf2, x);
    if(max_err > 0.5) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      std::cerr << "err val:" << graErr(qsf2, x) << std::endl;
      std::cerr <<  "x: " << x << std::endl;
      std::cerr << "a: " << a << std::endl;
      std::cerr << "c:" << c << std::endl;
      suc = false;
      break;
    }
  } while(--time);
  if(suc) { std::cout << "[INFO]" << __FUNCTION__ << "passed!" << std::endl; }
}

int main(int argc, char *argv[])
{
  testValSpecific();
  testValRandom(100);
  testGarErr(100);
  return 0;
}

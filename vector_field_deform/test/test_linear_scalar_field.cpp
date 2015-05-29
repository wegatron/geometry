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

void testValSpecafic()
{
  Eigen::Vector3d x, u, c;
  bool suc = true;
  { // group 1
    x<< 1.1, 1.3, 1.8;
    u<< 0.5, 0.6, 0.7;
    c<< 2.4, 1.2, 4.1;
    LinearScalarField lsf(u.data(), c.data());
    if( fabs(lsf.val(x.data()) +2.2)>EPS ) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      suc = false;
    }
  }
  { // group 2
    x<< 10.12, 24.4, 45.8;
    u<< 23.923, 123.123, 193.0;
    c<< 923.0, 123.0, 52.23;
    LinearScalarField lsf(u.data(), c.data());
    if( fabs(lsf.val(x.data()) + 35219.7460400000)>EPS ) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      suc = false;
    }
  }
  if(suc) { std::cout << "[INFO] " << __FUNCTION__ << " passed!" << std::endl; }
}

void testValRandom( size_t time)
{
  bool suc = true;
  Eigen::Vector3d x, u, c;

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_real_distribution<double> distribution(-1000.0,1000.0);
 do {
    for(size_t i=0; i<3; ++i) {
      x[i] = distribution(generator);
      u[i] = distribution(generator);
      c[i] = distribution(generator);
    }
    LinearScalarField lfs(u.data(), c.data());
    if(fabs(lfs.val(x.data()) - u.dot(x-c)) > EPS) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      std::cerr <<  "x: " << x << std::endl;
      std::cerr << "u: " << u << std::endl;
      std::cerr << "c:" << c << std::endl;
      suc = false;
      break;
    }
 } while(--time);
  if(suc) { std::cout << "[INFO]"  << __FUNCTION__ << " passed!" << std::endl; }
}

void testGraErr(size_t time)
{
  bool suc = true;
  Eigen::VectorXd x(3), u(3), c(3), g(3);

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_real_distribution<double> distribution(-1000.0,1000.0);
  do {
    for(size_t i=0; i<3; ++i) {
      x[i] = distribution(generator);
      u[i] = distribution(generator);
      c[i] = distribution(generator);
    }
    LinearScalarField lfs(u.data(), c.data());
    double max_err = graErr(lfs, x);
    if(max_err > 5e-4) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      std::cerr << "err val:" << graErr(lfs, x) << std::endl;
      std::cerr <<  "x: " << x << std::endl;
      std::cerr << "u: " << u << std::endl;
      std::cerr << "c:" << c << std::endl;
      suc = false;
      break;
    }
  } while(--time);
  if(suc) { std::cout << "[INFO]"  << __FUNCTION__ << " passed!" << std::endl; }
}


int main(int argc, char *argv[])
{
  testValSpecafic();
  testValRandom(100);
  testGraErr(10);
  return 0;
}

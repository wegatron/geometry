#include "../scalar_field.h"

#include <iostream>
#include <chrono>
#include <cmath>
#include <Eigen/Dense>

#include "jac_hes_err.h"

using namespace std;
using namespace zsw;

#define EPS 1e-6

void testValSpecific()
{
  bool suc = true;
  double x;
  //group 1
  {
    BlendFunc blf(12.5,20.6);
    x = 15.1;
    if(fabs(blf.val(&x) - 0.100441936100081) > EPS) {
      std::cerr << "val:" << blf.val(&x) << std::endl;
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      suc = false;
    }
  }

  //group 1
  {
    BlendFunc blf(56.123,87.16);
    x = 67.777;
    if(fabs(blf.val(&x) - 0.152125767593332) > EPS) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
    }
  }

  //group 3 boundary 0
  {
    BlendFunc blf(26.123,47.16);
    x = 26.123;
    if(fabs(blf.val(&x)) > EPS) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
    }
  }

  //group 4 boundary 1
  {
    BlendFunc blf(26.123,47.16);
    x = 47.16;
    if(fabs(blf.val(&x)-1) > EPS) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
    }
  }
  if(suc) { std::cout << "[INFO]" << __FUNCTION__ << "passed!" << std::endl; }
}


void testValRandom(size_t time)
{
  bool suc = true;
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_real_distribution<double> distribution_ri(0.0,1000.0);
  std::uniform_real_distribution<double> distribution_r(1000.0,1500.0);
  std::uniform_real_distribution<double> distribution_ro(1500.0,3000.0);

  double x, r[2];
  do{
     r[0] = distribution_ri(generator);
     r[1] = distribution_ro(generator);
     x = distribution_r(generator);
     BlendFunc blf(r[0], r[1]);
     double ans = (4*pow((x-r[0]),3)*(r[1]-x)+pow(x-r[0],4)) / pow(r[1]-r[0],4);
     if(fabs(blf.val(&x) - ans) > EPS) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      std::cout << "ri, ro, x:" << r[0]  << " " << r[1] << " " << x << std::endl;
       suc = false;
     }
  }while(--time);
  if(suc) { std::cout << "[INFO]" << __FUNCTION__ << "passed!" << std::endl; }
}

void testGraErr(size_t time)
{
  bool suc = true;
  Eigen::VectorXd x(1);
  double r[2];

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_real_distribution<double> distribution_ri(0.0,1000.0);
  std::uniform_real_distribution<double> distribution_r(1000.0,1500.0);
  std::uniform_real_distribution<double> distribution_ro(1500.0,3000.0);
  do {
    x[0] = distribution_r(generator);
    r[0] = distribution_ri(generator);
    r[1] = distribution_ro(generator);
    BlendFunc blf(r[0], r[1]);
    double max_err = graErr(blf, x);
    if(max_err > 5e-4) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      std::cerr << "err val:" << max_err << std::endl;
      std::cerr <<  "x: " << x << std::endl;
      std::cerr << "ri, ro:" << r[0] << " " << r[1] << std::endl;
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

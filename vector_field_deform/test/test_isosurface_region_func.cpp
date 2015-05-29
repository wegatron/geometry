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
  {
    Eigen::Vector3d b,c,x;
    b << 34.5, 23.3, 0.0;
    c << 2.9, -6.8, 8.0;
    x << 45, 5, 12;
    double r[2] = {0.0, 23.0};
    IsosurfacesRegionFunc ifs_region(b.data(), c.data(), r[0], r[1]);
    if(fabs(ifs_region.val(x.data()) - 1727.39) > EPS) {
      std::cerr << "[ERROR]" << __FILE__ << __LINE__ << std::endl;
      suc = false;
    }
  }
  if(suc) { std::cout << "[INFO] " << __FUNCTION__ << " passed!" << std::endl; }
}

int main(int argc, char *argv[])
{
  testValSpecific();
  return 0;
}

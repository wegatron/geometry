#define BOOST_AUTO_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include "../triangulation_fix.h"

using namespace std;

BOOST_AUTO_TEST_SUITE(test_tet_quality)

BOOST_AUTO_TEST_CASE(test_regular_tet)
{
  Eigen::Matrix<zsw::Scalar,3,1> v0,v1,v2,v3;
  v0<<0,0,0; v1<<1,1,0;
  v2<<1,0,1; v3<<0,1,1;
  zsw::Scalar val = zsw::calcTetQuality(v0,v1,v2,v3);
  BOOST_CHECK_MESSAGE(fabs(val-1)<1e-3, "regular tet val:"+std::to_string(val));
}

BOOST_AUTO_TEST_CASE(test_val_diff0)
{
  Eigen::Matrix<zsw::Scalar,3,1> v0,v1,v2,v3;
  v0<<0,0,0; v1<<2,0,0;
  v2<<1,1,0; v3<<1,0.5,3;
  zsw::Scalar val0 = zsw::calcTetQuality(v0,v1,v2,v3);
  v3<<1,0.5,5;
  zsw::Scalar val1 = zsw::calcTetQuality(v0,v1,v2,v3);
  std::cerr << val0 << ":" << val1 << std::endl;
  BOOST_CHECK(val0>val1);
}

BOOST_AUTO_TEST_CASE(test_check_val)
{
  Eigen::Matrix<zsw::Scalar,3,1> v0,v1,v2,v3;
  v0<<0.588376, -0.010356, -1.178430;
  v1<<0.606244, -0.024685, -1.1918;
  v2<<0.596656, -0.034007, -1.18226;
  v3<<0.610707, -0.040254, -1.188310;
  zsw::Scalar val0 = zsw::calcTetQuality(v0,v1,v2,v3);
  std::cerr << "v0:" << v0.transpose() << std::endl;
  std::cerr << "v1:" << v1.transpose() << std::endl;
  std::cerr << "v2:" << v2.transpose() << std::endl;
  std::cerr << "v3:" << v3.transpose() << std::endl;
  std::cerr << val0 << std::endl;
}
BOOST_AUTO_TEST_SUITE_END()

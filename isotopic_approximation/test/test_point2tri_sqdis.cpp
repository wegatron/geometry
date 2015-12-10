#define BOOST_AUTO_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include "../basic_op.h"

using namespace std;

BOOST_AUTO_TEST_SUITE(test_point2tri_dis)

BOOST_AUTO_TEST_CASE(test_0)
{
  Eigen::Matrix<zsw::Scalar,3,1> v0,v1,v2,vt;
  v0 << 0,0,0;
  v1 << 2,0,0;
  v2 << 1,2,0;
  vt << 1,1,1;
  zsw::Scalar sq_dis = zsw::calcPoint2TriSquaredDis(vt, v0, v1,v2);
  BOOST_CHECK_MESSAGE(fabs(sq_dis-1)<1e-3, "Dis error is:"+std::to_string(fabs(sq_dis-1)));
}

BOOST_AUTO_TEST_SUITE_END()

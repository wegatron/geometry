#define BOOST_AUTO_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include "../triangulation2.h"

using namespace std;

BOOST_AUTO_TEST_SUITE(KernelRegionJudger)

BOOST_AUTO_TEST_CASE(single_constraint)
{
  zsw::KernelRegionJudger krj;
  Eigen::Matrix<zsw::Scalar,3,1> v0; v0 << 0,0,1;
  Eigen::Matrix<zsw::Scalar,3,1> v1; v1 << 1,0,0;
  Eigen::Matrix<zsw::Scalar,3,1> v2; v2 << 0,1,0;
  Eigen::Matrix<zsw::Scalar,3,1> vr;  vr << 1,1,1;
  krj.addConstraint(v0,v1,v2,vr);

  Eigen::Matrix<zsw::Scalar,3,1> vc0;
  vc0<<0.1,0.1,0.1;
  BOOST_CHECK(!krj.judge(vc0));

  Eigen::Matrix<zsw::Scalar,3,1> vc1;
  vc1<<1,0,1;
  BOOST_CHECK(krj.judge(vc1));

  Eigen::Matrix<zsw::Scalar,3,1> vc2;
  vc2<<-1,0,0;
  BOOST_CHECK(!krj.judge(vc2));
}

BOOST_AUTO_TEST_CASE(mutiple_constraint)
{
  zsw::KernelRegionJudger krj;
  Eigen::Matrix<zsw::Scalar,3,1> v0; v0 << 0,0,1;
  Eigen::Matrix<zsw::Scalar,3,1> v1; v1 << 1,0,0;
  Eigen::Matrix<zsw::Scalar,3,1> v2; v2 << 0,1,0;
  Eigen::Matrix<zsw::Scalar,3,1> vr;  vr << 1,1,1;
  krj.addConstraint(v0,v1,v2,vr);

  v0<<1,1,0;
  v1<<0,1,1;
  v2<<1,0,1;
  vr<<0,0,0;
  krj.addConstraint(v0,v1,v2,vr);

  v0<<0,1,1;
  v1<<0,0,1;
  v2<<0,1,0;
  vr<<1,1,0;
  krj.addConstraint(v0,v1,v2,vr);

  Eigen::Matrix<zsw::Scalar,3,1> vc0;
  vc0<<0.1,0.1,0.1;
  BOOST_CHECK(!krj.judge(vc0));

  Eigen::Matrix<zsw::Scalar,3,1> vc1;
  vc1<<1,1,1;
  BOOST_CHECK(!krj.judge(vc1));

  Eigen::Matrix<zsw::Scalar,3,1> vc2;
  vc2<<1,0.5,1;
  BOOST_CHECK(!krj.judge(vc2));

  Eigen::Matrix<zsw::Scalar,3,1> vc3;
  vc3<<0.5,0.5,0.5;
  BOOST_CHECK(krj.judge(vc3));

  Eigen::Matrix<zsw::Scalar,3,1> vc4;
  vc4<<1,0.5,0.1;
  BOOST_CHECK(krj.judge(vc4));
}

BOOST_AUTO_TEST_SUITE_END()

#define BOOST_AUTO_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include "../basic_op.h"

using namespace std;

BOOST_AUTO_TEST_SUITE(test_point2tri_dis)

static double dist2_point2line(const Eigen::Matrix<zsw::Scalar,3,1> &p,
                               const Eigen::Matrix<zsw::Scalar,3,1> &d, //d is normalized
                               const Eigen::Matrix<zsw::Scalar,3,1> &q,
                               double &t)
{
  t = d.dot(q - p);
  return (q - (p + t * d)).squaredNorm();
}

static double dist2_point2lineseg(const Eigen::Matrix<zsw::Scalar,3,1> &p0,
                                  const Eigen::Matrix<zsw::Scalar,3,1> &p1,
                                  const Eigen::Matrix<zsw::Scalar,3,1> &q,
                                  Eigen::Matrix<zsw::Scalar,3,1> &proj_pt,
                                  double &t)
{
  const Eigen::Matrix<zsw::Scalar,3,1> &d = (p1 - p0).normalized();

  double dist2 = dist2_point2line(p0, d, q, t);

  if (t < 0){
    t = 0;
    proj_pt = p0;
    dist2 = (q - p0).squaredNorm();
  }
  else if (t > 1){
    t = 1;
    proj_pt = p1;
    dist2 = (q - p1).squaredNorm();
  }
  else{
    proj_pt = p0 + t*d;
  }

  return dist2;
}

static int dist2_point2triangle(const Eigen::Matrix<zsw::Scalar,3,1> &p0,
                                const Eigen::Matrix<zsw::Scalar,3,1> &p1,
                                const Eigen::Matrix<zsw::Scalar,3,1> &p2,
                                const Eigen::Matrix<zsw::Scalar,3,1> &fn,
                                const Eigen::Matrix<zsw::Scalar,3,1> &q,
                                double threshold,
                                Eigen::Matrix<zsw::Scalar,3,1> &proj_pt,
                                double &dist2)
{
  const double dist = -(q - p0).dot(fn);

  if (fabs(dist) > threshold)
    return -1;

  proj_pt = q + dist * fn; //project point

  bool flag = true;
  Eigen::Matrix<zsw::Scalar,3,1> nor = fn.cross(p1 - p0);
  if (nor.dot(p2 - p0) * nor.dot(proj_pt - p0) < 0){
    flag = false;
  }
  else{
    nor = fn.cross(p2 - p1);
    if (nor.dot(p0 - p1) * nor.dot(proj_pt - p1) < 0){
      flag = false;
    }
    else{
      nor = fn.cross(p0 - p2);
      if (nor.dot(p1 - p2) * nor.dot(proj_pt - p2) < 0){
        flag = false;
      }
    }
  }

  if (!flag){//outside the triangle, search the boundary
    double t;
    double tmp_dist2 = dist2_point2lineseg(p0, p1, q, proj_pt, t);
    dist2 = tmp_dist2;
    int    eid = 1;
    Eigen::Matrix<zsw::Scalar,3,1> pt;
    tmp_dist2 = dist2_point2lineseg(p1, p2, q, pt,t);
    if (dist2 > tmp_dist2){
      dist2 = tmp_dist2;
      proj_pt = pt;
      ++eid;
    }
    tmp_dist2 = dist2_point2lineseg(p2, p0, q, pt, t);
    if (dist2 > tmp_dist2){
      dist2 = tmp_dist2;
      proj_pt = pt;
      ++eid;
    }
    return eid;
  }
  dist2 = (q - proj_pt).squaredNorm();
  return 0;
}

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

BOOST_AUTO_TEST_CASE(test_1)
{
  Eigen::Matrix<zsw::Scalar,3,1> v0, v1, v2, vt;
  v0 = Eigen::Matrix<zsw::Scalar,3,1>::Random();
  v1 = Eigen::Matrix<zsw::Scalar,3,1>::Random();
  v2 = Eigen::Matrix<zsw::Scalar,3,1>::Random();
  vt = Eigen::Matrix<zsw::Scalar,3,1>::Random();
  zsw::Scalar sq_dis = zsw::calcPoint2TriSquaredDis(vt,v0,v1,v2);
  zsw::Scalar sq_dis2;
  Eigen::Matrix<zsw::Scalar,3,1> fn, proj_pt;
  fn = (v1-v0).cross(v2-v0); fn.normalize();
  dist2_point2triangle(v0,v1,v2, fn, vt, 1000, proj_pt, sq_dis2);
  if(fabs(sq_dis-sq_dis2)>1e-3) {
    std::cerr << "v0 = [" << v0.transpose() << std::endl;
    std::cerr << "v1 = [" << v1.transpose() << std::endl;
    std::cerr << "v2 = [" << v2.transpose() << std::endl;
    std::cerr << "vt = [" <<  vt.transpose() << std::endl;
  }
  BOOST_CHECK(fabs(sq_dis-sq_dis2)<1e-3);
}

BOOST_AUTO_TEST_SUITE_END()

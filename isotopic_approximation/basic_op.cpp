#include "basic_op.h"

#if 0
zsw::Scalar zsw::calcPoint2TriSquaredDis(const Eigen::Matrix<zsw::Scalar,3,1> &pt,
                                         const Eigen::Matrix<zsw::Scalar,3,1> &v0,
                                         const Eigen::Matrix<zsw::Scalar,3,1> &v1,
                                         const Eigen::Matrix<zsw::Scalar,3,1> &v2)
{
  // calc the projection point
  Eigen::Matrix<zsw::Scalar,3,2> A;
  A.block<3,1>(0,0)=v1-v0; A.block<3,1>(0,1)=v2-v0;
  Eigen::Matrix<zsw::Scalar,2,3> AT = A.transpose();
  Eigen::PartialPivLU<Eigen::Matrix<zsw::Scalar,2,2>> pplu;
  pplu.compute(AT*A);
  Eigen::Matrix<zsw::Scalar,2,1> ans = pplu.solve(AT*(pt-v0));
  Eigen::Matrix<zsw::Scalar,3,1> project_p =A*ans+v0;
  zsw::Scalar project_f_square_dis=(pt-project_p).squaredNorm();
  if(ans[0]>-zsw::const_val::eps && ans[1]>-zsw::const_val::eps && ans[0]+ans[1]<1+zsw::const_val::eps) {
    return project_f_square_dis;
  }

  // three edge condition, judge if is in the inner side of edge
  Eigen::Matrix<zsw::Scalar,3,1> a = (v1-v0);
  Eigen::Matrix<zsw::Scalar,3,1> na = a.normalize();
  Eigen::Matrix<zsw::Scalar,3,1> b = (v2-v0);
  Eigen::Matrix<zsw::Scalar,3,1> nb = b.normalize();
  zsw::Scalar jv;
  if(ans[1]<=0) {
    jv = ans[0]*a.norm() + ans[1]*na.dot(b);
    if(jv<0) {
      return (pt-v0).squaredNorm();
    } else if(jv>1) {
      return (pt-v1).squaredNorm();
    } else {
      return ((project_p-v0).cross(na)).squaredNorm()+project_f_square_dis;
    }
  }
  if(ans[0]<=0) {
    jv = ans[1]*b.norm() + ans[0]*nb.dot(a);
    if(jv<0) {
      return (pt-v0).squaredNorm();
    } else if(jv>1) {
      return (pt-v2).squaredNorm();
    } else {
      return ((project_p-v0).cross(nb)).squaredNorm()+project_f_square_dis;
    }
  }
  zsw::Scalar y[2] = {1-ans[0]-ans[1], ans[0]};
  assert(y[0] < 0); if(y[0] > 0) { std::cerr << "ERROR!!!!" << __FILE__ << std::endl; exit(__LINE__); }
  Eigen::Matrix<zsw::Scalar,3,1> c = a - b;
  Eigen::Matrix<zsw::Scalar,3,1> nc = c.normalize();
  jv = -y[0]*b.dot(nc)+y[1]*c.norm();
  if(jv<0) {
    return (pt-v2).squaredNorm();
  } else if(jv>1) {
    return (pt-v1).squaredNorm();
  } else {
    return ((project_p-v2).cross(nc)).squaredNorm()+project_f_square_dis;
  }
}
#endif


zsw::Scalar zsw::calcPoint2TriSquaredDis(const Eigen::Matrix<zsw::Scalar,3,1> &pt,
                                         const Eigen::Matrix<zsw::Scalar,3,1> &v0,
                                         const Eigen::Matrix<zsw::Scalar,3,1> &v1,
                                         const Eigen::Matrix<zsw::Scalar,3,1> &v2)
{
  // calc the projection point
  Eigen::Matrix<zsw::Scalar,3,2> A;
  A.block<3,1>(0,0)=v1-v0; A.block<3,1>(0,1)=v2-v0;
  Eigen::Matrix<zsw::Scalar,2,3> AT = A.transpose();
  Eigen::PartialPivLU<Eigen::Matrix<zsw::Scalar,2,2>> pplu;
  pplu.compute(AT*A);
  Eigen::Matrix<zsw::Scalar,2,1> ans = pplu.solve(AT*(pt-v0));
  if(ans[0]>0 && ans[1]>0 && ans[0]+ans[1]<1) {
    return (pt-v0 - A*ans).squaredNorm();
  }

  // init ret
  zsw::Scalar ret=min((pt-v0).squaredNorm(), (pt-v1).squaredNorm());
  ret = min(ret, (pt-v2).squaredNorm());

  {    // project to v0 v1
    Eigen::Matrix<zsw::Scalar,3,1> a = v1-v0;
    Eigen::Matrix<zsw::Scalar,3,1> na = a.normalize();
    zsw::Scalar jv=(pt-v0).dot(na);
    if(jv>0 && jv<a.norm()) { // on segement v0 v1
      ret = min(ret, ((pt-v0).cross(na)).squaredNorm());
    }
  }

  {   // project to v0 v2
    Eigen::Matrix<zsw::Scalar,3,1> b = v2-v0;
    Eigen::Matrix<zsw::Scalar,3,1> nb = b.normalize();
    zsw::Scalar jv=(pt-v0).dot(nb);
    if(jv>0 && jv<b.norm()) {
      ret = min(ret, ((pt-v0).cross(nb)).squaredNorm());
    }
  }

  {  // project to v1 v2
    Eigen::Matrix<zsw::Scalar,3,1> c = v2-v1;
    Eigen::Matrix<zsw::Scalar,3,1> nc = c.normalize();
    zsw::Scalar jv=(pt-v1).dot(nc);
    if(jv>0 && jv<c.norm()) {
      ret = min(ret, ((pt-v1).cross(nc)).squaredNorm());
    }
  }
  return ret;
}

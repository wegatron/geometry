#include "sampling.h"

#include <iostream>
#include <Eigen/Geometry>
#include <zswlib/const_val.h>

void zsw::Sampler::sampleSigmaDense(const zsw::mesh::TriMesh &tm, const zsw::Scalar sigma, std::vector<Eigen::Matrix<zsw::Scalar, 3, 1>> &samples)
{
  for(zsw::mesh::TriMesh::ConstFaceIter f_it=tm.faces_begin(); f_it!=tm.faces_end(); ++f_it) {
    Eigen::Matrix<zsw::Scalar, 3, 3> tri_points;
    int i=0;
    for(zsw::mesh::TriMesh::ConstFaceVertexIter fv_it=tm.cfv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
      tri_points.block<3,1>(0,i++) = tm.point(*fv_it);
    }
    sampleTriangle(tri_points, sigma, samples);
  }
  for(zsw::mesh::TriMesh::ConstVertexIter v_it=tm.vertices_begin(); v_it!=tm.vertices_end(); ++v_it) {
    samples.push_back(tm.point(*v_it));
  }
}

void zsw::Sampler::sampleTriangle(Eigen::Matrix<zsw::Scalar, 3, 3> tri_points, const zsw::Scalar sigma, std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &samples)
{
  Eigen::Matrix<zsw::Scalar, 3, 1> translate;
  Eigen::Matrix<zsw::Scalar, 3, 3> rotate;
  calcLocalCoordinate(tri_points, translate, rotate);
  Eigen::Matrix<zsw::Scalar, 3, 3> inv_rotate = rotate.inverse();
  static const Eigen::Matrix<zsw::Scalar, 1, 3> one = Eigen::Matrix<zsw::Scalar,3,1>::Ones();
  tri_points += translate * one;
  tri_points = rotate * tri_points;
  const bool reflection = tri_points(1,2)<0;
  if(reflection) { tri_points(1,2) = -tri_points(1,2);}
  // assert
  assert((tri_points.block<3,1>(0,0)).squaredNorm() < zsw::const_val::eps);
  assert((tri_points.block<1,3>(2,0)).squaredNorm() < zsw::const_val::eps);
  assert(fabs(tri_points(1,1)) < zsw::const_val::eps);
  // sample local triangle
  const zsw::Scalar ky0 = tri_points(0,2)/tri_points(1, 2); // (x_2 -x_0) / (y_2-y_0)
  const zsw::Scalar ky1 = (tri_points(0,2)-tri_points(0,1))/(tri_points(1, 2) - tri_points(1,1));   // (x_2-x_1) / (y_2-y_1)

  std::vector<Eigen::Matrix<zsw::Scalar, 3, 1>> cur_samples;
  const zsw::Scalar step = sigma*1.4142135623730950488;
  zsw::Scalar tmp_y = step;
  zsw::Scalar bound_x[2];
  zsw::Scalar init_x[2];
  init_x[0] = ky0>0 ? 0 : step*ky0;
  init_x[1] = ky1>0 ? tri_points(0,1)+step*ky1 : tri_points(0,1);
  const zsw::Scalar step_x[2] = { step*ky0, step*ky1};
  bound_x[0] = init_x[0]; bound_x[1] = init_x[1];
  for(;tmp_y<tri_points(1,2); tmp_y+=step) {
    for(zsw::Scalar tmp_x=bound_x[0]; tmp_x<bound_x[1]; tmp_x+=step) {
      // resolve point
      Eigen::Matrix<zsw::Scalar,3,1> tmp_sample;
      tmp_sample << tmp_x+sigma*0.7071067811865475244, tmp_y-sigma*0.7071067811865475244, 0;
      if(resolvePoint(tri_points, tmp_sample)) {
        if(reflection) { tmp_sample(1) = -tmp_sample(1); }
        tmp_sample = inv_rotate*tmp_sample-translate;
        samples.push_back(tmp_sample);
      }
    }
    bound_x[0] += step_x[0];
    bound_x[1] += step_x[1];
  }
  // sample the angle area of point 2
  {
    bound_x[0] = std::min(bound_x[0]+step_x[0], tri_points(0,2));
    bound_x[1] = std::max(bound_x[1]+step_x[1], tri_points(0,2));
    for(zsw::Scalar tmp_x=bound_x[0]; tmp_x<bound_x[1]; tmp_x+=step) {
      Eigen::Matrix<zsw::Scalar,3,1> tmp_sample;
      tmp_sample << tmp_x+sigma*0.7071067811865475244, tmp_y-sigma*0.7071067811865475244, 0;
      if(resolvePoint(tri_points, tmp_sample)) {
        if(reflection) { tmp_sample(1) = -tmp_sample(1); }
        tmp_sample = inv_rotate*tmp_sample-translate;
        samples.push_back(tmp_sample);
      }
    }
  }
}

void zsw::Sampler::calcLocalCoordinate(const Eigen::Matrix<zsw::Scalar, 3, 3> &tri_points, Eigen::Matrix<zsw::Scalar, 3, 1> &translate, Eigen::Matrix<zsw::Scalar, 3, 3> &rotate)
{
  translate = -tri_points.block<3,1>(0,0);
  Eigen::Matrix<zsw::Scalar,3,1> m;
  m<<0,0,1;
  Eigen::Matrix<zsw::Scalar,3,1> n = (tri_points.block<3,1>(0,1) - tri_points.block<3,1>(0,0)).cross(tri_points.block<3,1>(0,2)-tri_points.block<3,1>(0,0));
  // rotate to xy plane
  Eigen::Quaternion<zsw::Scalar> qn = Eigen::Quaternion<zsw::Scalar>::FromTwoVectors(n, m);
  rotate = qn.toRotationMatrix();

  // let x1 one the x axis
  m<<1,0,0;
  n=rotate*(tri_points.block<3,1>(0,1)-tri_points.block<3,1>(0,0));
  qn = Eigen::Quaternion<zsw::Scalar>::FromTwoVectors(n,m);
  rotate = qn.toRotationMatrix() * rotate;
}

bool zsw::Sampler::resolvePoint(const Eigen::Matrix<zsw::Scalar,3,3> &tri_points, Eigen::Matrix<zsw::Scalar,3,1> &sample_point)
{
  // adjust the sample point into the triangle
  if(!sameSide2D(tri_points.block<3,1>(0,0), tri_points.block<3,1>(0,1), tri_points.block<3,1>(0,2), sample_point)) {
    projectToLine(tri_points.block<3,1>(0,0), tri_points.block<3,1>(0,1), sample_point);
    if((sample_point-tri_points.block<3,1>(0,0)).dot(tri_points.block<3,1>(0,1)-sample_point) < -zsw::const_val::eps) { // sample point is not on v0-v1
      return false;
      // if((sample_point-tri_points.block<3,1>(0,0)).dot(tri_points.block<3,1>(0,1)-tri_points.block<3,1>(0,0)) > 0) {
      //   sample_point = tri_points.block<3,1>(0,1);
      // } else { sample_point = tri_points.block<3,1>(0,0); }
    }
  } else if(!sameSide2D(tri_points.block<3,1>(0,0), tri_points.block<3,1>(0,2), tri_points.block<3,1>(0,1), sample_point)) {
    projectToLine(tri_points.block<3,1>(0,0), tri_points.block<3,1>(0,2), sample_point);
    if((sample_point-tri_points.block<3,1>(0,0)).dot(tri_points.block<3,1>(0,2)-sample_point) < -zsw::const_val::eps) { // sample point is not on v0-v2
      return false;
      // if((sample_point-tri_points.block<3,1>(0,0)).dot(tri_points.block<3,1>(0,2)-tri_points.block<3,1>(0,0)) > 0) {
      //   sample_point = tri_points.block<3,1>(0,2);
      // } else { sample_point = tri_points.block<3,1>(0,0); }
    }
  } else if(!sameSide2D(tri_points.block<3,1>(0,1), tri_points.block<3,1>(0,2), tri_points.block<3,1>(0,0), sample_point)) {
    projectToLine(tri_points.block<3,1>(0,1), tri_points.block<3,1>(0,2), sample_point);
    if((sample_point-tri_points.block<3,1>(0,1)).dot(tri_points.block<3,1>(0,2)-sample_point) < -zsw::const_val::eps) { // sample point is not on v1-v2
      return false;
      // if((sample_point-tri_points.block<3,1>(0,1)).dot(tri_points.block<3,1>(0,2)-tri_points.block<3,1>(0,1)) > 0) {
      //   sample_point = tri_points.block<3,1>(0,2);
      // } else { sample_point = tri_points.block<3,1>(0,1); }
    }
  }
  return true;
}

bool zsw::sameSide2D(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> &v1, const Eigen::Matrix<zsw::Scalar,3,1> &vr, const Eigen::Matrix<zsw::Scalar,3,1> &vt)
{
  Eigen::Matrix<zsw::Scalar,3,3> v_check;
  v_check.block<3,1>(0,0) = v1-v0;
  v_check.block<3,1>(0,1) = vr-v0;
  v_check.block<3,1>(0,2) = vt-v0;
  assert(fabs(v_check.determinant()) < zsw::const_val::eps);
  Eigen::Matrix<zsw::Scalar,3,1> vn0 = v_check.block<3,1>(0,0).cross(v_check.block<3,1>(0,1));
  Eigen::Matrix<zsw::Scalar,3,1> vn1 = v_check.block<3,1>(0,0).cross(v_check.block<3,1>(0,2));
  return vn0.dot(vn1)>zsw::const_val::eps;
}

void zsw::Sampler::projectToLine(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> &v1, Eigen::Matrix<zsw::Scalar,3,1> &sample_point)
{
  Eigen::Matrix<zsw::Scalar,3,1> ab = v1-v0;
  Eigen::Matrix<zsw::Scalar,3,1> ac = sample_point-v0;
  Eigen::Matrix<zsw::Scalar,3,1> cb = v1-sample_point;
  Eigen::Matrix<zsw::Scalar,3,1> n=(ab.cross(ac)).cross(ab);  n.normalize();
  sample_point += (cb.dot(n)) * n;
}

#if 0
void zsw::sampleTet(const zsw::Scalar dense, const Eigen::Matrix<zsw::Scalar, 3, 4> &tet_points, std::vector<Point> &sample_points)
{
  Eigen::Matrix<zsw::Scalar, 3, 2> bbox;
  calcBbox<zsw::Scalar, 3, 4>(tet_points, bbox);
  for(zsw::Scalar v0=bbox(0,0); v0<bbox(0,1); v0+=dense) {
    for(zsw::Scalar v1=bbox(1,0); v1<bbox(1,1); v1+=dense) {
      for(zsw::Scalar v2=bbox(2,0); v2<bbox(2,1); v2+=dense) {
        Point pt(v0,v1,v2);
        if(inTet(pt, tet_points)) {
          sample_points.push_back(pt);
        }
      }
    }
  }
}
#endif

bool zsw::inTet(const Point &pt, const Eigen::Matrix<zsw::Scalar,3,4> &tet_points)
{
  Eigen::Matrix<zsw::Scalar,3,1> vt;
  vt << pt[0], pt[1], pt[2];
  for(size_t i=0; i<4; ++i) {
    Eigen::Matrix<zsw::Scalar,3,1> v[3];
    size_t id=0;
    for(size_t pi=0; pi<4; ++pi) {
      if(pi==i) { continue; }
      v[id++] = tet_points.block<3,1>(0,pi);
    }
    if(!sameSide3D(v[0], v[1], v[2], tet_points.block<3,1>(0,i), vt)) {
      return false;
    }
  }
  return true;
}

bool zsw::sameSide3D(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> v1, const Eigen::Matrix<zsw::Scalar,3,1> &v2, const Eigen::Matrix<zsw::Scalar,3,1> &vr, const Eigen::Matrix<zsw::Scalar,3,1> &vt)
{
  Eigen::Matrix<zsw::Scalar,3,1> n=(v1-v0).cross(v2-v0);
  zsw::Scalar value[2]={n.dot(vr-v0), n.dot(vt-v0)};
  return (value[0]>=0 && value[1]>=0) || (value[0]<=0 && value[1]<=0);
}

void zsw::sampleTriangle(const Eigen::Matrix<zsw::Scalar, 3, 3> &tri_points, const zsw::Scalar r,
                         const zsw::Scalar val_exp, std::vector<zsw::JudgePoint> &samples)
{
  Eigen::Matrix<zsw::Scalar,3,1> v0=tri_points.block<3,1>(0,0);
  Eigen::Matrix<zsw::Scalar,3,1> v1=tri_points.block<3,1>(0,1);
  Eigen::Matrix<zsw::Scalar,3,1> v2=tri_points.block<3,1>(0,2);
  Eigen::Matrix<zsw::Scalar,3,1> n=v1-v0;
  Eigen::Matrix<zsw::Scalar,3,1> m=v2-v0;
  zsw::Scalar step_n=r/n.norm();
  zsw::Scalar step_m=r/m.norm();
  bool sn_flg=true;
  for(zsw::Scalar sn=0; sn_flg; sn+=step_n) {
    if(sn>1) { sn=1; sn_flg=false; }
    const Eigen::Matrix<zsw::Scalar,3,1> tmp_n=sn*n;
    const zsw::Scalar max_sm=1-sn-zsw::const_val::eps;
    bool sm_flg=true;
    for(zsw::Scalar sm=0; sm_flg; sm+=step_m) {
      if(sm>max_sm) { sm=max_sm; sm_flg=false; }
      samples.push_back({v0+tmp_n+sm*m, val_exp, val_exp});
    }
  }
}

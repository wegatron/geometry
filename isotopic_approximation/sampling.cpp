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
      resolvePoint(tri_points, tmp_sample);
      if(reflection) { tmp_sample(1) = -tmp_sample(1); }
      tmp_sample = inv_rotate*tmp_sample-translate;
      samples.push_back(tmp_sample);
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
      resolvePoint(tri_points, tmp_sample);
      if(reflection) { tmp_sample(1) = -tmp_sample(1); }
      tmp_sample = inv_rotate*tmp_sample-translate;
      samples.push_back(tmp_sample);
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
  n=tri_points.block<3,1>(0,1)-tri_points.block<3,1>(0,0);
  qn = Eigen::Quaternion<zsw::Scalar>::FromTwoVectors(n,m);
  rotate = qn.toRotationMatrix() * rotate;
}

void zsw::Sampler::resolvePoint(const Eigen::Matrix<zsw::Scalar,3,3> &tri_points, Eigen::Matrix<zsw::Scalar,3,1> &sample_point)
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;

}

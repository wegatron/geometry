#include <iostream>
#include <fstream>
#include "../sampling.h"

using namespace std;

void test_calcLocalCoordinate()
{
  zsw::Sampler sampler;
  Eigen::Matrix<zsw::Scalar,3,3> tri_points;
  tri_points<< 1,0,0,
    0,1,0,
    0,0,1;
  Eigen::Matrix<zsw::Scalar,3,1> translate;
  Eigen::Matrix<zsw::Scalar,3,3> rotate;
  sampler.calcLocalCoordinate(tri_points, translate, rotate);

  Eigen::Matrix<zsw::Scalar,3,3> local_points = tri_points;
  local_points += translate * Eigen::Matrix<zsw::Scalar,1,3>::Ones();
  local_points = rotate * local_points;

  std::cerr << "local_points:\n" << local_points << std::endl;

  std::cerr << "points transform back:\n" << (rotate.inverse() * local_points) - translate * Eigen::Matrix<zsw::Scalar,1,3>::Ones() << std::endl;
}

void test_sampleTriangle()
{
  zsw::Sampler sampler;
  Eigen::Matrix<zsw::Scalar,3,3> tri_points;
  tri_points<< 1,0,0,
    0,1,0,
    0,0,1;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> samples;
  sampler.sampleTriangle(tri_points, 0.1, samples);
  // write obj file
  ofstream ofs("/home/wegatron/tmp.obj", std::ofstream::out);
  ofs << "f 1 2 3" << std::endl;
  for(int i=0; i<3; ++i) {
    ofs << "v " << tri_points.block<3,1>(0,i).transpose() << std::endl;
  }
  for(Eigen::Matrix<zsw::Scalar,3,1> &tmp_sample : samples) {
    ofs << "v " << tmp_sample.transpose() << std::endl;
  }
  ofs.close();
}

int main(int argc, char *argv[])
{
  test_sampleTriangle();
  //test_calcLocalCoordinate();
  return 0;
}

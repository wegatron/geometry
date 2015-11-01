#include <fstream>

#include <zswlib/mesh/vtk.h>
#include <zswlib/error_ctrl.h>
#include "../sampling.h"

using namespace std;

void test_sampleTriangle()
{
  Eigen::Matrix<zsw::Scalar,3,3> tri_points;
  tri_points<< 1,0,0,
    0,1,0,
    0,0,1;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> samples;
  zsw::sampleTriangle(tri_points, 0.1, samples);
  // write obj file
  std::ofstream ofs("/home/wegatron/tmp.obj", std::ofstream::out);
  ofs << "f 1 2 3" << std::endl;
  for(int i=0; i<3; ++i) {
    ofs << "v " << tri_points.block<3,1>(0,i).transpose() << std::endl;
  }
  for(Eigen::Matrix<zsw::Scalar,3,1> &tmp_sample : samples) {
    ofs << "v " << tmp_sample.transpose() << std::endl;
  }
  ofs.close();
}

void test_sample_tet()
{
  Eigen::Matrix<zsw::Scalar,3,4> tet_points;
  tet_points << 2,0,0,1,
                          0,4,-1,1,
                          0,0,0,4;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> sample_points;
  zsw::sampleTet(tet_points, 0.4, sample_points);
  std::ofstream ofs;
  OPEN_STREAM("/home/wegatron/test_sample_tet.vtk", ofs, std::ofstream::out, return);
  vector<zsw::Scalar> pt_data(12, 0);
  copy(tet_points.data(), tet_points.data()+12, &pt_data[0]);
  for(const Eigen::Matrix<zsw::Scalar,3,1> &pt : sample_points) {
    pt_data.push_back(pt[0]);
    pt_data.push_back(pt[1]);
    pt_data.push_back(pt[2]);
  }
  size_t tets[4] = {0,1,2,3};
  tet2vtk(ofs, &pt_data[0], 4+sample_points.size(), &tets[0], 1);
}

int main(int argc, char *argv[])
{
  //test_sampleTriangle();
  test_sample_tet();
  return 0;
}

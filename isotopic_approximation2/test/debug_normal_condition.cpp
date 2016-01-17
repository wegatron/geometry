#include <iostream>
#include <fstream>

#include "../constraint.h"

using namespace std;

void readJptsFromVtk(const std::string &file_path, std::vector<Eigen::Matrix<zsw::Scalar,3,1>>  &jpts)
{
  ifstream ifs(file_path);
  size_t n;
  ifs>>n;
  jpts.reserve(n);
  Eigen::Matrix<zsw::Scalar,3,1> pt;
  for(size_t i=0; i<n; ++i) {
    ifs >> pt[0] >> pt[1] >> pt[2];
    jpts.push_back(pt);
  }
}

int main(int argc, char *argv[])
{
  // inner jpts
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> inner_jpts; // input
  readJptsFromVtk("/home/wegatron/tmp/normal_cond_debug/inner_jpts.vtk", inner_jpts);
  std::shared_ptr<zsw::Flann<zsw::Scalar>> inner_kdtree_ptr(new zsw::Flann<zsw::Scalar>(inner_jpts[0].data(), inner_jpts.size()));
  // outer jpts
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> outer_jpts; // input
  readJptsFromVtk("/home/wegatron/tmp/normal_cond_debug/outer_jpts.vtk", outer_jpts);
  std::shared_ptr<zsw::Flann<zsw::Scalar>> outer_kdtree_ptr(new zsw::Flann<zsw::Scalar>(outer_jpts[0].data(), outer_jpts.size()));
  // tris_pts
  Eigen::Matrix<zsw::Scalar,3,4> tri_pts; // input
  tri_pts<<63.1166,76.0815,73.8558,65.7652,
    -14.6293,-2.0456,-15.1596,-17.2288,
    125.97,125.865,114.733,119.895;
  Eigen::Matrix<zsw::Scalar,3,4> scaled_tri_pts;
  Eigen::Matrix<zsw::Scalar,3,1> bc=0.25*(tri_pts.block<3,1>(0,0)+tri_pts.block<3,1>(0,1)+tri_pts.block<3,1>(0,2)
                                          +tri_pts.block<3,1>(0,3));
  Eigen::Matrix<zsw::Scalar,1,4> tmp_one=Eigen::Matrix<zsw::Scalar,1,4>::Ones();
  scaled_tri_pts=0.7*tri_pts+0.3*bc*tmp_one;
  Eigen::Matrix<zsw::Scalar,4,1> vals; // input
  vals<<1,1,1,-1;
  std::string file_path="/home/wegatron/tmp/normal_cond_debug/denug.vtk";
  if(!zsw::normalCondition(vals, scaled_tri_pts, tri_pts,inner_jpts, outer_jpts, inner_kdtree_ptr, outer_kdtree_ptr,
                           true, &file_path)) {
    std::cout << "normal condition failed!" << std::endl;
  }
  return 0;
}

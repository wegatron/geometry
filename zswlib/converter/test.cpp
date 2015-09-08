#include <iostream>
#include <vector>

#include <boost/foreach.hpp>

#include "matrix_conv.h"

using namespace std;

template<typename VAL_TYPE>
void test_func(std::vector<VAL_TYPE> vec)
{
  for (int i=0; i<vec.size(); ++i)
    cout << vec[i] << " ";
}

int main(int argc, char *argv[])
{
  // std::vector<Eigen::Vector3d> nodes;
  // for (int i=0; i<10; ++i) {
  //   Eigen::Vector3d tmp_v(i, i+1, i+2);
  //   nodes.push_back(tmp_v);
  // }
  // zjucad::matrix::matrix<double> zju_mat;
  // zsw::vvec3ToZjumat(nodes, zju_mat);

  // cout << "zju_mat:" << endl;
  // for (int i=0; i<zju_mat.size(1); ++i) {
  //   for (int j=0; j<zju_mat.size(2); ++j) {
  //     cout << zju_mat(i,j) << " ";
  //   }
  //   cout << endl;
  // }

  // cout << "convert back:" << endl;
  // nodes.clear();
  // zsw::zjumat3ToVvec(zju_mat, nodes);
  // BOOST_FOREACH(auto tmp_v, nodes) {
  //   cout << tmp_v.transpose() << " ";
  // }

  std::vector<Eigen::Vector4i> tets;
  for (int i=0; i< 10; ++i) {
    Eigen::Vector4i tmp_v(i, i+1, i+2, i+3);
    tets.push_back(tmp_v);
  }

  zjucad::matrix::matrix<int> zju_mat;
  zsw::vvec4ToZjumat(tets, zju_mat);

  cout << "zju_mat:" << endl;
  for (int i=0; i<zju_mat.size(1); ++i) {
    for (int j=0; j<zju_mat.size(2); ++j) {
      cout << zju_mat(i,j) << " ";
    }
    cout << endl;
  }

  cout << "convert back:" << endl;
  tets.clear();
  zsw::zjumat4ToVvec(zju_mat, tets);
  BOOST_FOREACH(auto tmp_v, tets) {
    cout << tmp_v.transpose() << ";";
  }
  return 0;
}

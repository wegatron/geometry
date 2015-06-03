#include<iostream>

#include "../bilateral_normal_filter.h"

using namespace std;
using namespace zsw;

void filterTrimesh(const string &obj_prefix)
{
  #pragma omp  parallel for
  for (size_t i=1; i<10; ++i) {
    string file_path = obj_prefix+".obj";
    jtf::mesh::tri_mesh trimesh(file_path.c_str());
    BilateralNormalFilter bnf;
    bnf.setSt(i);
    bnf.filter(trimesh);
    writeVtk(obj_prefix+to_string(i)+".vtk", trimesh.trimesh_.mesh_, trimesh.trimesh_.node_);
  }
}

int main(int argc, char *argv[])
{
  // filterTrimesh("/home/wegatron/workspace/geometry/result/bilateral_filter/zsw_test_NII/cube/cube");
  filterTrimesh("/home/wegatron/workspace/geometry/result/bilateral_filter/zsw_test_NII/2/Tooth_15");
  std::cout << __FILE__ << __LINE__ << std::endl;
  filterTrimesh("/home/wegatron/workspace/geometry/result/bilateral_filter/zsw_test_NII/3/UpperJaw_after");
  std::cout << __FILE__ << __LINE__ << std::endl;
  filterTrimesh("/home/wegatron/workspace/geometry/result/bilateral_filter/zsw_test_NII/4/buda");
  std::cout << __FILE__ << __LINE__ << std::endl;
  filterTrimesh("/home/wegatron/workspace/geometry/result/bilateral_filter/zsw_test_NII/6/UpperJaw");
// #pragma omp parallel sections
//   {
//     {
//       filterTrimesh("/home/wegatron/workspace/geometry/result/bilateral_filter/zsw_test_NII/cube/cube");
//     }
// #pragma omp section
//     {
//       filterTrimesh("/home/wegatron/workspace/geometry/result/bilateral_filter/zsw_test_NII/2/Tooth_15");
//     }
// #pragma omp section
//     {
//       filterTrimesh("/home/wegatron/workspace/geometry/result/bilateral_filter/zsw_test_NII/3/UpperJaw_after");
//     }
// #pragma omp section
//     {
//     filterTrimesh("/home/wegatron/workspace/geometry/result/bilateral_filter/zsw_test_NII/4/buda");
//   }
// #pragma omp section
//     {
//       filterTrimesh("/home/wegatron/workspace/geometry/result/bilateral_filter/zsw_test_NII/6/UpperJaw");
//     }
//   }
  return 0;
}

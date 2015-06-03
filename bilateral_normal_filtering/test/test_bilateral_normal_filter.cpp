#include <iostream>
#include <chrono>

#include "../bilateral_normal_filter.h"

using namespace std;
using namespace zsw;

void filterTrimesh(const string &obj_prefix)
{
  // #pragma omp  parallel for
  // for (size_t i=1; i<10; ++i) {
  //   string file_path = obj_prefix+".obj";
  //   jtf::mesh::tri_mesh trimesh(file_path.c_str());
  //   BilateralNormalFilter bnf;
  //   bnf.setSt(i);
  //   bnf.filter(trimesh);
  //   writeVtk(obj_prefix+to_string(i)+".vtk", trimesh.trimesh_.mesh_, trimesh.trimesh_.node_);
  // }
  auto t1 = std::chrono::high_resolution_clock::now();

  string file_path = obj_prefix+".obj";
  jtf::mesh::tri_mesh trimesh(file_path.c_str());
  BilateralNormalFilter bnf;
  bnf.filter(trimesh);

  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "load and filter took "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
            << " milliseconds\n" << std::endl;

  writeVtk(obj_prefix+to_string(3)+".vtk", trimesh.trimesh_.mesh_, trimesh.trimesh_.node_);
}

int main(int argc, char *argv[])
{
  // // filterTrimesh("/home/wegatron/workspace/geometry/result/bilateral_filter/zsw_test_NII/cube/cube");
  // filterTrimesh("/home/wegatron/workspace/geometry/result/bilateral_filter/zsw_test_NII/2/Tooth_15");
  // std::cout << __FILE__ << __LINE__ << std::endl;
  // filterTrimesh("/home/wegatron/workspace/geometry/result/bilateral_filter/zsw_test_NII/3/UpperJaw_after");
  // std::cout << __FILE__ << __LINE__ << std::endl;
  // filterTrimesh("/home/wegatron/workspace/geometry/result/bilateral_filter/zsw_test_NII/4/buda");
  // std::cout << __FILE__ << __LINE__ << std::endl;
  // filterTrimesh("/home/wegatron/workspace/geometry/result/bilateral_filter/zsw_test_NII/6/UpperJaw");
  if(argc != 2) {
    std::cerr << "usage: bnf [obj_file_prefix]" << std::endl;
    return __LINE__;
  }
  filterTrimesh(argv[1]);
  return 0;
}

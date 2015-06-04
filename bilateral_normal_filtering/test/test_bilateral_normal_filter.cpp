#include <iostream>
#include <chrono>

#include "../bilateral_normal_filter.h"

using namespace std;
using namespace zsw;

void filterTrimesh(const string &obj_prefix)
{
  auto t1 = std::chrono::high_resolution_clock::now();

  string file_path = obj_prefix+".obj";
  jtf::mesh::tri_mesh trimesh(file_path.c_str());
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "load and mesh and edge2cell "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
            << " milliseconds\n" << std::endl;

  BilateralNormalFilter bnf;
  bnf.filter(trimesh);
  auto t3 = std::chrono::high_resolution_clock::now();
  std::cout << "filter took "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count()
            << " milliseconds\n" << std::endl;

  writeVtk(obj_prefix+to_string(3)+".vtk", trimesh.trimesh_.mesh_, trimesh.trimesh_.node_);
}

int main(int argc, char *argv[])
{
  if(argc != 2) {
    std::cerr << "usage: bnf [obj_file_prefix]" << std::endl;
    return __LINE__;
  }
  filterTrimesh(argv[1]);
  return 0;
}

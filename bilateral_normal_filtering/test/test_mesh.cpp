#include <iostream>

#include "jtflib/mesh/trimesh.h"

using namespace std;

int main(int argc, char *argv[])
{
  jtf::mesh::tri_mesh trimesh("/home/wegatron/tmp/Tooth_15.obj");
  std::pair<size_t, size_t> result = trimesh.ea_->query(1049, 1148);
  std::cout << result.first << " " << result.second << std::endl;
  return 0;
}

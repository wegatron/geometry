#include<iostream>

#include "../bilateral_normal_filter.h"

using namespace std;
using namespace zsw;

void filterTrimesh(const string &obj_prefix)
{
  for (size_t i=1; i<10; ++i) {
    // jtf::mesh::tri_mesh trimesh("/home/wegatron/tmp/Tooth_15.obj");
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
  // filterTrimesh("/home/wegatron/workspace/geometry/result/bilateral_filter/cube/cube");
#pragma omp parallel sections
  {
    {
      filterTrimesh("/home/wegatron/workspace/geometry/result/bilateral_filter/newTest/3/UpperJaw_after");
    }
#pragma omp section
    {
      filterTrimesh("/home/wegatron/workspace/geometry/result/bilateral_filter/newTest/4/buda");
    }
#pragma omp section
    {
    filterTrimesh("/home/wegatron/workspace/geometry/result/bilateral_filter/newTest/2/Tooth_15");
  }
#pragma omp section
    {
      filterTrimesh("/home/wegatron/workspace/geometry/result/bilateral_filter/newTest/6/UpperJaw");
    }
  }
  return 0;
}

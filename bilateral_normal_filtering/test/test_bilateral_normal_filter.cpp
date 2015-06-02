#include<iostream>

#include "../bilateral_normal_filter.h"

using namespace std;
using namespace zsw;

int main(int argc, char *argv[])
{
  jtf::mesh::tri_mesh trimesh("/home/wegatron/tmp/cube_input.obj");
  BilateralNormalFilter bnf;
  bnf.filter(trimesh);
  jtf::mesh::save_obj("/home/wegatron/tmp/cube_output.obj", trimesh.trimesh_.mesh_, trimesh.trimesh_.node_);
  return 0;
}

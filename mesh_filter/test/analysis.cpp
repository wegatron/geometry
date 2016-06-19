#include "../res_analysis.h"

#include <iostream>
#include <fstream>

#include <boost/lexical_cast.hpp>

void dv(const std::string &file_dir, const size_t st)
{
  jtf::mesh2::meshes mesh_a;
  if(jtf::mesh2::load_obj(std::string(file_dir+"out_put"+boost::lexical_cast<std::string>(st)+".obj").c_str(), mesh_a.mesh_, mesh_a.node_)) {
    std::cerr << "can't open file " + file_dir+ "out_put.obj for read" << std::endl;
    return;
  }

  jtf::mesh2::meshes mesh_b;
  if(jtf::mesh2::load_obj(std::string(file_dir+"smooth_1.obj").c_str(), mesh_b.mesh_, mesh_b.node_)) {
    std::cerr << "can't open file " + file_dir + "smooth_1.obj for read" << std::endl;
    return;
  }

  zsw::mesh::diffMeshVertex(mesh_a, mesh_b, file_dir+"diff"+boost::lexical_cast<std::string>(st)+".vtk");
}

void loop(const std::string &file_dir)
{
  for(size_t i=3; i<8; ++i) {
    std::cout << "analysis one frame " << i << std::endl;
    dv(file_dir, i);
  }
}

int main(int argc, char *argv[])
{
  if(argc != 2) {
    std::cout << "usage: analysis [file_dir]" << std::endl;
    return __LINE__;
  }
  loop(argv[1]);
  return 0;
}

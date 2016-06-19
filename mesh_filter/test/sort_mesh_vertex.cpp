#include "../res_analysis.h"

#include <iostream>
#include <boost/lexical_cast.hpp>
#include <jtflib/mesh2/mesh.h>


int main(int argc, char *argv[])
{
  if(argc != 2) {
    std::cerr << "usage: [sort_mesh_v] [prefix]" << std::endl;
  }
  {
    jtf::mesh2::meshes mesh_a;
    std::string prefix(argv[1]);
    if(jtf::mesh2::load_obj(std::string(prefix+"input.obj").c_str(), mesh_a.mesh_, mesh_a.node_)) {
      std::cerr << "can't open file " + prefix+ "input.obj for read" << std::endl;
      return __LINE__;
    }
    zsw::mesh::sortMeshVertex(mesh_a.mesh_, mesh_a.node_);
    jtf::mesh2::save_obj(std::string(prefix+"sinput.obj").c_str(), mesh_a.mesh_, mesh_a.node_);
  }

  {
    jtf::mesh2::meshes mesh_a;
    std::string prefix(argv[1]);
    if(jtf::mesh2::load_obj(std::string(prefix+"expected.obj").c_str(), mesh_a.mesh_, mesh_a.node_)) {
      std::cerr << "can't open file " + prefix+ "expected.obj for read" << std::endl;
      return __LINE__;
    }
    zsw::mesh::sortMeshVertex(mesh_a.mesh_, mesh_a.node_);
    jtf::mesh2::save_obj(std::string(prefix+"sexpected.obj").c_str(), mesh_a.mesh_, mesh_a.node_);
  }
  return 0;
}

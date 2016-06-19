#include "../res_analysis.h"

#include <iostream>
#include <fstream>

#include "../vtk.h"

void testDv(int argc, char *argv[])
{
  if(argc != 4) {
    std::cout << "usage: [diffv] mesh_a mesh_b out_put" << std::endl;
    return;
  }

  jtf::mesh2::meshes mesh_a;
  if(jtf::mesh2::load_obj(argv[1], mesh_a.mesh_, mesh_a.node_)) {
    std::cerr << "can't open file " + std::string(argv[1]) + " for read" << std::endl;
    return;
  }

  jtf::mesh2::meshes mesh_b;
  if(jtf::mesh2::load_obj(argv[2], mesh_b.mesh_, mesh_b.node_)) {
    std::cerr << "can't open file " + std::string(argv[2]) + " for read" << std::endl;
    return;
  }
  zsw::mesh::diffMeshVertex(mesh_a, mesh_b, argv[3]);
}

void testSortMesh()
{
  jtf::mesh2::meshes mesh_a;
  if(jtf::mesh2::load_obj("E:/workspace/geometry/bilateral_normal_filtering/result/2/Tooth_15.obj", mesh_a.mesh_, mesh_a.node_)) {
    std::cerr << "can't open file E:/workspace/geometry/bilateral_normal_filtering/result/2/Tooth_15.obj for read" << std::endl;
    return;
  }
  zsw::mesh::sortMeshVertex(mesh_a.mesh_, mesh_a.node_);
  jtf::mesh2::save_obj("E:/workspace/geometry/bilateral_normal_filtering/result/2/Tooth_15_sorted.obj",mesh_a.mesh_, mesh_a.node_);
}

void testData2Vtk()
{
  jtf::mesh2::meshes mesh_a;
  if(jtf::mesh2::load_obj("E:/workspace/geometry/bilateral_normal_filtering/result/cube/cube.obj", mesh_a.mesh_, mesh_a.node_)) {
    std::cerr << "can't open file E:/workspace/geometry/bilateral_normal_filtering/result/2/cube/cube.obj for read" << std::endl;
    return;
  }
  std::ofstream ofs("E:/workspace/geometry/bilateral_normal_filtering/result/cube/cube_debug.vtk");
  if(!ofs) {
    std::cerr << "can't open file:E:/workspace/geometry/bilateral_normal_filtering/result/2/cube/cube_debug.obj for write!" << std::endl;
  }

  tri2vtk(ofs, mesh_a.node_.data(), mesh_a.node_.cols(),
          mesh_a.mesh_.data(), mesh_a.mesh_.cols());
  std::vector<double> val(mesh_a.node_.cols(), 0.0);
  point_data(ofs, val.begin(), val.size(), "diff_val");
}

int main(int argc, char *argv[])
{
  // testSortMesh();
  testDv(argc, argv);
  //testData2Vtk();
  return 0;
}

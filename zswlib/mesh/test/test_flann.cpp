#include <iostream>
#include <jtflib/mesh2/mesh.h>
#include <zswlib/mesh/zsw_flann.h>

int main(int argc, char *argv[])
{
  jtf::mesh2::meshes tmp_mesh;
  jtf::mesh2::load_obj("E:/workspace/geometry/bilateral_normal_filtering/result/2/expected.obj", tmp_mesh.mesh_, tmp_mesh.node_);
  zsw::Flann<double> flann(tmp_mesh.node_.data(), (size_t)tmp_mesh.node_.cols());
  Eigen::Matrix<double, 3, Eigen::Dynamic> query = tmp_mesh.node_;
  std::vector<size_t> indices;
  std::vector<double> dist;
  flann.queryNearest(query, indices, dist);
  for(size_t i=0; i<indices.size(); ++i) {
    std::cout << indices[i] << "/" << dist[i] << " ";
  }
  std::cout << std::endl;
  return 0;
}

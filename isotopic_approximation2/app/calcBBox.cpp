#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>

using namespace std;

int main(int argc, char *argv[])
{
  std::string file_path(argv[1]);
  zsw::mesh::TriMesh input_mesh;
  if(!OpenMesh::IO::read_mesh(input_mesh, file_path)) {
    std::cerr << "can't open file " << file_path << std::endl;
    abort();
  }

  Eigen::Matrix<zsw::Scalar,3,1> min_pos=input_mesh.point(*input_mesh.vertices_begin());
  Eigen::Matrix<zsw::Scalar,3,1> max_pos=min_pos;

  for(auto vit=input_mesh.vertices_begin(); vit!=input_mesh.vertices_end(); ++vit) {
    auto tmp_pos=input_mesh.point(*vit);
    for(size_t i=0; i<3; ++i) {
      if(tmp_pos[i]<min_pos[i]) { min_pos[i]=tmp_pos[i]; }
      else if(tmp_pos[i]>max_pos[i]) {max_pos[i]=tmp_pos[i]; }
    }
  }
  std::cout << "bbox\n min:" << min_pos.transpose() << "\nmax:" << max_pos.transpose() << std::endl;
  zsw::Scalar edge_length[3]={max_pos[0]-min_pos[0], max_pos[1]-min_pos[1], max_pos[2]-min_pos[2]};
  std::cout << "edge length:\n" << edge_length[0] << "\n" << edge_length[1] << "\n" << edge_length[2] << std::endl;
  return 0;
}

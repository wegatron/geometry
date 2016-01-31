#include <fstream>
#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/mesh/mesh_type.h>

using namespace std;


void genNoise(const zsw::Scalar min_val, const zsw::Scalar max_val, zsw::mesh::TriMesh &in_mesh)
{
  if(!in_mesh.has_vertex_normals()) {
    in_mesh.request_face_normals();
    in_mesh.request_vertex_normals();
    in_mesh.update_normals();
  }

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(min_val, max_val);

  for(zsw::mesh::TriMesh::VertexIter vit=in_mesh.vertices_begin(); vit!=in_mesh.vertices_end(); ++vit) {
    double dis_val = distribution(generator);
    const Eigen::Matrix<zsw::Scalar, 3, 1> &ep = in_mesh.point(*vit);
    // if(ep[1]>-24.7659) { continue; }
    Eigen::Matrix<zsw::Scalar, 3, 1> offset = in_mesh.normal(*vit) * dis_val;
    in_mesh.set_point(*vit, ep+offset);
  }
}

int main(int argc, char *argv[])
{
  if(argc !=5) {
    std::cout << "usage: gen_noise [in_filepath] [out_filepath] [min_val] [max_val]" << std::endl;
    return 0;
  }

  const std::string in_filepath(argv[1]);
  const std::string out_filepath(argv[2]);
  const double min_val = atof(argv[3]);
  const double max_val = atof(argv[4]);
  std::cout << "gen_noise " << in_filepath << "->"
            << out_filepath << " " << min_val <<  " " << max_val << std::endl;
  zsw::mesh::TriMesh in_mesh;
  if(!OpenMesh::IO::read_mesh(in_mesh, in_filepath)) {
    std::cerr << "[ERROR] can't read mesh file:" << in_filepath << std::endl;
    abort();
  }
  genNoise(min_val, max_val, in_mesh);
  if(!OpenMesh::IO::write_mesh(in_mesh, out_filepath)) {
    std::cerr << "[ERROR] can't write mesh to file:" << out_filepath << std::endl;
    abort();
  }
  return 0;
}

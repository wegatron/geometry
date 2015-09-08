#include <jtflib/mesh/io.h>
#include "../zsw_io.h"

int main(int argc, char *argv[])
{
  if(argc!=2) {
    std::cout << "usage fileType [file_path]" << std::endl;
  }
  std::string file(argv[1]);
  zjucad::matrix::matrix<int> tris;
  zjucad::matrix::matrix<double> vertex;
  zsw::load_stl_binary(file, tris, vertex);
  return 0;
}

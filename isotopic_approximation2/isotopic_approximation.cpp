#include <iostream>
#include <fstream>
#include <zswlib/error_ctrl.h>
#include "isotopic_approximation.h"
#include "basic_op.h"
#include "bound_sphere.h"

using namespace std;

namespace zsw{

  void Approximation::init(const zsw::Scalar &surf_sample_r, const zsw::Scalar &tet_sample_r,
                           const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bi_vertices,
                           const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bo_vertices)
  {
    surf_sample_r_=surf_sample_r;
    tet_sample_r_=tet_sample_r;
    Eigen::Matrix<zsw::Scalar,3,2> bbox;
    calcBBOX(bo_vertices, bbox);
    zsw::Scalar scale = 0.5*(bbox.block<3,1>(0,0) - bbox.block<3,1>(0,1)).norm();
    Eigen::Matrix<zsw::Scalar,3,1> transform = 0.5*(bbox.block<3,1>(0,0)+bbox.block<3,1>(0,1));
    zsw::BoundSphere bs("/home/wegatron/workspace/geometry/data/bound_sphere.obj", scale, transform);
    const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bs_vertices = bs.getVertices();
    std::vector<std::pair<Point, VertexInfo>> vertices;
    for(const Eigen::Matrix<zsw::Scalar,3,1> &v : bs_vertices) {
      vertices.push_back({Point(v[0], v[1], v[2]), VertexInfo(zsw::BBOX_POINT, v, 0.0)});
    }
    for(const Eigen::Matrix<zsw::Scalar,3,1> &v : bi_vertices) {
      vertices.push_back({Point(v[0],v[1],v[2]), VertexInfo(zsw::INNER_POINT, v, 0.0)});
    }
    for(const Eigen::Matrix<zsw::Scalar,3,1> &v : bo_vertices) {
      vertices.push_back({Point(v[0],v[1],v[2]), VertexInfo(zsw::OUTER_POINT, v, 0.0)});
    }
    tw_.reset(new zsw::TriangulationWapper(vertices));
  }

  void Approximation::simpTolerance()
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  }

  void Approximation::mutuallTessellation()
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  }

  void Approximation::simpZeroSurface()
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  }

  void Approximation::writeZeroSurface(const std::string &filepath) const
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  }

#if 0
  void Approximation::writeTetMesh(const std::string &filepath) const
  {
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    const TTds &tds=tw_->getTds();
    ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
        << tds.number_of_vertices() << " float" << std::endl;
    size_t v_index=0;
    CGAL::Unique_hash_map<TTds::Vertex_handle, std::size_t > v_map;
    for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
      ofs << *vit << std::endl;
      v_map[vit]=v_index++;
    }
    const size_t cells_number=tds.number_of_cells();
    std::vector<TTds::Cell_iterator> chds(cells_number);
    size_t ci=0;
    for(TTds::Cell_iterator cit=tds.cells_begin(); ci<cells_number; ++cit) { chds[ci++] = cit; }

    ofs << "CELLS "<< cells_number << " " << cells_number*5 <<std::endl;
    for(size_t i=0; i<cells_number; ++i) {
      ofs << "4 " << v_map[chds[i]->vertex(0)] << " " <<
        v_map[chds[i]->vertex(1)] << " " <<
        v_map[chds[i]->vertex(2)] << " " <<
        v_map[chds[i]->vertex(3)] << std::endl;
    }
    ofs << "CELL_TYPES " << cells_number << std::endl;
    for(size_t i=0; i<cells_number; ++i) {      ofs << "10" << std::endl;    }
    ofs.close();
    tds.print_cells(std::cout, v_map);
  }
#else

  void Approximation::writeTetMesh(const std::string &filepath,
                                   std::vector<std::function<bool(const TTds::Cell_handle)>> ignore_tet_funcs) const
  {
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    const DelaunayTriangulation &delaunay=tw_->getDelaunay();
    ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
        << delaunay.number_of_vertices()+1 << " float" << std::endl;

    size_t v_index=0;
    for(auto vit=delaunay.all_vertices_begin(); vit!=delaunay.all_vertices_end(); ++vit) {
      ofs << *vit << std::endl;
      vit->info().index_=v_index++;
    }

    stringstream ss;
    size_t valid_cells_number=0;
    for(auto cit=delaunay.all_cells_begin(); cit!=delaunay.all_cells_end(); ++cit) {
      bool ignore=false;
      for(auto igf : ignore_tet_funcs) {        if(igf(cit)) { ignore=true; break; }      }
      if(ignore) { continue; }
      ss << "4 " << cit->vertex(0)->info().index_ << " " <<
        cit->vertex(1)->info().index_ << " " <<
        cit->vertex(2)->info().index_ << " " <<
        cit->vertex(3)->info().index_ << std::endl;
      ++valid_cells_number;
    }
    ofs << "CELLS "<< valid_cells_number << " " << valid_cells_number*5 <<std::endl;
    ofs << ss.str();
    ofs << "CELL_TYPES " << valid_cells_number << std::endl;
    for(size_t i=0; i<valid_cells_number; ++i) {      ofs << "10" << std::endl;    }
    ofs.close();
  }

#endif
}

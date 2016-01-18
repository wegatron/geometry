#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_set>
#include <zswlib/error_ctrl.h>
#include "isotopic_approximation.h"
#include "basic_op.h"
#include "bound_sphere.h"
#include "sampling.h"
#include "smoother.h"

using namespace std;

namespace zsw{
  void Approximation::writeJudgePoints(const std::string &filepath) const
  {
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
        << jpts_.size() << " float" << std::endl;
    for(const JudgePoint &jpt : jpts_) {
      ofs << jpt.pt_.transpose() << std::endl;
    }
    ofs << "CELLS " << jpts_.size() << " " << jpts_.size()*2 << std::endl;
    for(size_t i=0; i<jpts_.size(); ++i) {      ofs << "1 " << i << std::endl;    }
    ofs << "CELL_TYPES " << jpts_.size() << std::endl;
    for(size_t i=0; i<jpts_.size(); ++i) { ofs << "1" << std::endl; }
  }

  void Approximation::writeZeroSurface(const std::string &filepath) const
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  }

  void Approximation::writeTetMesh(const std::string &filepath,
                                   std::vector<std::function<bool(const TTds::Cell_handle)>> ignore_tet_funcs) const
  {
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    const TTds &tds=tw_->getTds();
    ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
        << tds.number_of_vertices()-1 << " float" << std::endl;

    CGAL::Unique_hash_map<TTds::Vertex_handle,size_t> v_map;
    size_t v_index=0;
    for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
      if(vit->info().pt_type_==zsw::INVALID_POINT) { continue; }
      ofs << *vit << std::endl;
      v_map[vit] = v_index++;
    }

    stringstream ss;
    size_t valid_cells_number=0;
    for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) {
      bool ignore=false;
      for(auto igf : ignore_tet_funcs) {        if(igf(cit)) { ignore=true; break; }      }
      if(ignore || ignore_invalid(cit)) { continue; }
      ss << "4 " << v_map[cit->vertex(0)] << " " <<
        v_map[cit->vertex(1)] << " " <<
        v_map[cit->vertex(2)] << " " <<
        v_map[cit->vertex(3)] << std::endl;
      ++valid_cells_number;
    }
    ofs << "CELLS "<< valid_cells_number << " " << valid_cells_number*5 <<std::endl;
    ofs << ss.str();
    ofs << "CELL_TYPES " << valid_cells_number << std::endl;
    for(size_t i=0; i<valid_cells_number; ++i) {      ofs << "10" << std::endl;    }
    ofs.close();
  }

  void Approximation::writeAdjcentCells(const std::string &filepath, const TTds::Edge &e) const
  {
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    const TTds &tds=tw_->getTds();
    ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
        << tds.number_of_vertices() << " float" << std::endl;

    CGAL::Unique_hash_map<TTds::Vertex_handle,size_t> v_map;
    size_t v_index=0;
    for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
      ofs << *vit << std::endl;
      v_map[vit] = v_index++;
    }

    std::unordered_set<TTds::Cell_handle, CGAL::Handle_hash_function> cell_set;
    {
      std::list<TTds::Cell_handle> cells;
      tds.incident_cells(e.first->vertex(e.second), std::back_inserter(cells));
      for(auto chd : cells) { cell_set.insert(chd); }
    }

    {
      std::list<TTds::Cell_handle> cells;
      tds.incident_cells(e.first->vertex(e.third), std::back_inserter(cells));
      for(auto chd : cells) { cell_set.insert(chd); }
    }

    ofs << "CELLS " << cell_set.size() << " " << cell_set.size()*5 << std::endl;
    for(auto cit=cell_set.begin(); cit!=cell_set.end(); ++cit) {
      ofs << "4 " << v_map[(*cit)->vertex(0)] << " " <<
        v_map[(*cit)->vertex(1)] << " " <<
        v_map[(*cit)->vertex(2)] << " " <<
        v_map[(*cit)->vertex(3)] << std::endl;
    }
    ofs << "CELL_TYPES " << cell_set.size() << std::endl;
    for(size_t ci=0; ci<cell_set.size(); ++ci) { ofs << "10" << std::endl; }
    ofs.close();
  }

  void Approximation::writeJudgePoints(const std::string &filepath, const vector<const JudgePoint*> &jpts) const
  {
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
        << jpts.size() << " float" << std::endl;
    for(const JudgePoint *jpt : jpts) {
      ofs << jpt->pt_.transpose() << std::endl;
    }
    ofs << "CELLS " << jpts.size() << " " << jpts.size()*2 << std::endl;
    for(size_t i=0; i<jpts.size(); ++i) {      ofs << "1 " << i << std::endl;    }
    ofs << "CELL_TYPES " << jpts.size() << std::endl;
    for(size_t i=0; i<jpts.size(); ++i) { ofs << "1" << std::endl; }
  }

  void writePoints(const std::string &filepath, const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &pts)
  {
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
        << pts.size() << " float" << std::endl;
    for(const Eigen::Matrix<zsw::Scalar,3,1> &pt : pts) {
      ofs << pt.transpose() << std::endl;
    }
    ofs << "CELLS " << pts.size() << " " << pts.size()*2 << std::endl;
    for(size_t i=0; i<pts.size(); ++i) {      ofs << "1 " << i << std::endl;    }
    ofs << "CELL_TYPES " << pts.size() << std::endl;
    for(size_t i=0; i<pts.size(); ++i) { ofs << "1" << std::endl; }
  }

  void Approximation::writeAdjcentCells(const std::string &filepath, const std::vector<Chd> &chds) const
  {
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    const TTds &tds=tw_->getTds();
    ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
        << tds.number_of_vertices() << " float" << std::endl;

    CGAL::Unique_hash_map<TTds::Vertex_handle,size_t> v_map;
    size_t v_index=0;
    for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
      ofs << *vit << std::endl;
      v_map[vit] = v_index++;
    }

    ofs << "CELLS " << chds.size() << " " << chds.size()*5 << std::endl;
    for(auto cit=chds.begin(); cit!=chds.end(); ++cit) {
      ofs << "4 " << v_map[(*cit)->vertex(0)] << " " <<
        v_map[(*cit)->vertex(1)] << " " <<
        v_map[(*cit)->vertex(2)] << " " <<
        v_map[(*cit)->vertex(3)] << std::endl;
    }
    ofs << "CELL_TYPES " << chds.size() << std::endl;
    for(size_t ci=0; ci<chds.size(); ++ci) { ofs << "10" << std::endl; }
    ofs.close();
  }
}

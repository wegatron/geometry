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
      ofs << jpt.pt_cur_.transpose() << std::endl;
    }
    ofs << "CELLS " << jpts_.size() << " " << jpts_.size()*2 << std::endl;
    for(size_t i=0; i<jpts_.size(); ++i) {      ofs << "1 " << i << std::endl;    }
    ofs << "CELL_TYPES " << jpts_.size() << std::endl;
    for(size_t i=0; i<jpts_.size(); ++i) { ofs << "1" << std::endl; }
  }

  void Approximation::writeZeroSurface(const std::string &filepath) const
  {
    const TTds &tds = tw_->getTds();
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    CGAL::Unique_hash_map<TTds::Vertex_handle,size_t> v_map;
    size_t v_index=1;
    for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
      if(vit->info().pt_type_ != zsw::ZERO_POINT) { continue; }
      ofs << "v " << *vit << std::endl;
      v_map[vit] = v_index++;
    }

    TTds::Vertex_handle zero_vit[4];
    for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) {
      size_t zero_cnt = 0;
      size_t inner_cnt = 0;
      Eigen::Matrix<zsw::Scalar,3,1> tmp_vs[4];
      Eigen::Matrix<zsw::Scalar,3,3> mat;
      for(size_t i=0; i<4; ++i) {
        if(cit->vertex(i)->info().pt_type_ == zsw::ZERO_POINT) { zero_vit[zero_cnt] = cit->vertex(i);
          tmp_vs[zero_cnt][0] = zero_vit[zero_cnt]->point()[0];
          tmp_vs[zero_cnt][1] = zero_vit[zero_cnt]->point()[1];
          tmp_vs[zero_cnt][2] = zero_vit[zero_cnt]->point()[2];
          ++zero_cnt;
        }
        else if(cit->vertex(i)->info().pt_type_ == zsw::INNER_POINT) {
          tmp_vs[3][0] = cit->vertex(i)->point()[0];
          tmp_vs[3][1] = cit->vertex(i)->point()[1];
          tmp_vs[3][2] = cit->vertex(i)->point()[2];
          inner_cnt++;
        }
      }
      if(zero_cnt != 3 || inner_cnt !=1) { continue; }
      mat.block<3,1>(0,0) = tmp_vs[1] - tmp_vs[0];
      mat.block<3,1>(0,1) = tmp_vs[2] - tmp_vs[0];
      mat.block<3,1>(0,2) = tmp_vs[3] - tmp_vs[0];
      if(mat.determinant() < 0) {
        ofs << "f " << v_map[zero_vit[0]] << " " << v_map[zero_vit[1]] << " " << v_map[zero_vit[2]] << std::endl;
      } else {
        ofs << "f " << v_map[zero_vit[1]] << " " << v_map[zero_vit[0]] << " " << v_map[zero_vit[2]] << std::endl;
      }
    }
    ofs.close();
  }

    void Approximation::writeTetMesh(const std::string &filepath,
                                     std::vector<std::function<bool(const TTds::Cell_handle)>> ignore_tet_funcs,
                                     const TTds *tds_ptr, bool cur_pts) const
    {
      std::ofstream ofs;
      OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
      const TTds &tds=(tds_ptr==nullptr) ? tw_->getTds() : *tds_ptr;
      ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
          << tds.number_of_vertices()-1 << " float" << std::endl;

      CGAL::Unique_hash_map<TTds::Vertex_handle,size_t> v_map;
      size_t v_index=0;
      for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
        if(vit->info().pt_type_==zsw::INVALID_POINT) { continue; }
        if(cur_pts) { ofs << *vit << std::endl; }
        else { ofs << vit->info().pos_c_.transpose() << std::endl; }
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
#if 0
      ofs << "# Append index info:" << std::endl;
      v_index=0;
      for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
        if(vit->info().pt_type_==zsw::INVALID_POINT) { continue; }
        ofs << "#" << v_index << ":" << vit->info().index_  << " type=" << vit->info().pt_type_ << std::endl;
        ++v_index;
      }
#endif
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
        ofs << jpt->pt_cur_.transpose() << std::endl;
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

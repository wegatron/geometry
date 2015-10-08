#include "triangulation.h"

//#include <algorithm>
#include <fstream>
#include <zswlib/mesh/vtk.h>
#include <zswlib/error_ctrl.h>
#include <zswlib/zsw_clock_c11.h>
#include "sampling.h"

using namespace std;

#define DEBUG 1

namespace zsw
{
  TetMesh::TetMesh(const vector<Point> bz_points, const vector<Point> &bo_points, const vector<Point> &bi_points,
                   const zsw::Scalar sample_dense)
  {
    sample_dense_ = sample_dense;
    // add points
    for(const Point &pt : bz_points) { vertices_.push_back({0, pt, -1, {}}); }
    for(const Point &pt : bo_points) { vertices_.push_back({1, pt, -1, {}}); }
    for(const Point &pt : bi_points) {  vertices_.push_back({-1, pt, -1, {}}); }

    // 3d triangulation
    std::vector<pair<Point, size_t>> tmo; // outer tetmesh and inner tetmesh
    size_t cnt = -1;
    for(const Point &pt : bz_points) {  tmo.push_back(pair<Point, size_t>(pt, ++cnt)); }
    vector<pair<Point, size_t>> tmi(tmo);

    for(const Point &pt : bo_points) { tmo.push_back(pair<Point, size_t>(pt, ++cnt)); }
    for(const Point &pt : bi_points) { tmi.push_back(pair<Point, size_t>(pt, ++cnt)); }

    Delaunay to(tmo.begin(), tmo.end());
    Delaunay ti(tmi.begin(), tmi.end());

    ASSURE(to.is_valid(), exit(__LINE__));
    ASSURE(ti.is_valid(), exit(__LINE__));

    // add tets
    size_t tet_id=0;
    addTets(to, tet_id);
    //std::cerr << "[INFO] stage 1 finished! tet size:" << tet_id << std::endl;
    addTets(ti, tet_id);
    addEdges(ti,to);
  }

  void TetMesh::addEdges(const Delaunay &ti, const Delaunay &to)
  {
    // add edges
    for(Delaunay::Finite_edges_iterator eit=to.finite_edges_begin();
        eit!=to.finite_edges_end(); ++eit) {
      edges_.push_back({eit->second, eit->third});
    }

    for(Delaunay::Finite_edges_iterator eit=ti.finite_edges_begin();
        eit!=ti.finite_edges_end(); ++eit) {
      if(vertices_[eit->second].pt_type_==0 && vertices_[eit->third].pt_type_==0) { continue; }
      edges_.push_back({eit->second, eit->third});
    }

    for(Edge & te : edges_) {
      set<size_t> te_fv;
      set<size_t> tet_id_set;
      vector<size_t> tet_ids;
      for(size_t tet_id : vertices_[te.vind0_].tet_ids_) { tet_id_set.insert(tet_id); }
      for(size_t tet_id : vertices_[te.vind1_].tet_ids_) {
        if(tet_id_set.find(tet_id)!=tet_id_set.end()) { tet_ids.push_back(tet_id); }
      }

      for(size_t tet_id : tet_ids) {
        Tet &tt=tets_[tet_id];
        if(tt.vind0_ != te.vind0_ && tt.vind0_ !=te.vind1_) { te_fv.insert(tt.vind0_); }
        if(tt.vind1_ != te.vind0_ && tt.vind1_ !=te.vind1_) { te_fv.insert(tt.vind1_); }
        if(tt.vind2_ != te.vind0_ && tt.vind2_ !=te.vind1_) { te_fv.insert(tt.vind2_); }
        if(tt.vind3_ != te.vind0_ && tt.vind3_ !=te.vind1_) { te_fv.insert(tt.vind3_); }
      }
      te.fv_.resize(te_fv.size());
      copy(te_fv.begin(), te_fv.end(), te.fv_.begin());
    }
  }
  void TetMesh::addTets(const Delaunay &td, size_t &tet_id)
  {
    for(Delaunay::Finite_cells_iterator cit=td.finite_cells_begin();
        cit!=td.finite_cells_end(); ++cit) {
      Tet tmp_tet = {cit->vertex(0)->info(), cit->vertex(1)->info(), cit->vertex(2)->info(), cit->vertex(3)->info()};
      tets_.push_back(tmp_tet);
      vertices_[tmp_tet.vind0_].tet_ids_.push_back(tet_id);
      vertices_[tmp_tet.vind1_].tet_ids_.push_back(tet_id);
      vertices_[tmp_tet.vind2_].tet_ids_.push_back(tet_id);
      vertices_[tmp_tet.vind3_].tet_ids_.push_back(tet_id);

      // static zsw::common::ClockC11 clock;
      // static int count = 0;
      // count++;
      // if(count%100 == 0) {       std::cerr << "cost:" << clock.time() << std::endl;      }

      //take sample points
      Eigen::Matrix<zsw::Scalar, 3, 4> sample_tet_points;
      sample_tet_points(0,0) = vertices_[tmp_tet.vind0_].pt_[0];
      sample_tet_points(1,0) = vertices_[tmp_tet.vind0_].pt_[1];
      sample_tet_points(2,0) = vertices_[tmp_tet.vind0_].pt_[2];

      sample_tet_points(0,1) = vertices_[tmp_tet.vind1_].pt_[0];
      sample_tet_points(1,1) = vertices_[tmp_tet.vind1_].pt_[1];
      sample_tet_points(2,1) = vertices_[tmp_tet.vind1_].pt_[2];

      sample_tet_points(0,2) = vertices_[tmp_tet.vind2_].pt_[0];
      sample_tet_points(1,2) = vertices_[tmp_tet.vind2_].pt_[1];
      sample_tet_points(2,2) = vertices_[tmp_tet.vind2_].pt_[2];

      sample_tet_points(0,3) = vertices_[tmp_tet.vind3_].pt_[0];
      sample_tet_points(1,3) = vertices_[tmp_tet.vind3_].pt_[1];
      sample_tet_points(2,3) = vertices_[tmp_tet.vind3_].pt_[2];

      // add sample points into sample_points_
      zsw::sampleTet(sample_dense_, sample_tet_points, sample_points_);
      ++tet_id;
    }
  }


    void TetMesh::simplify()
    {
      // collapse ZEdges
      bool collapsable = true;
      while(collapsable) {
        collapsable = false;
        // here the edges_'s vector size should not change
        for(Edge &edge : edges_) {
          for(; vertices_[edge.vind0_].father_!=-1; edge.vind0_=vertices_[edge.vind0_].father_);
          for(; vertices_[edge.vind1_].father_!=-1; edge.vind1_=vertices_[edge.vind1_].father_);
          if(edge.vind0_==edge.vind1_) { continue; }
          if(vertices_[edge.vind0_].pt_type_==0 && vertices_[edge.vind1_].pt_type_==0) {
            if(collapseZEdge(edge)) { collapsable=true; break; }
          }
        }
      }
    }

    bool TetMesh::collapseZEdge(Edge &edge)
    {
      // for(; pf_[edge.first]!=-1; edge.first=pf_[edge.first]);
      // for(; pf_[edge.second]!=-1; edge.second=pf_[edge.second]);

      // // using clean tets
      // std::vector<size_t> tet_ids;
      // for(size_t id : tet_points_[edge.first].tet_ids_) {
      //   if(isValidTet(tets_[id])) {     tet_ids.push_back(id);      }
      // }
      // tet_points_[edge.first].tet_ids_=tet_ids;

      // tet_ids.clear();
      // for(size_t id : tet_points_[edge.second].tet_ids_) {
      //   if(isValidTet(tets_[id])) {     tet_ids.push_back(id);      }
      // }
      // tet_points_[edge.second].tet_ids_=tet_ids;

      // // kernel_region_judge
      // KernelRegionJudger krj;
      // for(size_t id : tet_points_[edge.first].tet_ids_) {
      //   updateKrj(id);
      // }
      // for(size_t id : tet_points_[edge.second].tet_ids_) {
      //   updateKrj(id);
      // }

      // bool isfind=false;
      // Point ret_pt;
      // for(size_t id : tet_points_[edge.first].tet_ids_) {
      //   for(Point pt : tets_[id].sample_points_) {
      //     if(krj.judge(pt)) { isfind=true; ret_pt=pt; break; }
      //   }
      // }
      // vector<int> new_tet_ids;
      // if(!isfind) {
      //   for(size_t id : tet_points_[edge.second].tet_ids_) {
      //     new_tet_ids.push_back(id);
      //     for(Point pt : tets_[id].sample_points_) {
      //       if(krj.judge(pt)) { isfind=true; ret_pt=pt; break; }
      //     }
      //   }
      // }

      // if(isfind) {      // resolve
      //   tet_points_[edge.first].point_type_=0;
      //   tet_points_[edge.first].pt_data_=ret_pt;
      //   pf_[edge.second]=edge.first;
      // }
      // return isfind;
      return false;
    }

    void TetMesh::writeVtk(const std::string &filepath)
    {
//       std::ofstream ofs;
//       OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
//       std::vector<zsw::Scalar> pts_data;
//       std::vector<size_t> tets_data;
//       cleanPoints();
//       for(const TetPoint &tet_point : tet_points_) {
//         pts_data.push_back(tet_point.pt_data_[0]);
//         pts_data.push_back(tet_point.pt_data_[1]);
//         pts_data.push_back(tet_point.pt_data_[2]);
//       }

// #ifdef DEBUG
//       // add sampling points in the tet
//       for(Tet &tet : tets_) {
//         if(!isValidTet(tet)) { continue; }
//         if(tet_points_[tet.pt_ids_[0]].point_type_==tet_points_[tet.pt_ids_[1]].point_type_ &&
//            tet_points_[tet.pt_ids_[1]].point_type_==tet_points_[tet.pt_ids_[2]].point_type_ &&
//            tet_points_[tet.pt_ids_[2]].point_type_==tet_points_[tet.pt_ids_[3]].point_type_) { continue; } // useless tet the same type
//         for(const Point &pt : tet.sample_points_) {
//           pts_data.push_back(pt[0]);
//           pts_data.push_back(pt[1]);
//           pts_data.push_back(pt[2]);
//         }
//       }
// #endif
//       for(Tet &tet : tets_) {
//         if(!isValidTet(tet)) { continue; }
//         if(tet_points_[tet.pt_ids_[0]].point_type_==tet_points_[tet.pt_ids_[1]].point_type_ &&
//            tet_points_[tet.pt_ids_[1]].point_type_==tet_points_[tet.pt_ids_[2]].point_type_ &&
//            tet_points_[tet.pt_ids_[2]].point_type_==tet_points_[tet.pt_ids_[3]].point_type_) { continue; } // useless tet the same type

//         tets_data.push_back(tet.pt_ids_[0]);
//         tets_data.push_back(tet.pt_ids_[1]);
//         tets_data.push_back(tet.pt_ids_[2]);
//         tets_data.push_back(tet.pt_ids_[3]);
//       }
//       size_t n_pts = pts_data.size()/3;
//       size_t n_tets = tets_data.size()/4;
//       std::cout << "[INFO] point size:" << n_pts << std::endl;
//       std::cout << "[INFO] tet size:" << n_tets << std::endl;
//       tet2vtk(ofs, &pts_data[0], n_pts, &tets_data[0], n_tets);
    }

    void TetMesh::writeZeroSetSurface(const std::string &filepath)
    {
      std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
    }
  }

#include "triangulation.h"

//#include <algorithm>
#include <fstream>
#include <zswlib/mesh/vtk.h>
#include <zswlib/error_ctrl.h>


using namespace std;

namespace zsw
{
  TetMesh::TetMesh(const vector<Point> bz_points, const vector<Point> &bo_points, const vector<Point> &bi_points)
  {
    std::vector<pair<Point, size_t>> tmo; // outer tetmesh and inner tetmesh
    size_t cnt = -1;
    for(const Point &pt : bz_points) {
      tmo.push_back(pair<Point, size_t>(pt, ++cnt));
    }
    vector<pair<Point, size_t>> tmi(tmo);

    for(const Point &pt : bo_points) {
      tmo.push_back(pair<Point, size_t>(pt, ++cnt));
    }

    for(const Point &pt : bi_points) {
      tmi.push_back(pair<Point, size_t>(pt, ++cnt));
    }

    Delaunay to(tmo.begin(), tmo.end());
    Delaunay ti(tmi.begin(), tmi.end());
    // fill data
    for(const Point &pt : bz_points) {
      tet_points_.push_back({true, pt, {}});
    }
    for(const Point &pt : bo_points) {
      tet_points_.push_back({false, pt, {}});
    }
    for(const Point &pt : bi_points) {
      tet_points_.push_back({false, pt, {}});
    }

    ASSURE(to.is_valid(), exit(__LINE__));
    ASSURE(ti.is_valid(), exit(__LINE__));

    size_t tet_id = -1;
    for(Delaunay::Finite_cells_iterator cit=to.finite_cells_begin();
        cit!=to.finite_cells_end(); ++cit) {
      size_t p_ids[4] = {cit->vertex(0)->info(), cit->vertex(1)->info(), cit->vertex(2)->info(), cit->vertex(3)->info()};
      tets_.push_back({{p_ids[0], p_ids[1], p_ids[2], p_ids[3]}, {}});
      ++tet_id;
      tet_points_[p_ids[0]].tet_ids_.push_back(tet_id);
      tet_points_[p_ids[1]].tet_ids_.push_back(tet_id);
      tet_points_[p_ids[2]].tet_ids_.push_back(tet_id);
      tet_points_[p_ids[3]].tet_ids_.push_back(tet_id);
      // @todo take sample points
      std::cerr << " @TODO take sample points in the tet!" << std::endl;
    }

    // edges
    for(Delaunay::Finite_edges_iterator eit=to.finite_edges_begin();
        eit!=to.finite_edges_end(); ++eit) {
      // pair<size_t, size_t> edge(eit->vertex(0)->info(), eit->vertex(1)->info());
      // if(edge.first > edge.second) { swap(edge.first, edge.second); }
      edges_.push_back(pair<size_t, size_t>(eit->second, eit->third));
      std::cerr << "@TODO check edge is half edge or not!" << std::endl;
    }

  }

  void TetMesh::simplify()
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  }

  void TetMesh::cleanPoints()
  {
    // set new index
    std::vector<TetPoint> n_tet_points;
    std::vector<size_t> index(tet_points_.size(), -1);
    size_t id = -1;
    for(size_t i=0; i<pf_.size(); ++i) {
      if(pf_[i] == -1) {
        index[i] = ++id;
        n_tet_points.push_back(tet_points_[i]);
      }
    }

    std::vector<size_t> id_map(tet_points_.size(), -1);
    for(size_t i=0; i<pf_.size(); ++i) {
      size_t tmp = pf_[i];
      while(pf_[tmp] != -1) { tmp = pf_[tmp]; }
      id_map[i] = index[tmp];
    }

    // using id_map to update edges and tets
    for(Tet &tet : tets_) {
      for(size_t i=0; i<4; ++i) {
        tet.pt_ids_[i] = id_map[tet.pt_ids_[i]];
      }
    }

    for(std::pair<size_t, size_t> &edge : edges_) {
      edge.first = id_map[edge.first];
      edge.second = id_map[edge.second];
    }

    tet_points_ = n_tet_points;
    pf_.resize(tet_points_.size());
    fill(pf_.begin(), pf_.end(), -1);
  }

  void TetMesh::writeVtk(const std::string &filepath)
  {
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    std::vector<zsw::Scalar> pts_data;
    std::vector<size_t> tets_data;
    cleanPoints();
    for(const TetPoint &tet_point : tet_points_) {
      pts_data.push_back(tet_point.pt_data_[0]);
      pts_data.push_back(tet_point.pt_data_[1]);
      pts_data.push_back(tet_point.pt_data_[2]);
    }
    for(const Tet &tet : tets_) {
      if(tet.pt_ids_[0]==tet.pt_ids_[1] || tet.pt_ids_[1] == tet.pt_ids_[2]  || tet.pt_ids_[2] == tet.pt_ids_[3]) { continue;  }
      tets_data.push_back(tet.pt_ids_[0]+1);
      tets_data.push_back(tet.pt_ids_[1]+1);
      tets_data.push_back(tet.pt_ids_[2]+1);
      tets_data.push_back(tet.pt_ids_[3]+1);
    }
    size_t n_pts = pts_data.size()/3;
    size_t n_tets = tets_data.size()/4;
    tet2vtk(ofs, &pts_data[0], n_pts, &tets_data[0], n_tets);
  }

  void TetMesh::writeZeroSetSurface(const std::string &filepath)
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  }

  bool TetMesh::isValidTet(Tet &tet)
  {
    for(size_t i=0; i<4; ++i) {
      for(; pf_[tet.p_ids_[i]]!=-1; tet.p_ids_[i]=pf_[tet.p_ids_[i]]);
    }
    sort(&(tet.pt_ids_[0]), &(tet.pt_ids_[0])+4);
    return (tet.pt_ids_[0]!=tet.pt_ids_[1]) &&  (tet.pt_ids_[1]!=tet.pt_ids_[2])
      &&  (tet.pt_ids_[2]!=tet.pt_ids_[3]);
  }

  bool TetMesh::isValidEdge(pair<size_t, size_t> &edge)
  {
    for(; pf_[edge.first]!=-1; edge.first=pf_[edge.first]);
    for(; pf_[edge.second]!=-1; edge.second=pf_[edge.second]);
    return edge.first!=edge.second;
  }
}

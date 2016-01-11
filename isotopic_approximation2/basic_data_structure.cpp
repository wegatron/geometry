#include "basic_data_structure.h"

#include <fstream>
#include <unordered_map>

#include <zswlib/error_ctrl.h>

namespace zsw{

  TriangulationWapper::TriangulationWapper(const std::vector<std::pair<Point, VertexInfo>> &vertices)
    : delaunay_triangulation_(vertices.begin(), vertices.end()), tds_(delaunay_triangulation_.tds())  {
    // remove infinite_vertex and cells
    //    auto vh = delaunay_triangulation_.infinite_vertex();
    //    delaunay_triangulation_.remove();
    //    std::list<TTds::Cell_handle> invcells;
    //    tds_.incident_cells(vh, std::back_inserter(invcells));
    //    tds_.delete_cells(invcells.begin(), invcells.end()) ;
    //    tds_.delete_vertex(vh);
  }

  bool TriangulationWapper::isSatisfyLinkCondition(const TTds::Edge &edge) const
  {
    TTds::Facet_circulator fit = tds_.incident_facets(edge);
    TTds::Facet_circulator done(fit);

    do {
      std::unordered_map<size_t, size_t> count_map;
      for(size_t i=0; i<4; ++i) {
        if(i==fit->second) { continue; }

        assert(tds_.is_cell(edge.first));
        assert(fit->first->vertex(i) != TTds::Vertex_handle());
        assert(tds_.is_vertex(fit->first->vertex(i)));
        assert(tds_.is_valid());

        std::vector<Vhd> adj_vhds;
        tds_.adjacent_vertices(fit->first->vertex(i), std::back_inserter(adj_vhds));
        for(Vhd vhd : adj_vhds) {
          auto iter = count_map.find(vhd->info().index_);
          if(iter==count_map.end()) { count_map[vhd->info().index_]=1; }
          else { count_map[vhd->info().index_] = iter->second+1; }
        }
      }
      size_t cnt=0;
      for(auto iter=count_map.begin(); iter!=count_map.end(); ++iter) {
        if(iter->second==3) { ++cnt; }
      }
      if(cnt!=2) { return false; }
    } while(++fit!=done);
    return true;
  }

  void TriangulationWapper::calcBoundTris(const TTds::Edge &edge, std::vector<VertexTriple> &bound_tris,
                                          std::vector<Vhd> &opposite_vs) const
  {
    assert(tds_.is_edge(edge.first,edge.second,edge.third));
    int ind[2]={edge.second, edge.third};
    for(size_t vi=0; vi<2; ++vi) {
      auto cur_vhd=edge.first->vertex(ind[vi]);
      std::list<TTds::Cell_handle> cells;
      tds_.incident_cells(cur_vhd, std::back_inserter(cells));
      for(auto cit : cells) {
        Vhd vhds[3];
        size_t cnt=0;
        for(size_t i=0; i<4; ++i) {
          if(cit->vertex(i)!=edge.first->vertex(ind[0]) && cit->vertex(i)!=edge.first->vertex(ind[1])) {
            vhds[cnt++]=cit->vertex(i);
          }
        }
        if(cnt==3) {
          bound_tris.push_back({vhds[0], vhds[1], vhds[2]});
          opposite_vs.push_back(cur_vhd);
        }
      }
    }
  }

  void TriangulationWapper::collapseEdge(TTds::Edge &edge, Vhd merge_vhd, const Eigen::Matrix<zsw::Scalar,3,1> &pt)
  {
    assert(tds_.is_valid());
    assert(tds_.is_edge(edge.first, edge.second, edge.third));
    Vhd remove_v = (edge.first->vertex(edge.second)==merge_vhd) ? edge.first->vertex(edge.third) : edge.first->vertex(edge.second);
    std::vector<TTds::Cell_handle> hole;
    std::map<VertexTriple, Facet> outer_map;
    makeHole(remove_v, outer_map, hole);
    std::map<VertexTriple, Facet> vt_map(outer_map);
    for(auto oit=outer_map.begin(); oit!=outer_map.end(); ++oit) {
      if(oit->first.first==merge_vhd || oit->first.second==merge_vhd || oit->first.third==merge_vhd) { continue; }
      Chd nw_chd = tds_.create_cell();
      nw_chd->set_vertices(oit->first.first, oit->first.second, oit->first.third, merge_vhd);
      // for vertex triples
      for(size_t ai=0; ai<4; ++ai) {
        VertexTriple tmp_vt;
        Vhd *ptr=&tmp_vt.first;
        for(size_t j=0; j<4; ++j) {
          if(j==ai) { continue; }
          *ptr = nw_chd->vertex(j); ++ptr;
        }
        makeCanonical(tmp_vt);
        auto tmp_it = vt_map.find(tmp_vt);
        if(tmp_it!=vt_map.end()) {
          // glue two faces
          nw_chd->set_neighbor(ai, tmp_it->second.first);
          tmp_it->second.first->set_neighbor(tmp_it->second.second,nw_chd);
        } else { vt_map[tmp_vt]=std::make_pair(nw_chd, ai); }
      }
    }
    tds_.delete_cells(hole.begin(), hole.end());
    tds_.delete_vertex(remove_v);
    assert(tds_.is_valid());
  }

  void TriangulationWapper::makeHole(Vhd vhd, std::map<VertexTriple, Facet> &outer_map,
                                     std::vector<Chd> &hole)
  {
    tds_.incident_cells(vhd, std::back_inserter(hole));
    for(Chd &chd : hole) {
      // find neighbor
      int indv=chd->index(vhd);
      Chd oppo_chd=chd->neighbor(indv);
      VertexTriple tmp_vt;
      Vhd *ptr=&tmp_vt.first;
      for(int vi=0; vi<4; ++vi) { if(vi!=indv) { chd->vertex(vi)->set_cell(oppo_chd); *ptr=chd->vertex(vi); ++ptr; } }
      makeCanonical(tmp_vt);
      // find outer face
      Facet f(oppo_chd, oppo_chd->index(chd));
      outer_map.insert(std::make_pair(tmp_vt,f));
    }
  }

  void TriangulationWapper::insertInEdge(TTds::Edge &edge, const Point &pt)
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  }

  bool TriangulationWapper::isBoundaryEdge(const TTds::Edge &edge) const
  {
    if(!tds_.is_edge(edge.first, edge.second, edge.third)) { return false; }
    zsw::PointType oppo_type;
    if(edge.first->vertex(edge.second)->info().pt_type_==zsw::INNER_POINT
       && edge.first->vertex(edge.third)->info().pt_type_==zsw::INNER_POINT) {
      oppo_type= zsw::OUTER_POINT;
    } else if(edge.first->vertex(edge.second)->info().pt_type_==zsw::OUTER_POINT
              && edge.first->vertex(edge.third)->info().pt_type_==zsw::OUTER_POINT) {
      oppo_type = zsw::INNER_POINT;
    } else { return false; }
    CGAL::Container_from_circulator<TTds::Cell_circulator> cells(tds_.incident_cells(edge));
    bool ret=false;
    for(auto cell : cells) {
      if(cell.vertex(0)->info().pt_type_== oppo_type ||
         cell.vertex(1)->info().pt_type_==oppo_type ||
         cell.vertex(2)->info().pt_type_==oppo_type ||
         cell.vertex(3)->info().pt_type_==oppo_type) {
        ret=true;
        break;
      }
    }
    return ret;
  }

  void TriangulationWapper::writeVertex(const std::string &filepath, const std::vector<Vhd> &vs) const
  {
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
        << vs.size() << "float" << std::endl;
    for(Vhd vhd : vs) {      ofs << vhd->point()[0] << " " << vhd->point()[1] << " " << vhd->point()[2] << std::endl;    }
    ofs << "CELLS " << vs.size() << " " << vs.size()*2 << std::endl;
    for(size_t i=0; i<vs.size(); ++i) {      ofs << "1 " << i << std::endl;    }
    ofs << "CELL_TYPES " <<vs.size() << std::endl;
    for(size_t i=0; i<vs.size(); ++i) { ofs << "1" << std::endl; }
  }

  void TriangulationWapper::test() const
  {
    TTds &tds=tds_;
    std::unordered_map<std::string, TTds::Edge> edge_map;
    for(TTds::Edge_iterator eit=tds.edges_begin();
        eit!=tds.edges_end(); ++eit) {
      if( isBoundaryEdge(*eit)) {
        size_t key[2] = {eit->first->vertex(eit->second)->info().index_,
                         eit->first->vertex(eit->third)->info().index_};
        if(key[0]>key[1]) { std::swap(key[0], key[1]); }
        std::string key_str(std::to_string(key[0]) + "," + std::to_string(key[1]));
        if(edge_map.find(key_str)==edge_map.end()) { edge_map.insert(std::make_pair(key_str, *eit)); }
      }
    }
    std::cerr << "edge size:" << edge_map.size() << std::endl;
    while(!edge_map.empty()) {
      TTds::Edge e = edge_map.begin()->second; edge_map.erase(edge_map.begin());
      if(!tds.is_edge(e.first, e.second, e.third)) { continue; }
      std::cout << "try collapse edge:";
      std::cout << e.first->vertex(e.second)->info().index_ << " " <<
        e.first->vertex(e.third)->info().index_ << std::endl;
      isSatisfyLinkCondition(e);
    }
    std::cerr << "test passed!!!!!" << std::endl;
  }


  bool ignore_invalid(const TTds::Cell_handle cell) {
    return cell->vertex(0)->info().pt_type_==zsw::BBOX_POINT ||
      cell->vertex(1)->info().pt_type_==zsw::BBOX_POINT ||
      cell->vertex(2)->info().pt_type_==zsw::BBOX_POINT ||
      cell->vertex(3)->info().pt_type_==zsw::BBOX_POINT;
  }

  bool ignore_bbox(const TTds::Cell_handle cell) {
    return cell->vertex(0)->info().pt_type_==zsw::BBOX_POINT ||
      cell->vertex(1)->info().pt_type_==zsw::BBOX_POINT ||
      cell->vertex(2)->info().pt_type_==zsw::BBOX_POINT ||
      cell->vertex(3)->info().pt_type_==zsw::BBOX_POINT;
  }

  bool ignore_self_in(const TTds::Cell_handle cell) {
    return cell->vertex(0)->info().pt_type_==zsw::INNER_POINT &&
      cell->vertex(1)->info().pt_type_==zsw::INNER_POINT &&
      cell->vertex(2)->info().pt_type_==zsw::INNER_POINT &&
      cell->vertex(3)->info().pt_type_==zsw::INNER_POINT;
  }

  bool ignore_self_out(const TTds::Cell_handle cell) {
    return cell->vertex(0)->info().pt_type_==zsw::OUTER_POINT &&
      cell->vertex(1)->info().pt_type_==zsw::OUTER_POINT &&
      cell->vertex(2)->info().pt_type_==zsw::OUTER_POINT &&
      cell->vertex(3)->info().pt_type_==zsw::OUTER_POINT;
  }

  std::string cell2str(const TTds::Cell_handle cell)
  {
    std::stringstream ss;
    ss << "cell:" << cell->vertex(0)->info().index_ << " " <<
      cell->vertex(1)->info().index_ << " " <<
      cell->vertex(2)->info().index_ << " " <<
      cell->vertex(3)->info().index_;
    return ss.str();
  }

  void makeCanonical(VertexTriple &t)
  {
    int i = (t.first < t.second) ? 0 : 1;
    if(i==0) {
      i = (t.first < t.third) ? 0 : 2;
    } else {
      i = (t.second < t.third) ? 1 : 2;
    }
    Vhd tmp;
    switch(i){
    case 0: return;
    case 1:
      tmp = t.first;
      t.first = t.second;
      t.second = t.third;
      t.third = tmp;
      return;
    default:
      tmp = t.first;
      t.first = t.third;
      t.third = t.second;
      t.second = tmp;
    }
  }
}

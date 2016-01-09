#include "basic_data_structure.h"

#include <unordered_map>

namespace zsw{

  TriangulationWapper::TriangulationWapper(const std::vector<std::pair<Point, VertexInfo>> &vertices)
    : delaunay_triangulation_(vertices.begin(), vertices.end()), tds_(delaunay_triangulation_.tds())  {
    // remove infinite_vertex and cells
    auto vh = delaunay_triangulation_.infinite_vertex();
    std::list<TTds::Cell_handle> invcells;
    tds_.incident_cells(vh, std::back_inserter(invcells));
    tds_.delete_cells(invcells.begin(), invcells.end()) ;
    tds_.delete_vertex(vh);
  }

  bool TriangulationWapper::isSatisfyLinkCondition(const TTds::Edge &edge)
  {
    TTds::Facet_circulator fit = tds_.incident_facets(edge);
    TTds::Facet_circulator done(fit);
    do {
      std::unordered_map<size_t, size_t> count_map;
      for(size_t i=0; i<4; ++i) {
        if(i==fit->second) { continue; }
        std::list<Vhd> adj_vhds;
        tds_.adjacent_vertices(fit->first->vertex(i), std::back_inserter(adj_vhds));
        for(Vhd vhd : adj_vhds) {
          auto iter = count_map.find(vhd->info().index_);
          if(iter==count_map.end()) {            count_map[vhd->info().index_]=1;          }
          else { count_map[vhd->info().index_] = iter->second+1; }
        }
      }
      size_t cnt=0;
      for(auto iter=count_map.begin(); iter!=count_map.end(); ++iter) {
        if(iter->second==3) { ++cnt; }
      }
      if(cnt!=2) { return false; }
    } while(fit!=done);
    return true;
  }

  void TriangulationWapper::calcBoundTris(const TTds::Edge &edge, std::vector<Fhd> &bound_tris,
                                          std::vector<Vhd> &opposite_vs) const
  {
    for(size_t vi=0; vi<2; ++vi) {
      std::list<TTds::Cell_handle> cells;
      tds_.incident_cells(edge.first->vertex(vi), std::back_inserter(cells));
      for(auto cit : cells) {
        Vhd vhds[3];
        size_t cnt=0;
        for(size_t i=0; i<4; ++i) {
          if(cit->vertex(i)!=edge.first->vertex(0) && cit->vertex(i)!=edge.first->vertex(1)) {
            vhds[cnt++]=cit->vertex(i);
          }
        }
        if(cnt==3) {
          bound_tris.push_back({vhds[0], vhds[1], vhds[2]});
          opposite_vs.push_back(edge.first->vertex(vi));
        }
      }
    }
  }

    void TriangulationWapper::collapseEdge(TTds::Edge &edge, Vhd vhd, const Eigen::Matrix<zsw::Scalar,3,1> &pt)
    {
      std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
    }

    void TriangulationWapper::insertInEdge(TTds::Edge &edge, const Point &pt)
    {
      std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
    }

    bool TriangulationWapper::isBoundaryEdge(const TTds::Edge &edge) const
    {
      zsw::PointType oppo_type;
      edge.first;
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

  }

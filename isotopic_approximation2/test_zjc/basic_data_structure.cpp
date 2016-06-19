#include "basic_data_structure.h"

#include <fstream>
#include <unordered_map>

//#include <zswlib/error_ctrl.h>

namespace zsw{

  TriangulationWapper::TriangulationWapper(const std::vector<std::pair<Point, VertexInfo> > &vertices)
    : delaunay_triangulation_(vertices.begin(), vertices.end()) { }

  bool TriangulationWapper::isSatisfyLinkCondition(const TTds::Edge &edge)
  {
    TTds &tds_ = delaunay_triangulation_.tds();
    TTds::Edge tmp_edge = edge;
    TTds::Edge edge_ori = edge;
    TTds::Facet_circulator fit = tds_.incident_facets(tmp_edge);
    TTds::Facet_circulator done(fit);

    do {
      std::unordered_map<size_t, size_t> count_map;
      for(size_t i=0; i<4; ++i) {
        if(i==fit->second) { continue; }

        TTds::Cell_handle cell=edge.first;
        std::vector<Vhd> adj_vhds;

        std::cerr << "i=" << i << std::endl;
        assert(tds_.is_cell(edge.first));
        assert(cell==edge.first);
        std::cerr << edge.second  << " " << edge.third << std::endl;
        assert(edge==edge_ori);
        assert(fit->first->vertex(i) != TTds::Vertex_handle());
        assert(tds_.is_vertex(fit->first->vertex(i)));
        assert(tds_.is_valid());

        tds_.adjacent_vertices<TTds::False_filter>(fit->first->vertex(i), std::back_inserter(adj_vhds));

        //tds_.adjacent_vertices(fit->first->vertex(i), std::back_inserter(adj_vhds));

        adj_vhds.push_back(fit->first->vertex(i)); // todo delete
        std::cerr << "TODO delete vi!!!!" << std::endl;
//        writeVertex("/home/wegatron/tmp/vertex.vtk", adj_vhds);
        std::cerr << std::endl;
        assert(tds_.is_cell(cell));
        assert(cell==edge.first);
        std::cerr << edge.second  << " " << edge.third << std::endl;
        assert(edge==edge_ori);
        // for(Vhd vhd : adj_vhds) {
        //   auto iter = count_map.find(vhd->info().index_);
        //   if(iter==count_map.end()) {            count_map[vhd->info().index_]=1;          }
        //   else { count_map[vhd->info().index_] = iter->second+1; }
        // }
      }
      // size_t cnt=0;
      // for(auto iter=count_map.begin(); iter!=count_map.end(); ++iter) {
      //   if(iter->second==3) { ++cnt; }
      // }
      // if(cnt!=2) { return false; }
    } while(++fit!=done);
    return true;
  }

//  void TriangulationWapper::calcBoundTris(const TTds::Edge &edge, std::vector<Fhd> &bound_tris,
//                                          std::vector<Vhd> &opposite_vs) const
//  {
//    assert(tds_.is_edge(edge.first,edge.second,edge.third));
//    int ind[2]={edge.second, edge.third};
//    for(size_t vi=0; vi<2; ++vi) {
//      auto cur_vhd=edge.first->vertex(ind[vi]);
//      std::list<TTds::Cell_handle> cells;
//      tds_.incident_cells(cur_vhd, std::back_inserter(cells));
//      for(auto cit : cells) {
//        Vhd vhds[3];
//        size_t cnt=0;
//        for(size_t i=0; i<4; ++i) {
//          if(cit->vertex(i)!=edge.first->vertex(ind[0]) && cit->vertex(i)!=edge.first->vertex(ind[1])) {
//            vhds[cnt++]=cit->vertex(i);
//          }
//        }
//        if(cnt==3) {
//          bound_tris.push_back({vhds[0], vhds[1], vhds[2]});
//          opposite_vs.push_back(cur_vhd);
//        }
//      }
//    }
//  }

//  void TriangulationWapper::collapseEdge(TTds::Edge &edge, Vhd vhd, const Eigen::Matrix<zsw::Scalar,3,1> &pt)
//  {
//    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
//  }

//  void TriangulationWapper::insertInEdge(TTds::Edge &edge, const Point &pt)
//  {
//    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
//  }

//  bool TriangulationWapper::isBoundaryEdge(const TTds::Edge &edge) const
//  {
//    zsw::PointType oppo_type;
//    edge.first;
//    if(edge.first->vertex(edge.second)->info().pt_type_==zsw::INNER_POINT
//       && edge.first->vertex(edge.third)->info().pt_type_==zsw::INNER_POINT) {
//      oppo_type= zsw::OUTER_POINT;
//    } else if(edge.first->vertex(edge.second)->info().pt_type_==zsw::OUTER_POINT
//              && edge.first->vertex(edge.third)->info().pt_type_==zsw::OUTER_POINT) {
//      oppo_type = zsw::INNER_POINT;
//    } else { return false; }
//    CGAL::Container_from_circulator<TTds::Cell_circulator> cells(tds_.incident_cells(edge));
//    bool ret=false;
//    for(auto cell : cells) {
//      if(cell.vertex(0)->info().pt_type_== oppo_type ||
//         cell.vertex(1)->info().pt_type_==oppo_type ||
//         cell.vertex(2)->info().pt_type_==oppo_type ||
//         cell.vertex(3)->info().pt_type_==oppo_type) {
//        ret=true;
//        break;
//      }
//    }
//    return ret;
//  }

//  void TriangulationWapper::writeVertex(const std::string &filepath, const std::vector<Vhd> &vs) const
//  {
//    std::ofstream ofs;
//    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
//    ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
//        << vs.size() << "float" << std::endl;
//    for(Vhd vhd : vs) {      ofs << vhd->point()[0] << " " << vhd->point()[1] << " " << vhd->point()[2] << std::endl;    }
//    ofs << "CELLS " << vs.size() << " " << vs.size()*2 << std::endl;
//    for(size_t i=0; i<vs.size(); ++i) {      ofs << "1 " << i << std::endl;    }
//    ofs << "CELL_TYPES " <<vs.size() << std::endl;
//    for(size_t i=0; i<vs.size(); ++i) { ofs << "1" << std::endl; }
//  }

//  bool ignore_invalid(const TTds::Cell_handle cell) {
//    return cell->vertex(0)->info().pt_type_==zsw::BBOX_POINT ||
//      cell->vertex(1)->info().pt_type_==zsw::BBOX_POINT ||
//      cell->vertex(2)->info().pt_type_==zsw::BBOX_POINT ||
//      cell->vertex(3)->info().pt_type_==zsw::BBOX_POINT;
//  }

//  bool ignore_bbox(const TTds::Cell_handle cell) {
//    return cell->vertex(0)->info().pt_type_==zsw::BBOX_POINT ||
//      cell->vertex(1)->info().pt_type_==zsw::BBOX_POINT ||
//      cell->vertex(2)->info().pt_type_==zsw::BBOX_POINT ||
//      cell->vertex(3)->info().pt_type_==zsw::BBOX_POINT;
//  }

//  bool ignore_self_in(const TTds::Cell_handle cell) {
//    return cell->vertex(0)->info().pt_type_==zsw::INNER_POINT &&
//      cell->vertex(1)->info().pt_type_==zsw::INNER_POINT &&
//      cell->vertex(2)->info().pt_type_==zsw::INNER_POINT &&
//      cell->vertex(3)->info().pt_type_==zsw::INNER_POINT;
//  }

//  bool ignore_self_out(const TTds::Cell_handle cell) {
//    return cell->vertex(0)->info().pt_type_==zsw::OUTER_POINT &&
//      cell->vertex(1)->info().pt_type_==zsw::OUTER_POINT &&
//      cell->vertex(2)->info().pt_type_==zsw::OUTER_POINT &&
//      cell->vertex(3)->info().pt_type_==zsw::OUTER_POINT;
//  }

//  std::string cell2str(const TTds::Cell_handle cell)
//  {
//    std::stringstream ss;
//    ss << "cell:" << cell->vertex(0)->info().index_ << " " <<
//      cell->vertex(1)->info().index_ << " " <<
//      cell->vertex(2)->info().index_ << " " <<
//      cell->vertex(3)->info().index_;
//    return ss.str();
//  }


}
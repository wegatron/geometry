#include "basic_data_structure.h"

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

  bool TriangulationWapper::isSatisfyLinkCondition(const TTds::Edge &edge) const
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  }

  void TriangulationWapper::calcBoundTris(const TTds::Edge &edge, std::vector<Fhd> &bound_tris,
                                          std::vector<Vhd> &opposite_vs) const
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
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

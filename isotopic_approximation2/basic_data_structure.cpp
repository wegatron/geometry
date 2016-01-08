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

  bool TriangulationWapper::isSatisfyLinkCondition(const TTds::Edge &edge)
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  }

  // void TriangulationWapper::constructKernelRegionJudger(const Tds::Edge &edge, KernelRegionJudger &krj)
  // {
  //   std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  // }

  void TriangulationWapper::calcBoundTris(const TTds::Edge &edge, std::vector<Fhd> &bound_tris)
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  }

  void TriangulationWapper::collapseEdge(TTds::Edge &edge, const Point &pt)
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  }

  void TriangulationWapper::insertInEdge(TTds::Edge &edge, const Point &pt)
  {
    std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
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

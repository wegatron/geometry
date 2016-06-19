#ifndef BASIC_DATA_STRUCTURE_H
#define BASIC_DATA_STRUCTURE_H

#include <Eigen/Dense>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

namespace zsw{

    typedef double Scalar;

  enum PointType
  {
    INVALID_POINT=0,
    BBOX_POINT=1,
    ZERO_POINT=2,
    OUTER_POINT=4,
    INNER_POINT=8
  };

  struct VertexInfo
  {
    size_t index_;
    PointType pt_type_;
    // vertex's corresponding position in original space(in refine stage) or deformed space(in simp stage)
    Eigen::Matrix<zsw::Scalar,3,1> pos_c_;
    zsw::Scalar max_dis_; // max distance travel
    VertexInfo() {
      index_=-1;
      pt_type_=INVALID_POINT;
      pos_c_=Eigen::Matrix<zsw::Scalar,3,1>::Zero();
      max_dis_=0.0;
    }
    VertexInfo(size_t index, PointType pt_type,
               Eigen::Matrix<zsw::Scalar,3,1> pos_c,
               zsw::Scalar max_dis) {
      index_=index;
      pt_type_=pt_type;
      pos_c_=pos_c;
      max_dis_=max_dis;
    }
  };

  struct JudgePoint
  {
    Eigen::Matrix<zsw::Scalar,3,1> pt_;
    zsw::Scalar val_exp_;
    zsw::Scalar val_cur_;
  };

  typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
  typedef CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, K>    Vb;
  typedef CGAL::Triangulation_data_structure_3<Vb>                    Tds;
  //Use the Fast_location tag. Default or Compact_location works too.
  typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> DelaunayTriangulation;
  typedef DelaunayTriangulation::Point                                             Point;
  typedef CGAL::Triangulation_data_structure_3<CGAL::Triangulation_hierarchy_vertex_base_3<CGAL::Triangulation_vertex_base_with_info_3<zsw::VertexInfo, CGAL::Epick, CGAL::Triangulation_vertex_base_3<CGAL::Epick, CGAL::Triangulation_ds_vertex_base_3<CGAL::Triangulation_data_structure_3<CGAL::Triangulation_vertex_base_with_info_3<zsw::VertexInfo, CGAL::Epick> > > > > >, CGAL::Triangulation_ds_cell_base_3<void>, CGAL::Sequential_tag>  TTds;
  typedef TTds::Vertex_handle Vhd;
//  typedef CGAL::Triple<Vhd, Vhd, Vhd> Fhd;
  typedef std::pair<Point, VertexInfo> P_Info;

  class TriangulationWapper final
  {
  public:
    TriangulationWapper(const std::vector<std::pair<Point, VertexInfo> > &vertices);
    bool isSatisfyLinkCondition(const TTds::Edge &edge);
//    void calcBoundTris(const TTds::Edge &edge, std::vector<Fhd> &bound_tris, std::vector<Vhd> &opposite_vs) const;
//    void collapseEdge(TTds::Edge &edge, Vhd vhd, const Eigen::Matrix<zsw::Scalar,3,1> &pt);
//    void insertInEdge(TTds::Edge &edge, const Point &pt);
//    bool isBoundaryEdge(const TTds::Edge &edge) const;

//    void writeVertex(const std::string &filepath, const std::vector<Vhd> &vs) const;

    DelaunayTriangulation delaunay_triangulation_;
  };


//  bool ignore_invalid(const TTds::Cell_handle cell);

//  bool ignore_bbox(const TTds::Cell_handle cell);

//  bool ignore_self_in(const TTds::Cell_handle cell);

//  bool ignore_self_out(const TTds::Cell_handle cell);

//  std::string cell2str(const TTds::Cell_handle cell);

}
#endif /* BASIC_DATA_STRUCTURE_H */

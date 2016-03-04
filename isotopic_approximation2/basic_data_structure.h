#ifndef BASIC_DATA_STRUCTURE_H
#define BASIC_DATA_STRUCTURE_H

#include <unordered_set>
#include <Eigen/Dense>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <zswlib/config.h>

namespace zsw{

  enum PointType
  {
    INVALID_POINT=0,
    BBOX_POINT=1,
    ZERO_POINT=2,
    OUTER_POINT=4,
    INNER_POINT=8
  };

  struct Plane
  {
    Eigen::Matrix<zsw::Scalar,3,1> normal_;
    Eigen::Matrix<zsw::Scalar,3,1> v0_;
    zsw::Scalar d_;
  };

  struct VertexInfo
  {
    size_t index_;
    PointType pt_type_;
    Eigen::Matrix<zsw::Scalar,3,1> pos_ori_;
    zsw::Scalar max_dis_; // max distance travel
    VertexInfo() {
      index_=-1;
      pt_type_=INVALID_POINT;
      pos_ori_=Eigen::Matrix<zsw::Scalar,3,1>::Zero();
      max_dis_=0.0;
    }
    VertexInfo(size_t index, PointType pt_type,
               Eigen::Matrix<zsw::Scalar,3,1> pos_ori,
               zsw::Scalar max_dis) {
      index_=index;
      pt_type_=pt_type;
      pos_ori_=pos_ori;
      max_dis_=max_dis;
    }
  };

  struct CellInfo
  {
    size_t v_index_[4];
    bool satisfy_normal_cond_;
    CellInfo() {
      v_index_[0]=v_index_[1]=v_index_[2]=v_index_[3]=-1;
      satisfy_normal_cond_=false;
    }
  };

  struct JudgePoint
  {
    Eigen::Matrix<zsw::Scalar,3,1> ptd_; // point position in deformed space
    Eigen::Matrix<zsw::Scalar,3,1> pto_; // point position in origional space
    const zsw::Scalar val_exp_;
    zsw::Scalar val_cur_;
  };

  typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
  typedef CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, K>    Vb;
  typedef CGAL::Triangulation_cell_base_with_info_3<CellInfo, K> Cb;
  typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                    Tds;

  typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> DelaunayTriangulation;
  typedef DelaunayTriangulation::Point                                             Point;

  typedef CGAL::Triangulation_data_structure_3<CGAL::Triangulation_hierarchy_vertex_base_3<CGAL::Triangulation_vertex_base_with_info_3<zsw::VertexInfo, CGAL::Epick, CGAL::Triangulation_vertex_base_3<CGAL::Epick, CGAL::Triangulation_ds_vertex_base_3<Tds > > > >, Cb, CGAL::Sequential_tag> TTds;

  typedef TTds::Vertex_handle Vhd;
  typedef TTds::Facet Facet;
  typedef TTds::Cell_handle Chd;

  typedef CGAL::Triple<Vhd, Vhd, Vhd> VertexTriple;
  typedef std::pair<Point, VertexInfo> PointData;

  class TriangulationWapper final
  {
  public:
    TriangulationWapper(const std::vector<std::pair<Point, VertexInfo>> &vertices);
    Vhd addPointInDelaunaySafe(const Eigen::Matrix<zsw::Scalar,3,1> &pt,
                               VertexInfo &vertex_info,
                               std::vector<Chd> &chds);

    Vhd addPointInDelaunay(const Eigen::Matrix<zsw::Scalar,3,1> &pt,
                           VertexInfo &vertex_info,
                           std::vector<Chd> &chds);

    bool isBoundaryEdge(const TTds::Edge &edge) const;
    bool isZeroEdge(const TTds::Edge &e) const;
    bool isBZEdge(const TTds::Edge &e) const;
    bool isValidCell(Chd chd) const;
    bool isTolCell(Chd chd) const;
    bool isBBoxInnerCell(Chd chd) const;
    bool isValid() const;

    bool isSatisfyLinkCondition(const TTds::Edge &edge) const;
    void calcBoundTris(const TTds::Edge &edge, std::vector<VertexTriple> &bound_tris, std::vector<Vhd> &opposite_vs) const;
    void calcBoundTrisAdvance(const TTds::Edge &edge, std::vector<VertexTriple> &bound_tris, std::vector<Vhd> &opposite_vs) const;
    void calcAdjZeroSupportPlanes(const TTds::Edge &edge, std::vector<Plane> &adj_zero_support_planes) const;

    // void initCellKeySet(std::unordered_set<std::string> &cell_key_set) const;

    void collapseEdge(TTds::Edge &edge, Vhd vhd, const Eigen::Matrix<zsw::Scalar,3,1> &pt);
    Vhd insertInEdge(TTds::Edge &edge, const Point &pt, const PointType pt_type,
                     TTds &tds);

    const DelaunayTriangulation &getDelaunay() { return delaunay_triangulation_; }
    const TTds &getTds() const { return tds_; }
    TTds &getTds()  { return tds_; }
    void setTds(TTds &tds) { tds_=tds; }

    void makeHole(Vhd vhd, std::map<VertexTriple, Facet> &outer_map,
                  std::vector<Chd> &hole);

    void test() const;
    void writeVertex(const std::string &filepath, const std::vector<Vhd> &vs) const;
  private:
    DelaunayTriangulation delaunay_triangulation_;
    TTds &tds_;
    size_t next_v_id_;
  };

  bool ignore_invalid(const TTds::Cell_handle cell);

  bool ignore_bbox(const TTds::Cell_handle cell);

  bool ignore_self_in(const TTds::Cell_handle cell);

  bool ignore_self_out(const TTds::Cell_handle cell);

  bool ignore_out(const TTds::Cell_handle cell);

  // can construct a tolerance cell or a cell with zero point
  bool isConstructTZCell(const PointType pt_type0, const PointType pt_type1,
                          const PointType pt_type2, const PointType pt_type3,
                          Eigen::Matrix<zsw::Scalar,4,1> &val);

  std::string cell2str(const Chd cell);

  std::string cell2key(const Chd chd);

  std::string edge2key(const TTds::Edge &e);

  void makeCanonical(VertexTriple &t);

  void writeCellsAndPoints(const std::string &filepath,
                           std::vector<Eigen::Matrix<zsw::Scalar,3,4>> cells,
                           std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &pts);
}
#endif /* BASIC_DATA_STRUCTURE_H */

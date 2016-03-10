#include "basic_data_structure.h"

#include <fstream>
#include <unordered_map>
#include <zswlib/error_ctrl.h>

#include "isotopic_debug.h"

namespace zsw{

  TriangulationWapper::TriangulationWapper(const std::vector<std::pair<Point, VertexInfo>> &vertices)
    : delaunay_triangulation_(vertices.begin(), vertices.end()), tds_(delaunay_triangulation_.tds())  {
    next_v_id_=vertices.size();
  }

  Chd findCellLink(const TTds &tds, const size_t index0, const size_t index1)
  {
    for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) {
      size_t cnt=0;
      for(size_t i=0; i<4; ++i) {
        cnt+=(cit->vertex(i)->info().index_==index0);
        cnt+=(cit->vertex(i)->info().index_==index1);
      }
      if(cnt==2) { return cit; }
    }
    std::cerr << "can't find cell link " << index0 << "-" << index1 << std::endl;
    abort();
    return nullptr;
  }

  void writeCells(const std::string &filepath, const Chd *chd_ptr, const size_t n)
  {
    std::unordered_map<size_t, Vhd> v_map;
    for(size_t i=0; i<n; ++i) {
      for(size_t j=0; j<4; ++j) {
        v_map[chd_ptr[i]->vertex(j)->info().index_]=chd_ptr[i]->vertex(j);
      }
    }
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
        << v_map.size() << " float" << std::endl;
    std::unordered_map<size_t,size_t> ind_map;
    size_t ind=0;
    for(auto it=v_map.begin(); it!=v_map.end(); ++it) {
      ind_map.insert(std::make_pair(it->first, ind++));
      ofs << it->second->point()[0] << " " << it->second->point()[1] << " " << it->second->point()[2] << std::endl;
    }
    ofs << "CELLS " << n << " " << n*5 << std::endl;
    for(size_t i=0; i<n; ++i) {
      ofs << "4 " << ind_map[chd_ptr[i]->vertex(0)->info().index_] << " "
          << ind_map[chd_ptr[i]->vertex(1)->info().index_] << " "
          << ind_map[chd_ptr[i]->vertex(2)->info().index_] << " "
          << ind_map[chd_ptr[i]->vertex(3)->info().index_] << std::endl;
    }
    ofs << "CELL_TYPES " << n << std::endl;
    for(size_t ci=0; ci<n; ++ci) { ofs << "10"<< std::endl; }
    ofs.close();
  }

  void writeFace(const std::string &filepath, TTds::Facet_circulator fit)
  {
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS 3 float" << std::endl;
    for(size_t i=0; i<4; ++i) {
      if(i!=fit->second) {
        ofs << fit->first->vertex(i)->point()[0] << " "
            << fit->first->vertex(i)->point()[1] << " "
            << fit->first->vertex(i)->point()[2] << std::endl;
      }
    }
    ofs << "CELLS 1 4\n3 0 1 2\nCELL_TYPES 1\n5" << std::endl;
    ofs.close();
  }

  void debugLinkConditionWriteout(const std::string &output_dir, const TTds &tds, TTds::Facet_circulator fit,
                                  const std::unordered_map<size_t,size_t> &count_map)
  {
    struct CellTriple
    {
      Chd first,second,third;
    };
    writeFace(output_dir+"face.vtk",fit);
    std::vector<CellTriple> cells;
    size_t index[3];
    size_t *ptr=index;
    for(size_t i=0; i<4; ++i) {      if(i!=fit->second) { *ptr=fit->first->vertex(i)->info().index_; ++ptr; }    }
    for(auto iter=count_map.begin(); iter!=count_map.end(); ++iter) {
      if(iter->second==3) {
        CellTriple ct;
        ct.first=findCellLink(tds, iter->first, index[0]);
        ct.second=findCellLink(tds, iter->first, index[1]);
        ct.third=findCellLink(tds, iter->first, index[2]);
        cells.push_back(ct);
      }
    }
    size_t ind=0;
    for(const CellTriple &ct : cells) {
      // write cell
      writeCells(output_dir+std::to_string(ind++)+".vtk", &(ct.first), 3);
    }
  }

  void writeCellWithIndex(const TTds &tds, const std::string &file_path, Vhd vhd0, Vhd vhd1, const std::string &key_str)
  {
    std::cerr << "key_str:" << key_str << std::endl;
    std::cerr << "edge pt:" << vhd0->point()[0] << " " << vhd0->point()[1] << " "
              << vhd0->point()[2] << std::endl;
    std::cerr << "edge pt:" << vhd1->point()[0] << " " << vhd1->point()[1] << " "
              << vhd1->point()[2] << std::endl;
    std::istringstream is(key_str);
    size_t key[2];
    is>>key[0]>>key[1];
    Vhd ovhds[2];
    Vhd *ptr=ovhds;
    for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
      if(vit->info().index_==key[0] || vit->info().index_==key[1]) { *ptr=vit; ++ptr; }
    }

    std::cerr << "o pt:" << ovhds[0]->point()[0] << " " << ovhds[0]->point()[1] << " "
              << ovhds[0]->point()[2] << std::endl;
    std::cerr << "o pt:" << ovhds[1]->point()[0] << " " << ovhds[1]->point()[1] << " "
              << ovhds[1]->point()[2] << std::endl;

    std::vector<Chd> chds;
    int vind[4];
    Chd out_chd;
    assert(!tds.is_cell(vhd0,vhd1,ovhds[0],ovhds[1], out_chd, vind[0], vind[1], vind[2],vind[3]));
    for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) {
      size_t cnt=0;
      for(size_t i=0; i<4; ++i) {
        cnt+=(cit->vertex(i)==vhd0);
        cnt+=(cit->vertex(i)==ovhds[0]);
        cnt+=(cit->vertex(i)==ovhds[1]);
      }
      if(cnt>=3) { chds.push_back(cit); break; }
    }
    for(auto cit=tds.cells_begin(); cit!=tds.cells_end(); ++cit) {
      size_t cnt=0;
      for(size_t i=0; i<4; ++i) {
        cnt+=(cit->vertex(i)==vhd1);
        cnt+=(cit->vertex(i)==ovhds[0]);
        cnt+=(cit->vertex(i)==ovhds[1]);
      }
      if(cnt>=3) { chds.push_back(cit); break; }
    }
    writeCells(file_path, chds.data(), chds.size());
  }

  void writeEdge(const TTds &tds, const std::string &file_path, Vhd vhd0, Vhd vhd1)
  {
    Chd chd;
    int i,j;
    if(tds.is_edge(vhd0, vhd1, chd, i, j)) {      writeCells(file_path, &chd, 1);    }
    else { std::cerr << "is not edge!!!!!!!!" << std::endl; }
  }

  bool TriangulationWapper::isSatisfyLinkCondition(const TTds::Edge &edge) const
  {
    std::unordered_set<std::string> edge_key_set[2];
    const Vhd e_vhd[2] = {edge.first->vertex(edge.second), edge.first->vertex(edge.third)};
    std::vector<Chd> cells;
    for(size_t evi=0; evi<2; ++evi) {
      cells.clear();
      tds_.incident_cells(e_vhd[evi], std::back_inserter(cells));
      for(Chd chd : cells) {
        Vhd tmp_vhds[3];
        size_t cnt=0;
        for(size_t i=0; i<4; ++i) {
          if(chd->vertex(i)!=e_vhd[0] && chd->vertex(i)!=e_vhd[1]) {tmp_vhds[cnt++]=chd->vertex(i);}
        }
        if(cnt!=3) { continue; }
        size_t key[3]={tmp_vhds[0]->info().index_,tmp_vhds[1]->info().index_, tmp_vhds[2]->info().index_};
        std::sort(&key[0],&key[0]+3);
        edge_key_set[evi].insert(std::to_string(key[0])+" "+std::to_string(key[1]));
        edge_key_set[evi].insert(std::to_string(key[0])+" "+std::to_string(key[2]));
        edge_key_set[evi].insert(std::to_string(key[1])+" "+std::to_string(key[2]));
      }
    }
    size_t lk_a_b_cnt=0;
    for(const std::string &key_str : edge_key_set[0]) {
      if(edge_key_set[1].find(key_str)!=edge_key_set[1].end()) { ++lk_a_b_cnt; }
    }

    size_t lk_ab_cnt=0;
    TTds::Cell_circulator cit=tds_.incident_cells(edge);
    TTds::Cell_circulator done(cit);
    do{
      ++lk_ab_cnt;
      ++cit;
    }while(cit!=done);
#if 0
    {
      if(lk_a_b_cnt==lk_ab_cnt) { return true; }
      std::unordered_set<std::string> e_keys;
      TTds::Cell_circulator check_cit=tds_.incident_cells(edge);
      TTds::Cell_circulator check_done(cit);
      size_t inci_i=0;
      do{
        size_t tmp_key[2];
        size_t *ptr=tmp_key;
        for(size_t i=0; i<4; ++i) {
          if(check_cit->vertex(i)==e_vhd[0] || check_cit->vertex(i)==e_vhd[1]) {continue;}
          *ptr=check_cit->vertex(i)->info().index_; ++ptr;
        }
        if(tmp_key[0]>tmp_key[1]) { std::swap(tmp_key[0], tmp_key[1]); }
        e_keys.insert(std::to_string(tmp_key[0])+" "+std::to_string(tmp_key[1]));
        Chd tmp_chd=check_cit;
        writeCells("/home/wegatron/tmp/linkcond_debug/debug_incident"+std::to_string(inci_i++)+".vtk", &tmp_chd, 1);
        ++check_cit;
      } while(check_cit!=check_done);

      for(const std::string &key_str : edge_key_set[0]) {
        if(edge_key_set[1].find(key_str)!=edge_key_set[1].end()) {
          if(e_keys.find(key_str)==e_keys.end()) {
            writeCellWithIndex(tds_, "/home/wegatron/tmp/linkcond_debug/debug.vtk", e_vhd[0], e_vhd[1], key_str);
            assert(tds_.is_edge(edge.first, edge.second, edge.third));
            assert(edge.first->vertex(edge.second)==e_vhd[0] && edge.first->vertex(edge.third)==e_vhd[1]);
          }
        }
      }
    }
#endif
    return lk_a_b_cnt==lk_ab_cnt;
  }

#if 0
  bool TriangulationWapper::isSatisfyLinkCondition(const TTds::Edge &edge) const
  {
    TTds::Facet_circulator fit = tds_.incident_facets(edge);
    TTds::Facet_circulator done(fit);
    do {
      std::unordered_map<size_t, size_t> count_map;
      for(size_t i=0; i<4; ++i) {
        if(i==fit->second) { continue; }
        std::vector<Vhd> adj_vhds;
        tds_.adjacent_vertices(fit->first->vertex(i), std::back_inserter(adj_vhds));
        for(Vhd vhd : adj_vhds) {
          assert(vhd->info().index_!=fit->first->vertex(i)->info().index_);
          auto iter=count_map.find(vhd->info().index_);
          if(iter==count_map.end()) { count_map[vhd->info().index_]=1; }
          else { count_map[vhd->info().index_] = iter->second+1; }
        }
      }
      size_t cnt=0;
      for(auto iter=count_map.begin(); iter!=count_map.end(); ++iter) {
        if(iter->second==3) { ++cnt;}
      }
      if(cnt!=2) {
#if 1
        debugLinkConditionWriteout("/home/wegatron/tmp/linkcond_debug/", tds_, fit, count_map);
        std::cout << __FILE__ << __LINE__ << std::endl;
        abort();
#endif
        return false;
      }
    } while(++fit!=done);
    return true;
  }
#endif

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

  void TriangulationWapper::calcAdjZeroSupportPlanes(const TTds::Edge &edge, std::vector<Plane> &adj_zero_support_planes) const
  {
    std::list<TTds::Cell_handle> tmp_cells[2];
    tds_.incident_cells(edge.first->vertex(edge.second), std::back_inserter(tmp_cells[0]));
    tds_.incident_cells(edge.first->vertex(edge.third), std::back_inserter(tmp_cells[1]));
    std::vector<TTds::Cell_handle> cells; cells.reserve(tmp_cells[0].size()+tmp_cells[1].size());
    std::set<std::string> cell_key;
    for(size_t i=0; i<2; ++i) {
      for(TTds::Cell_handle chd : tmp_cells[i]) {
        std::string key=cell2key(chd);
        if(cell_key.find(key)!=cell_key.end()) { continue; }
        cell_key.insert(key);
        cells.push_back(chd);
      }
    }
    Eigen::Matrix<zsw::Scalar,3,1> in_pts[4];
    Eigen::Matrix<zsw::Scalar,3,1> out_pts[4];
    Plane pl;
    for(auto cit : cells) {
      size_t in_j=0;
      size_t out_j=0;
      for(size_t i=0; i<4; ++i) {
        if(cit->vertex(i)->info().pt_type_==zsw::INNER_POINT) {
          in_pts[in_j][0]=cit->vertex(i)->point()[0];
          in_pts[in_j][1]=cit->vertex(i)->point()[1];
          in_pts[in_j][2]=cit->vertex(i)->point()[2];
          ++in_j;
        } else {
          out_pts[out_j][0]=cit->vertex(i)->point()[0];
          out_pts[out_j][1]=cit->vertex(i)->point()[1];
          out_pts[out_j][2]=cit->vertex(i)->point()[2];
          ++out_j;
        }
      }
      if(in_j+out_j!=4 || in_j==0 || in_j==4) { continue; } // not tol cell
      //calc zero plane
      if(in_j==1) { pl.normal_ = (out_pts[1]-out_pts[0]).cross(out_pts[2]-out_pts[0]); }
      else if(in_j==2) { pl.normal_ = (out_pts[1]-out_pts[0]).cross(in_pts[1]-in_pts[0]); }
      else { pl.normal_ = (in_pts[1]-in_pts[0]).cross(in_pts[2]-in_pts[0]); }
      pl.normal_.normalize();
      if(pl.normal_.dot(out_pts[0]-in_pts[0]) <0 ) { pl.normal_=-pl.normal_; }
      pl.v0_ = (in_pts[0]+out_pts[0])*0.5;
#if 0
      if(!checkZeroPlane(pl, cit)) { abort(); }
#endif
      adj_zero_support_planes.push_back(pl);
    }
  }

  void TriangulationWapper::collapseEdge(TTds::Edge &edge, Vhd merge_vhd, const Eigen::Matrix<zsw::Scalar,3,1> &pt)
  {
    // std::cout << __FILE__ << __LINE__ << std::endl;
    // if(!tds_.is_valid()) { std::cout << __FILE__ << __LINE__ << std::endl; abort(); }
    // std::cout << __FILE__ << __LINE__ << std::endl;
    assert(tds_.is_valid());
    assert(tds_.is_edge(edge.first, edge.second, edge.third));
    merge_vhd->set_point(Point(pt[0],pt[1],pt[2]));
    Vhd remove_v = (edge.first->vertex(edge.second)==merge_vhd) ? edge.first->vertex(edge.third) : edge.first->vertex(edge.second);
    std::vector<TTds::Cell_handle> hole;
    std::map<VertexTriple, std::pair<Facet, CGAL::Orientation>> outer_map;
    makeHole(remove_v, outer_map, hole);
    std::map<VertexTriple, std::pair<Facet, CGAL::Orientation>> vt_map(outer_map);

    // std::cerr << "edge : " << edge.first->vertex(edge.second)->info().index_ << " "
    //           << edge.first->vertex(edge.third)->info().index_  << std::endl;
    // std::cerr << "outer_map size:" << outer_map.size() << std::endl;
    for(auto oit=outer_map.begin(); oit!=outer_map.end(); ++oit) {
      if(oit->first.first==merge_vhd || oit->first.second==merge_vhd || oit->first.third==merge_vhd) { continue; }
      Chd nw_chd = tds_.create_cell();
      CGAL::Orientation o = orientation(oit->first.first->point(), oit->first.second->point(),
                                        oit->first.third->point(), merge_vhd->point());
      if(o == oit->second.second) {
        nw_chd->set_vertices(oit->first.first, oit->first.second, oit->first.third, merge_vhd);
      } else {
        nw_chd->set_vertices(oit->first.second, oit->first.first, oit->first.third, merge_vhd);
      }
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
          Facet tmp_face = tmp_it->second.first;
          nw_chd->set_neighbor(ai, tmp_face.first);
          tmp_face.first->set_neighbor(tmp_face.second,nw_chd);
        } else {
          vt_map[tmp_vt]=std::make_pair(std::make_pair(nw_chd, ai), CGAL::POSITIVE); // positive is invalid infomation
        }
      }
    }
    tds_.delete_cells(hole.begin(), hole.end());
    tds_.delete_vertex(remove_v);
    assert(tds_.is_valid());
    // std::cout << __FILE__ << __LINE__ << std::endl;
    // if(!tds_.is_valid(true)) { std::cout << __FILE__ << __LINE__ << std::endl; abort(); }
  }

  bool TriangulationWapper::isValid() const
  {
    // if(!tds_.is_valid()) { return false; }
    // bool flag=false;
    // for(Chd cit=tds_.cells_begin(); cit!=tds_.cells_end(); ++cit) {
    //   if(cit->vertex(0)->info().pt_type_==zsw::INVALID_POINT ||
    //      cit->vertex(1)->info().pt_type_==zsw::INVALID_POINT ||
    //      cit->vertex(2)->info().pt_type_==zsw::INVALID_POINT ||
    //      cit->vertex(3)->info().pt_type_==zsw::INVALID_POINT) { continue; }
    //   CGAL::Orientation o=orientation(cit->vertex(0)->point(), cit->vertex(1)->point(),
    //                                   cit->vertex(2)->point(), cit->vertex(3)->point());
    //   if(o == CGAL::NEGATIVE) {
    //     std::cout << __FILE__ << __LINE__ << std::endl;
    //     std::cout << cit->vertex(0)->point() << std::endl;
    //     std::cout << cit->vertex(1)->point() << std::endl;
    //     std::cout << cit->vertex(2)->point() << std::endl;
    //     std::cout << cit->vertex(3)->point() << std::endl;
    //     flag = true;
    //   }
    // }
    // if(flag) { abort(); }
    return tds_.is_valid();
  }

void TriangulationWapper::makeHole(Vhd vhd, std::map<VertexTriple, std::pair<Facet, CGAL::Orientation>> &outer_map,
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
    CGAL::Orientation o = orientation(chd->vertex(0)->point(), chd->vertex(1)->point(),
                                      chd->vertex(2)->point(), chd->vertex(3)->point());
    outer_map.insert( std::make_pair(tmp_vt, std::make_pair(f,o)) );
  }
}

  Vhd TriangulationWapper::insertInEdge(TTds::Edge &edge, const Point &pt,
                                        const PointType pt_type, TTds &tds)
  {
    Vhd vhd = tds.insert_in_edge(edge);
    vhd->set_point(pt);
    vhd->info().index_=next_v_id_++;
    vhd->info().pt_type_=pt_type;
    vhd->info().pos_c_<< pt[0], pt[1], pt[2];
    return vhd;
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

  Vhd TriangulationWapper::addPointInDelaunay(
                                                  const Eigen::Matrix<zsw::Scalar,3,1> &pt,
                                                  VertexInfo &vertex_info,
                                                  std::vector<Chd> &chds)
  {
    Point point(pt[0],pt[1], pt[2]);
    Vhd vhd=delaunay_triangulation_.insert(point);
    vertex_info.index_=next_v_id_++;
    vhd->info()=vertex_info;

    chds.clear();
    for(auto cit=tds_.cells_begin(); cit!=tds_.cells_end(); ++cit) {
      if(!isValidCell(cit)) { continue; }
      size_t cur_index[4]={cit->vertex(0)->info().index_, cit->vertex(1)->info().index_,
                           cit->vertex(2)->info().index_, cit->vertex(3)->info().index_};
      std::sort(cur_index, cur_index+4);
      bool flag=false;
      for(size_t i=0; i<4; ++i) {
        if(cur_index[i]==cit->info().v_index_[i]) { continue; }
        cit->info().v_index_[i]=cur_index[i];
        flag=true;
      }
      if(flag) { chds.push_back(cit); }
    }
    return vhd;
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

  bool TriangulationWapper::isZeroEdge(const TTds::Edge &edge) const
  {
    if(!tds_.is_edge(edge.first, edge.second, edge.third)) { return false; }
    return edge.first->vertex(edge.second)->info().pt_type_==zsw::ZERO_POINT
      && edge.first->vertex(edge.third)->info().pt_type_==zsw::ZERO_POINT;
  }

  bool TriangulationWapper::isBZEdge(const TTds::Edge &e) const
  {
    if(!tds_.is_edge(e.first, e.second, e.third)) { return false; }
    PointType pt_type[2]={
      e.first->vertex(e.second)->info().pt_type_,
      e.first->vertex(e.third)->info().pt_type_
    };
    size_t flag[3]={false, false, false};
    flag[0]=(pt_type[0]==zsw::ZERO_POINT || pt_type[1]==zsw::ZERO_POINT);
    flag[1]=(pt_type[0]==zsw::INNER_POINT || pt_type[1]==zsw::INNER_POINT);
    flag[2]=(pt_type[0]==zsw::OUTER_POINT || pt_type[1]==zsw::OUTER_POINT);
    return flag[0]&&(flag[1]||flag[2]);
  }

  bool TriangulationWapper::isValidCell(Chd chd) const
  {
    return  tds_.is_cell(chd)
      && chd->vertex(0)->info().pt_type_!=zsw::INVALID_POINT
      && chd->vertex(1)->info().pt_type_!=zsw::INVALID_POINT
      && chd->vertex(2)->info().pt_type_!=zsw::INVALID_POINT
      && chd->vertex(3)->info().pt_type_!=zsw::INVALID_POINT;
  }

  bool TriangulationWapper::isTolCell(Chd chd) const
  {
    if(!tds_.is_cell(chd)) { return false; }
    size_t i_cnt=0;
    size_t o_cnt=0;
    for(size_t i=0; i<4; ++i) {
      if(chd->vertex(i)->info().pt_type_==zsw::INNER_POINT) { ++i_cnt; }
      else if(chd->vertex(i)->info().pt_type_==zsw::OUTER_POINT || chd->vertex(i)->info().pt_type_==zsw::BBOX_POINT) { ++o_cnt; }
    }
    if(i_cnt>0 && o_cnt>0 && i_cnt+o_cnt==4) { return true; }
  }

  bool TriangulationWapper::isBBoxInnerCell(Chd chd) const
  {
    if(!tds_.is_cell(chd)) { return false; }
    return (chd->vertex(0)->info().pt_type_==zsw::INNER_POINT ||
            chd->vertex(1)->info().pt_type_==zsw::INNER_POINT ||
            chd->vertex(2)->info().pt_type_==zsw::INNER_POINT ||
            chd->vertex(3)->info().pt_type_==zsw::INNER_POINT) &&
      (chd->vertex(0)->info().pt_type_==zsw::BBOX_POINT ||
       chd->vertex(1)->info().pt_type_==zsw::BBOX_POINT ||
       chd->vertex(2)->info().pt_type_==zsw::BBOX_POINT ||
       chd->vertex(3)->info().pt_type_==zsw::BBOX_POINT);
  }

  void TriangulationWapper::resetVertexLastUpdate()
  {
    for(auto it = tds_.vertices_begin(); it != tds_.vertices_end(); ++it) {
      it->info().last_update_ = 0;
    }
  }


  bool isConstructTZCell(const PointType pt_type0, const PointType pt_type1,
                          const PointType pt_type2, const PointType pt_type3,
                          Eigen::Matrix<zsw::Scalar,4,1> &val)
  {
    size_t i_cnt=0;
    size_t o_cnt=0;
    size_t z_cnt=0;
    if(pt_type0==zsw::INNER_POINT) { ++i_cnt; val[0]=-1; }
    else if(pt_type0==zsw::OUTER_POINT) { ++o_cnt; val[0]=1; }
    else if(pt_type0==zsw::ZERO_POINT) { ++z_cnt; val[0]=0; }

    if(pt_type1==zsw::INNER_POINT) { ++i_cnt; val[1]=-1; }
    else if(pt_type1==zsw::OUTER_POINT) { ++o_cnt; val[1]=1; }
    else if(pt_type1==zsw::ZERO_POINT) { ++z_cnt; val[1]=0; }

    if(pt_type2==zsw::INNER_POINT) { ++i_cnt; val[2]=-1; }
    else if(pt_type2==zsw::OUTER_POINT) { ++o_cnt; val[2]=1; }
    else if(pt_type2==zsw::ZERO_POINT) { ++z_cnt; val[2]=0; }

    if(pt_type3==zsw::INNER_POINT) { ++i_cnt; val[3]=-1; }
    else if(pt_type3==zsw::OUTER_POINT) { ++o_cnt; val[3]=1; }
    else if(pt_type3==zsw::ZERO_POINT) { ++z_cnt; val[3]=0; }

    return (z_cnt!=0 && (i_cnt==0 || o_cnt==0)) || (i_cnt>0 && o_cnt>0 && i_cnt+o_cnt==4);
  }

  bool ignore_invalid(const TTds::Cell_handle cell)
  {
    return cell->vertex(0)->info().pt_type_==zsw::INVALID_POINT ||
      cell->vertex(1)->info().pt_type_==zsw::INVALID_POINT ||
      cell->vertex(2)->info().pt_type_==zsw::INVALID_POINT ||
      cell->vertex(3)->info().pt_type_==zsw::INVALID_POINT;
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

  bool ignore_out(const TTds::Cell_handle cell)
  {
    return cell->vertex(0)->info().pt_type_==zsw::OUTER_POINT ||
      cell->vertex(1)->info().pt_type_==zsw::OUTER_POINT ||
      cell->vertex(2)->info().pt_type_==zsw::OUTER_POINT ||
      cell->vertex(3)->info().pt_type_==zsw::OUTER_POINT;
  }

  std::string cell2str(const Chd cell)
  {
    std::stringstream ss;
    ss << "cell:" << cell->vertex(0)->info().index_ << " " <<
      cell->vertex(1)->info().index_ << " " <<
      cell->vertex(2)->info().index_ << " " <<
      cell->vertex(3)->info().index_;
    return ss.str();
  }

  std::string cell2key(const Chd chd)
  {
    size_t index[4]={chd->vertex(0)->info().index_,
                     chd->vertex(1)->info().index_,
                     chd->vertex(2)->info().index_,
                     chd->vertex(3)->info().index_};
    std::sort(&index[0], &index[0]+4);
    std::stringstream ss;
    ss << index[0] << "," << index[1] << ","  << index[2] << "," << index[3];
    return ss.str();
  }

  std::string edge2key(const TTds::Edge &e)
  {
    std::stringstream ss;
    size_t key[2] = {e.first->vertex(e.second)->info().index_,
                     e.first->vertex(e.third)->info().index_};
    if(key[0]>key[1]) { std::swap(key[0], key[1]); }
    return std::to_string(key[0]) + "," + std::to_string(key[1]);
  }

  void makeCanonical(VertexTriple &t)
  {
    std::sort(&t.first, &t.first+3, [](Vhd a, Vhd b){ return a->info().index_<b->info().index_; });
  }

  void writeCellsAndPoints(const std::string &filepath, std::vector<Eigen::Matrix<zsw::Scalar,3,4>> cells,
                           std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &pts)
  {
    std::ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
    const size_t n=cells.size()*4+pts.size();
    ofs << "# vtk DataFile Version 2.0\n TET\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "
        << n << " float" << std::endl;
    for(size_t i=0; i<cells.size(); ++i) {
      ofs << cells[i].block<3,1>(0,0).transpose() << std::endl;
      ofs << cells[i].block<3,1>(0,1).transpose() << std::endl;
      ofs << cells[i].block<3,1>(0,2).transpose() << std::endl;
      ofs << cells[i].block<3,1>(0,3).transpose() << std::endl;
    }
    for(size_t i=0; i<pts.size(); ++i) {
      ofs << pts[i].transpose() << std::endl;
    }
    ofs <<"CELLS " << cells.size()+pts.size() << " " << 5*cells.size()+2*pts.size() << std::endl;
    size_t  pt_n=0;
    for(size_t i=0; i<cells.size(); ++i) {
      ofs << "4 " << pt_n << " " << pt_n+1 << " " << pt_n+2 << " " << pt_n+3 << std::endl;
      pt_n+=4;
    }
    for(size_t i=0; i<pts.size(); ++i) {      ofs << "1 " << pt_n++ << std::endl;    }
    ofs << "CELL_TYPES " << cells.size()+pts.size() << std::endl;
    for(size_t i=0; i<cells.size();++i) {      ofs << "10" << std::endl;    }
    for(size_t i=0; i<pts.size(); ++i) { ofs<<"1"<< std::endl; }
    ofs.close();
  }

  void TriangulationWapper::swapVertex()
  {
    for(auto vit = tds_.vertices_begin(); vit != tds_.vertices_end(); ++vit) {
      if(vit->info().pt_type_ == INVALID_POINT) { continue; }
      Point pos_c(vit->info().pos_c_[0], vit->info().pos_c_[1], vit->info().pos_c_[2]);
      vit->info().pos_c_[0]=vit->point()[0];
      vit->info().pos_c_[1]=vit->point()[1];
      vit->info().pos_c_[2]=vit->point()[2];
      vit->set_point(pos_c);
    }
  }
}

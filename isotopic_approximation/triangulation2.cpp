#include "triangulation2.h"
#include <unordered_set>

#include <zswlib/mesh/vtk.h>
#include <zswlib/error_ctrl.h>
#include <zswlib/const_val.h>

#include "sampling.h"
#include "basic_op.h"
#include "debug.h"

#define ADD_VERTEX(pt_type, pt) do{                     \
    vertices_.push_back({true, pt_type, pt, {}, {}});   \
  }while(0)

#define ADD_TET(v0, v1, v2, v3) do{                                     \
    size_t t_id=tets_.size();                                           \
    vertices_[v0].tet_ids_.push_back(t_id);   vertices_[v1].tet_ids_.push_back(t_id); \
    vertices_[v2].tet_ids_.push_back(t_id);   vertices_[v3].tet_ids_.push_back(t_id); \
    tets_.push_back({true, {v0, v1, v2, v3}});                      \
  }while(0)

#define ADD_EDGE(v0, v1) do{                                            \
    size_t e_id=edges_.size();                                          \
    vertices_[v0].edge_ids_.push_back(e_id); vertices_[v1].edge_ids_.push_back(e_id); \
    edges_.push_back({true, {v0, v1}});                                 \
  }while(0)

#define REUSE_VERTEX(pt, pt_type, v_id) do{     \
    vertices_[v_id].valid_=true;                \
    vertices_[v_id].pt_=pt;                     \
    vertices_[v_id].pt_type_=pt_type;           \
    vertices_[v_id].tet_ids_.clear();           \
    vertices_[v_id].edge_ids_.clear();          \
  }while(0)

#define REUSE_EDGE(v0, v1, e_id) do{                    \
    vertices_[v0].edge_ids_.push_back(e_id);            \
    vertices_[v1].edge_ids_.push_back(e_id);            \
    edges_[e_id].valid_=true;                           \
    edges_[e_id].vid_[0]=v0; edges_[e_id].vid_[1]=v1;   \
  }while(0)

#define REUSE_TET(v0, v1, v2, v3, t_id) do {                            \
    vertices_[v0].tet_ids_.push_back(t_id);   vertices_[v1].tet_ids_.push_back(t_id); \
    vertices_[v2].tet_ids_.push_back(t_id);   vertices_[v3].tet_ids_.push_back(t_id); \
    tets_[t_id] = {true, v0, v1, v2, v3, {}};                          \
  } while(0)

#define CHECK_ADD_EDGE(v0, v1, isnew0, isnew1) do{      \
    if(isnew0 || isnew1) {                              \
      ADD_EDGE(v0, v1);                                 \
    }                                                   \
  }while(0)

#define CHECK_AND_ADD_ZERO_POINT(va, vb, nv, ev_map, isnew) do{         \
    std::pair<size_t, size_t> e= (va<vb) ? std::make_pair(va, vb) : std::make_pair(vb, va); \
    auto itr=ev_map.find(e);                                            \
    if(itr==ev_map.end()) {                                             \
      nv=vertices_.size();                                              \
      ADD_VERTEX(ZERO_POINT, (vertices_[va].pt_+vertices_[vb].pt_)/2.0); \
      ev_map[e]=nv;                                                     \
      isnew=true;                                                       \
    } else { isnew=false; nv=itr->second; }                             \
  }while(0)

#define INTET_CONSTRAINT(int_judge, v0, v1, v2, v3) do{                 \
    int_judge.addConstraint(vertices_[v1].pt_, vertices_[v2].pt_, vertices_[v3].pt_, vertices_[v0].pt_); \
    int_judge.addConstraint(vertices_[v0].pt_, vertices_[v2].pt_, vertices_[v3].pt_, vertices_[v1].pt_); \
    int_judge.addConstraint(vertices_[v0].pt_, vertices_[v1].pt_, vertices_[v3].pt_, vertices_[v2].pt_); \
    int_judge.addConstraint(vertices_[v0].pt_, vertices_[v1].pt_, vertices_[v2].pt_, vertices_[v3].pt_); \
  }while(0)

void zsw::KernelRegionJudger::addConstraint(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> &v1,
                       const Eigen::Matrix<zsw::Scalar,3,1> &v2, const Eigen::Matrix<zsw::Scalar,3,1> &vr)
{
    vec_v0.push_back(v0);
    Eigen::Matrix<zsw::Scalar,3,1> va=v1-v0;
    Eigen::Matrix<zsw::Scalar,3,1> vb=v2-v0;
    Eigen::Matrix<zsw::Scalar,3,1> vn=va.cross(vb);
    vn.normalized();
    if(vn.dot(vr-v0) < 0) { vn=-vn; }
    vec_vn.push_back(vn);
}

bool zsw::KernelRegionJudger::judge(const Eigen::Matrix<zsw::Scalar,3,1> &pt)
{
  for(size_t i=0; i<vec_v0.size(); ++i) {
    if(vec_vn[i].dot(pt-vec_v0[i]) < -zsw::const_val::eps) {
      return false;
    }
  }
  return true;
}

zsw::Triangulation::Triangulation(const zsw::Scalar r, std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bo_pts,
                                  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bi_pts)
{
  assert(bo_pts.size()!=0 && bi_pts.size()!=0);

  Eigen::Matrix<zsw::Scalar,3,2> bbox;
  calcBBOX(bo_pts, bbox);

  size_t pt_id=0;
  std::vector<std::pair<Point, size_t>> tet_points;
  for(const Eigen::Matrix<zsw::Scalar,3,1> &tmp : bo_pts) {
    tet_points.push_back({Point(tmp[0], tmp[1], tmp[2]), pt_id++});
    vertices_.push_back({true, OUTER_POINT, tmp, {}, {}});
  }

  for(const Eigen::Matrix<zsw::Scalar,3,1> &tmp : bi_pts) {
    tet_points.push_back({Point(tmp[0],tmp[1],tmp[2]), pt_id++});
    vertices_.push_back({true, INNER_POINT, tmp, {}, {}});
  }

  // add 8 bbox points
  for(size_t i=0; i<2; ++i) {
    for(size_t j=0; j<2; ++j) {
      for(size_t k=0; k<2; ++k) {
        tet_points.push_back({Point(bbox(0,i), bbox(1,j), bbox(2,k)), pt_id++});
        Eigen::Matrix<zsw::Scalar,3,1> tmp; tmp<<bbox(0,i), bbox(1,j), bbox(2,k);
        vertices_.push_back({true, BBOX_POINT, tmp, {}, {}});
      }
    }
  }

  Delaunay delaunay(tet_points.begin(), tet_points.end());
  init(r, delaunay);
}

void zsw::Triangulation::init(const zsw::Scalar r, Delaunay &delaunay)
{
  size_t t_id=0;
  for(Delaunay::Finite_cells_iterator cit=delaunay.finite_cells_begin();
      cit!=delaunay.finite_cells_end(); ++cit) {
    tets_.push_back({true, {cit->vertex(0)->info(), cit->vertex(1)->info(),
            cit->vertex(2)->info(), cit->vertex(3)->info()}, {}});
    // judge points in tets
    size_t bo_cnt=0, bi_cnt=0;
    Eigen::Matrix<zsw::Scalar,3,4> bo_points;
    Eigen::Matrix<zsw::Scalar,3,4> bi_points;
    Tet &tmp_tet=tets_.back();
    for(size_t tvid : tmp_tet.vid_) {
      if(vertices_[tvid].pt_type_==OUTER_POINT) {
        bo_points.block<3,1>(0,bo_cnt++)=vertices_[tvid].pt_;
      } else if(vertices_[tvid].pt_type_==INNER_POINT) {
        bi_points.block<3,1>(0,bi_cnt++)=vertices_[tvid].pt_;
      }
    }
    if(bo_cnt==3) { sampleTriangle(bo_points.block<3,3>(0,0), r, 1, tmp_tet.jpts_); }
    else if(bi_cnt==3) { sampleTriangle(bi_points.block<3,3>(0,0), r, -1, tmp_tet.jpts_); }
    // vertex info
    vertices_[tmp_tet.vid_[0]].tet_ids_.push_back(t_id);  vertices_[tmp_tet.vid_[1]].tet_ids_.push_back(t_id);
    vertices_[tmp_tet.vid_[2]].tet_ids_.push_back(t_id);  vertices_[tmp_tet.vid_[3]].tet_ids_.push_back(t_id++);
  }

  // edge info
  size_t e_id=0;
  for(Delaunay::Finite_edges_iterator eit=delaunay.finite_edges_begin();
      eit!=delaunay.finite_edges_end(); ++eit) {
    edges_.push_back({true, {eit->first->vertex(eit->second)->info(), eit->first->vertex(eit->third)->info()}});
    Edge &tmp_edge=edges_.back();
    assert(e_id<600000);
    vertices_[tmp_edge.vid_[0]].edge_ids_.push_back(e_id);
    vertices_[tmp_edge.vid_[1]].edge_ids_.push_back(e_id++);
  }
}

bool zsw::Triangulation::linkCondition(const Edge &e) const
{
  std::unordered_set<size_t> fv; // vertex construct a face with edge e
  std::set<size_t> adj_v0; // vertex link e.vid_[0]
  std::set<size_t> adj_v1; // vertex link e.vid_[1]
  for(size_t tid : vertices_[e.vid_[0]].tet_ids_) {
    bool isfv=false;
    for(size_t vid : tets_[tid].vid_) {
      if(vid == e.vid_[1]) { isfv=true; }
      adj_v0.insert(vid);
    }
    if(isfv) {    for(size_t vid : tets_[tid].vid_) {      fv.insert(vid);    }    }
  }
  for(size_t tid : vertices_[e.vid_[1]].tet_ids_) {
    for(size_t vid : tets_[tid].vid_) {
      adj_v1.insert(vid);
    }
  }
  size_t fv_cnt=fv.size();
  if(fv.find(e.vid_[0]) != fv.end()) { --fv_cnt; }
  if(fv.find(e.vid_[1]) != fv.end()) { --fv_cnt; }

  size_t cv_cnt=0;
  {
    std::set<size_t>::iterator it0=adj_v0.begin();
    std::set<size_t>::iterator it1=adj_v1.begin();
    while(it0!=adj_v0.end() && it1!=adj_v1.end()) {
      if(*it0 == *it1) { ++cv_cnt; ++it0; ++it1; }
      else if(*it0>*it1) { ++it1; }
      else { ++it0; }
    }
    if(adj_v0.find(e.vid_[0])!=adj_v0.end() && adj_v1.find(e.vid_[0])!=adj_v1.end()) { --cv_cnt; }
    if(adj_v0.find(e.vid_[1])!=adj_v0.end() && adj_v1.find(e.vid_[1])!=adj_v1.end()) { --cv_cnt; }
  }
  assert(cv_cnt >=fv_cnt);
  return cv_cnt==fv_cnt;
  return false;
}

void zsw::Triangulation::simpTolerance()
{
  for(Edge e : edges_) {
    // if e is invalid or is not bi and bo edge
    if(!e.valid_ || vertices_[e.vid_[0]].pt_type_!=vertices_[e.vid_[1]].pt_type_ ||
       (vertices_[e.vid_[0]].pt_type_!=OUTER_POINT && vertices_[e.vid_[0]].pt_type_!=INNER_POINT) ) {
      continue; }

    // link condition
    if(!linkCondition(e)) { continue; }

    std::unordered_set<size_t> tet_ids;
    for(size_t tid : vertices_[e.vid_[0]].tet_ids_) { tet_ids.insert(tid); }
    for(size_t tid : vertices_[e.vid_[1]].tet_ids_) { tet_ids.insert(tid); }

    KernelRegionJudger krj;
    std::list<JudgePoint> all_jpts;
    std::list<Eigen::Matrix<size_t,3,1>> bound_tris;
    std::list<size_t> debug_tet_ids;
    for(size_t tid : tet_ids) {
      all_jpts.splice(all_jpts.begin(), tets_[tid].jpts_);

      // add kernel region constraint
      size_t vcnt=0;
      Eigen::Matrix<size_t,4,1> tmp_bound_tri;
      for(size_t tvid : tets_[tid].vid_) {
        if(tvid == e.vid_[0] || tvid==e.vid_[1]) {          tmp_bound_tri[3]=tvid;        }
        else {          tmp_bound_tri[vcnt++]=tvid;        }
      }
      if(vcnt==3) {
        debug_tet_ids.push_back(tid);
        krj.addConstraint(vertices_[tmp_bound_tri[0]].pt_, vertices_[tmp_bound_tri[1]].pt_,
                          vertices_[tmp_bound_tri[2]].pt_, vertices_[tmp_bound_tri[3]].pt_);
        bound_tris.push_back(tmp_bound_tri.block<3,1>(0,0));
      }
    }

    // find the candicate points
    std::vector<JudgePoint> candicate_pts;
    for(const JudgePoint &jpt : all_jpts) {
      // judge if jpt in kernel region
      if(krj.judge(jpt.pt_)) { candicate_pts.push_back(jpt); }
    }

    std::cout << "edge:" << e.vid_[0] << " : "  << e.vid_[1] << std::endl;
    std::cout << "candicate_pts:" << candicate_pts.size() << std::endl;
    static int debug_cnt=0;
    if(++debug_cnt == 2) {
      for(size_t dt_id : vertices_[28].tet_ids_) {
        writeTet("/home/wegatron/tmp/simp_tol/dt_id"+std::to_string(dt_id)+".vtk", dt_id);
      }
      writeJudgePoints("/home/wegatron/tmp/simp_tol/judgepts", all_jpts);
      for(size_t dtid : debug_tet_ids) {      writeTet("/home/wegatron/tmp/simp_tol/dbt"+std::to_string(dtid)+".vtk", dtid);    }
      return;
    }

    // find the best point in the candicate_pts
    std::sort(candicate_pts.begin(), candicate_pts.end(),
              [](const JudgePoint &a, const JudgePoint &b){ return fabs(a.val_cur_-a.val_exp_)>fabs(b.val_cur_-b.val_exp_); });
    const JudgePoint *merge_point_ptr=nullptr;
    for(const JudgePoint &jpt : candicate_pts) {
      std::cout << __FILE__ << __LINE__ << std::endl;
      if(testCollapse(e, vertices_[e.vid_[0]].pt_type_, jpt.pt_, bound_tris, all_jpts)) { merge_point_ptr=&jpt; break; }
    }

    if(merge_point_ptr != nullptr) {
      std::cout << "collapse edge!!!" << std::endl;
      edgeCollapse(e, vertices_[e.vid_[0]].pt_type_, bound_tris, merge_point_ptr->pt_, all_jpts);
    } else { std::cout << "no poper merge point!"<< std::endl;    }
  }
}

void zsw::Triangulation::addZeroPoints(std::map<std::pair<size_t,size_t>, size_t, PairCompFunc> &ev_map)
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
}

bool pairComp(const std::pair<size_t, size_t> &a, const std::pair<size_t, size_t> &b)
{
  return a.first<b.first || a.second<b.second;
}

void zsw::Triangulation::mutualTessellation()
{
  // add zero point with edge
  std::map<std::pair<size_t,size_t>, size_t, PairCompFunc> ev_map(pairComp); // bi bo edge to vertex index map
  //addZeroPoints(ev_map);
  size_t tet_size=tets_.size();
  for(size_t t_id=0; t_id<tet_size; ++t_id) {
    Tet &tet=tets_[t_id];
    if(!tet.valid_) { continue; }
    size_t vo_cnt=0, vi_cnt=0;
    size_t vo[4], vi[4];
    for(size_t vid : tet.vid_) {
      if(vertices_[vid].pt_type_ == OUTER_POINT) { vo[vo_cnt++]=vid; }
      else if(vertices_[vid].pt_type_ == INNER_POINT){ vi[vi_cnt++]=vid; }
    }
    // bo : bi = 3 : 1
    if(vo_cnt==3 && vi_cnt==1) {
      tessellation3v1(vo[0],vo[1],vo[2],vi[0],tet,ev_map);
    }
    // bo : bi = 2 : 2
    else if(vo_cnt==2 && vo_cnt==vi_cnt) {
      tessellation2v2(vo[0],vo[1],vi[0],vi[1],tet, ev_map);
    }
    // bo : bi = 1 : 3
    else if(vo_cnt==1 && vi_cnt==3) {
      tessellation1v3(vo[0],vi[0],vi[1],vi[2],tet, ev_map);
    }
  }
}

void zsw::Triangulation::tessellation3v1(const size_t vo_0, const size_t vo_1,
                                         const size_t vo_2, const size_t vi_0,
                                         Tet &tet,
                                         std::map<std::pair<size_t,size_t>, size_t, PairCompFunc> &ev_map)
{
  // invalid old tet and edges
  invalidTet(tet);

  size_t inv_cnt=0;
  size_t inv_e_ids[3]={-1, -1, -1};
  for(size_t e_id : vertices_[vi_0].edge_ids_) {
    if(edges_[e_id].vid_[0]==vo_0 || edges_[e_id].vid_[0]==vo_1 || edges_[e_id].vid_[0]==vo_2
       || edges_[e_id].vid_[1]==vo_0 || edges_[e_id].vid_[1]==vo_1 || edges_[e_id].vid_[1]==vo_2)
      {
        inv_e_ids[inv_cnt]=e_id;
        if(++inv_cnt==3) { break; }
      }
  }
  assert(inv_cnt<=3);
  for(size_t ei=0; ei<inv_cnt; ++ei) {    invalidEdge(inv_e_ids[ei]);  }

  // add vertices
  size_t nv0, nv1, nv2;
  bool isnew0=false, isnew1=false, isnew2=false; // if nv* is new vertex
  CHECK_AND_ADD_ZERO_POINT(vo_0, vi_0, nv0, ev_map, isnew0);
  CHECK_AND_ADD_ZERO_POINT(vo_1, vi_0, nv1, ev_map, isnew1);
  CHECK_AND_ADD_ZERO_POINT(vo_2, vi_0, nv2, ev_map, isnew2);

  // add tets
  ADD_TET(vi_0, nv0, nv1, nv2);
  ADD_TET(nv0, vo_0, vo_2, vo_1);
  ADD_TET(nv1, nv0, vo_2, vo_1);
  ADD_TET(nv2, nv1, nv0, vo_2);

  // add edges
  CHECK_ADD_EDGE(nv0, vi_0, isnew0, false);  CHECK_ADD_EDGE(nv0, nv1, isnew0, isnew1);
  CHECK_ADD_EDGE(nv0, nv2, isnew0, isnew2);  CHECK_ADD_EDGE(nv0, vo_0, isnew0, false);
  CHECK_ADD_EDGE(nv0, vo_1, isnew0, false);  CHECK_ADD_EDGE(nv0, vo_2, isnew0, false);

  CHECK_ADD_EDGE(nv1, nv2, isnew1, isnew2);  CHECK_ADD_EDGE(nv1, vi_0, isnew1, false);
  CHECK_ADD_EDGE(nv1, vo_1, isnew1, false);  CHECK_ADD_EDGE(nv1, vo_2, isnew1, false);

  CHECK_ADD_EDGE(nv2, vi_0, isnew2, false);  CHECK_ADD_EDGE(nv2, vo_2, isnew2, false);
}

void zsw::Triangulation::tessellation2v2(const size_t vo_0, const size_t vo_1,
                                         const size_t vi_0, const size_t vi_1,
                                         Tet &tet,
                                         std::map<std::pair<size_t,size_t>, size_t, PairCompFunc> &ev_map)
{
  // invalid old tet and edges
  invalidTet(tet);
  tet.valid_=false;
  size_t inv_e_ids[4]={-1,-1,-1,-1};
  size_t inv_cnt=0;
  for(size_t e_id : vertices_[vo_0].edge_ids_) {
    if(edges_[e_id].vid_[0]==vi_0 || edges_[e_id].vid_[0]==vi_1
       || edges_[e_id].vid_[1]==vi_0 || edges_[e_id].vid_[1]==vi_1) {
      inv_e_ids[inv_cnt]=e_id;
      if(++inv_cnt==2) { break; }
    }
  }
  for(size_t e_id : vertices_[vo_1].edge_ids_) {
    if(edges_[e_id].vid_[0]==vi_0 || edges_[e_id].vid_[0]==vi_1
       || edges_[e_id].vid_[1]==vi_0 || edges_[e_id].vid_[1]==vi_1) {
      inv_e_ids[inv_cnt]=e_id;
      if(++inv_cnt==4) { break; }
    }
  }

  assert(inv_cnt<=4);
  for(size_t ei=0; ei<inv_cnt; ++ei) {    invalidEdge(inv_e_ids[ei]);  }
  // add vertices
  size_t nv0, nv1, nv2, nv3;
  bool isnew0=false, isnew1=false, isnew2=false, isnew3=false;
  CHECK_AND_ADD_ZERO_POINT(vo_0, vi_0, nv0, ev_map, isnew0);
  CHECK_AND_ADD_ZERO_POINT(vo_0, vi_1, nv1, ev_map, isnew1);
  CHECK_AND_ADD_ZERO_POINT(vo_1, vi_0, nv2, ev_map, isnew2);
  CHECK_AND_ADD_ZERO_POINT(vo_1, vi_1, nv3, ev_map, isnew3);

  // add tet
  ADD_TET(nv0, nv1, nv2, vo_0);  ADD_TET(nv0, nv1, nv2, vi_0);
  ADD_TET(nv1, nv2, vi_0, vi_1);  ADD_TET(nv1, nv2, nv3, vi_1);
  ADD_TET(nv2, nv1, vo_0, vo_1);  ADD_TET(nv2, nv1, nv3, vo_1);

  // add edge
  CHECK_ADD_EDGE(nv0, vo_0, isnew0, false); CHECK_ADD_EDGE(nv0, nv1, isnew0, isnew1);
  CHECK_ADD_EDGE(nv0, nv2, isnew0, isnew2); CHECK_ADD_EDGE(nv0, vi_0, isnew0, false);

  CHECK_ADD_EDGE(nv1, vo_0, isnew1, false); CHECK_ADD_EDGE(nv1, vi_0, isnew1, false);
  CHECK_ADD_EDGE(nv1, vi_1, isnew1, false); CHECK_ADD_EDGE(nv1, vo_1, isnew1, false);
  CHECK_ADD_EDGE(nv1, nv2, isnew1, isnew2); CHECK_ADD_EDGE(nv1, nv3, isnew1, isnew3);

  CHECK_ADD_EDGE(nv2, vo_0, isnew2, false); CHECK_ADD_EDGE(nv2, vo_1, isnew2, false);
  CHECK_ADD_EDGE(nv2, vi_0, isnew2, false); CHECK_ADD_EDGE(nv2, vi_1, isnew2, false);
  CHECK_ADD_EDGE(nv2, nv3, isnew2, isnew3);

  CHECK_ADD_EDGE(nv3, vo_1, isnew3, false); CHECK_ADD_EDGE(nv3, vi_1, isnew3, false);
}

void zsw::Triangulation::tessellation1v3(const size_t vo_0, const size_t vi_0,
                                         const size_t vi_1, const size_t vi_2,
                                         Tet &tet,
                                         std::map<std::pair<size_t,size_t>, size_t, PairCompFunc> &ev_map)
{
  // invalid old tet and edges
  invalidTet(tet);
  size_t inv_e_ids[3]={-1,-1,-1};
  size_t inv_cnt=0;
  for(size_t e_id : vertices_[vo_0].edge_ids_) {
    if(edges_[e_id].vid_[0]==vi_0 || edges_[e_id].vid_[0]==vi_1 || edges_[e_id].vid_[0]==vi_2
       || edges_[e_id].vid_[1]==vi_0 || edges_[e_id].vid_[1]==vi_0 || edges_[e_id].vid_[1]==vi_0) {
      inv_e_ids[inv_cnt]=e_id;
      if(++inv_cnt==3) { break; }
    }
  }

  assert(inv_cnt<=3);
  for(size_t ei=0; ei<inv_cnt; ++ei) { invalidEdge(inv_e_ids[ei]); }

  // add vertices
  size_t nv0, nv1, nv2;
  bool isnew0, isnew1, isnew2;
  CHECK_AND_ADD_ZERO_POINT(vo_0, vi_0, nv0, ev_map, isnew0);
  CHECK_AND_ADD_ZERO_POINT(vo_0, vi_1, nv1, ev_map, isnew1);
  CHECK_AND_ADD_ZERO_POINT(vo_0, vi_2, nv2, ev_map, isnew2);

  // add tets
  ADD_TET(vo_0, nv0, nv1, nv2);  ADD_TET(nv0, vi_0, vi_1, vi_2);
  ADD_TET(nv1, nv0, vi_1, vi_2);  ADD_TET(nv2, nv0, nv1, vi_2);

  // add edge
  CHECK_ADD_EDGE(nv0, nv1, isnew0, isnew1);  CHECK_ADD_EDGE(nv0, nv2, isnew0, isnew2);
  CHECK_ADD_EDGE(nv0, vo_0, isnew0, false);  CHECK_ADD_EDGE(nv0, vi_0, isnew0, false);
  CHECK_ADD_EDGE(nv0, vi_1, isnew0, false);  CHECK_ADD_EDGE(nv0, vi_2, isnew0, false);

  CHECK_ADD_EDGE(nv1, nv2, isnew1, isnew2);  CHECK_ADD_EDGE(nv1, vo_0, isnew1, false);
  CHECK_ADD_EDGE(nv1, vi_1, isnew1, false); CHECK_ADD_EDGE(nv1, vi_2, isnew1, false);

  CHECK_ADD_EDGE(nv2, vo_0, isnew2, false); CHECK_ADD_EDGE(nv2, vi_2, isnew2, false);
}

void zsw::Triangulation::invalidEdge(const size_t e_id)
{
  Edge &e=edges_[e_id];
  e.valid_=false;
  if(vertices_[e.vid_[0]].valid_) {
    vertices_[e.vid_[0]].edge_ids_.remove_if([e_id](const size_t &val){ return val==e_id; });
  }
  if(vertices_[e.vid_[1]].valid_) {
    vertices_[e.vid_[1]].edge_ids_.remove_if([e_id](const size_t &val){ return val==e_id; });
  }
}

void zsw::Triangulation::invalidTet(Tet &tet)
{
  tet.valid_=false;
  for(size_t v_id : tet.vid_) {
    if(vertices_[v_id].valid_) {
      vertices_[v_id].tet_ids_.remove_if([&tet, this](const size_t &t_id){ return &tet==&(this->tets_[t_id]); });
    }
  }
}

void zsw::Triangulation::writeTetMesh(const std::string &filepath, size_t mask) const
{
  std::ofstream ofs;
  OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
  std::vector<zsw::Scalar> pts_data;
  std::vector<size_t> tets_data;
  for(const Vertex &v : vertices_) {
    pts_data.push_back(v.pt_[0]);
    pts_data.push_back(v.pt_[1]);
    pts_data.push_back(v.pt_[2]);
  }

  for(const Tet &tet : tets_) {
    if(!tet.valid_) { continue; }
    bool ignore=false;
    for(size_t i=0; i<4;++i) {
      if(vertices_[tet.vid_[i]].pt_type_ & mask) { ignore=true; break; }
    }
    if(ignore) { continue; }
    tets_data.push_back(tet.vid_[0]);
    tets_data.push_back(tet.vid_[1]);
    tets_data.push_back(tet.vid_[2]);
    tets_data.push_back(tet.vid_[3]);
  }
  size_t n_pts = pts_data.size()/3;
  size_t n_tets = tets_data.size()/4;
  std::cout << "[INFO] point size:" << n_pts << std::endl;
  std::cout << "[INFO] tet size:" << n_tets << std::endl;
  tet2vtk(ofs, &pts_data[0], n_pts, &tets_data[0], n_tets);
}

void zsw::Triangulation::writeSurface(const std::string &filepath, PointType pt_tyte) const
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
}

bool zsw::Triangulation::testCollapse(const Edge &e, const PointType pt_type,
                                      const Eigen::Matrix<zsw::Scalar,3,1> &pt,
                                      const std::list<Eigen::Matrix<size_t,3,1>> &bound_tris,
                                      const std::list<JudgePoint> &jpts) const
{
  assert(pt_type!=zsw::BBOX_POINT);
  zsw::Scalar pt_val=0.0;
  if(pt_type==zsw::INNER_POINT) { pt_val=-1.0; }
  else if(pt_type==zsw::ZERO_POINT) { pt_val=0.0; }
  else { pt_val=1.0; }
  for(const Eigen::Matrix<size_t,3,1> &b_tr  : bound_tris) {
    Eigen::Matrix<zsw::Scalar,3,3> A;
    A.block<3,1>(0,0) = vertices_[b_tr[0]].pt_-pt;
    A.block<3,1>(0,1) = vertices_[b_tr[1]].pt_-pt;
    A.block<3,1>(0,2) = vertices_[b_tr[2]].pt_-pt;
    // on the same plane of one bound_tri
    if(fabs(A.determinant())<zsw::const_val::eps) { return false; }
    Eigen::PartialPivLU<Eigen::Matrix<zsw::Scalar,3,3>> pplu;
    pplu.compute(A);

    Eigen::Matrix<zsw::Scalar,3,1> nv;
    for(size_t i=0; i<3; ++i) {
      if(vertices_[b_tr[i]].pt_type_==zsw::INNER_POINT) { nv[i]=-1.0-pt_val; }
      else if(vertices_[b_tr[i]].pt_type_==zsw::ZERO_POINT) { nv[i]=0.0-pt_val; }
      else { nv[i]=1.0-pt_val; }
    }
    for(const JudgePoint &jpt : jpts) {
      //check if jpt in the tet and if jpt's error is in tolerance
      Eigen::Matrix<zsw::Scalar,3,1> ans = pplu.solve(jpt.pt_-pt);
      if(ans.squaredNorm()>1 || ans[0]<0 || ans[1]<0 || ans[2]<0) { continue; } // not in tet
      assert((A*ans-(jpt.pt_-pt)).squaredNorm()<zsw::const_val::eps);
      zsw::Scalar cur_val=pt_val+ans.dot(nv);
      if(fabs(cur_val-jpt.val_exp_) > 1.0) {
        return false;
      }
    }
  }
  return true;
}

void zsw::Triangulation::edgeCollapse(Edge &e, const PointType pt_type,
                                      const std::list<Eigen::Matrix<size_t,3,1>> &bound_tris,
                                      const Eigen::Matrix<zsw::Scalar,3,1> &pt,
                                      std::list<JudgePoint> &jpts)
{
  // invalid old tets and edges and vertex
  std::set<size_t> inv_tet_ids, inv_edge_ids;
  vertices_[e.vid_[0]].valid_=false;  vertices_[e.vid_[1]].valid_=false;
  for(size_t t_id : vertices_[e.vid_[0]].tet_ids_) { inv_tet_ids.insert(t_id); }
  for(size_t t_id : vertices_[e.vid_[1]].tet_ids_) { inv_tet_ids.insert(t_id); }
  for(size_t e_id : vertices_[e.vid_[0]].edge_ids_) { inv_edge_ids.insert(e_id); }
  for(size_t e_id : vertices_[e.vid_[1]].edge_ids_) { inv_edge_ids.insert(e_id); }

  for(size_t t_id : inv_tet_ids) { invalidTet(tets_[t_id]); }
  for(size_t e_id : inv_edge_ids) { invalidEdge(e_id); }

#if 0
  for(size_t ti=0; ti<tets_.size(); ++ti) {
    if(!tets_[ti].valid_) { continue; }
    for(size_t vi=0; vi<4; ++vi) {
      if(tets_[ti].vid_[vi]==e.vid_[0] || tets_[ti].vid_[vi]==e.vid_[1]) {
        std::cerr << "invalid tet failed!" << std::endl;
      }
    }
  }
#endif

  // add new vertex
  std::cerr << "pt_type:" << pt_type << std::endl;
  REUSE_VERTEX(pt, pt_type, e.vid_[0]);

  // add new tets
  zsw::Scalar pt_val=0.0;
  if(pt_type==zsw::INNER_POINT) { pt_val=-1.0; }
  else if(pt_type==zsw::ZERO_POINT) { pt_val=0.0; }
  else { pt_val=1.0; }
  assert(inv_tet_ids.size()>=bound_tris.size());
  auto tet_itr = inv_tet_ids.begin();
  for(const Eigen::Matrix<size_t,3,1> &b_tr : bound_tris) {
    std::cerr << "reuse tet:" << b_tr.transpose() << std::endl;
    REUSE_TET(e.vid_[0], b_tr[0], b_tr[1], b_tr[2], *tet_itr);
    // check if jpt in tet
    Eigen::Matrix<zsw::Scalar,3,3> A;
    A.block<3,1>(0,0) = vertices_[b_tr[0]].pt_-pt;
    A.block<3,1>(0,1) = vertices_[b_tr[1]].pt_-pt;
    A.block<3,1>(0,2) = vertices_[b_tr[2]].pt_-pt;
    Eigen::PartialPivLU<Eigen::Matrix<zsw::Scalar,3,3>> pplu;
    pplu.compute(A);
    Eigen::Matrix<zsw::Scalar,3,1> nv;
    for(size_t i=0; i<3; ++i) {
      if(vertices_[b_tr[i]].pt_type_==zsw::INNER_POINT) { nv[i]=-1.0-pt_val; }
      else if(vertices_[b_tr[i]].pt_type_==zsw::ZERO_POINT) { nv[i]=0.0-pt_val; }
      else { nv[i]=1.0-pt_val; }
    }

    // adjust jpts
    size_t judge_pts_cnt=0;
    for(JudgePoint &jpt : jpts) {
      //check if jpt in the tet and if jpt's error is in tolerance
      Eigen::Matrix<zsw::Scalar,3,1> ans = pplu.solve(jpt.pt_-pt);
      assert(((A*ans-(jpt.pt_-pt)).squaredNorm()<zsw::const_val::eps));

      if(ans.squaredNorm()>1+zsw::const_val::eps || ans[0]<-zsw::const_val::eps
         || ans[1]<-zsw::const_val::eps || ans[2]<-zsw::const_val::eps) { continue; } // not in tet
      jpt.val_cur_=pt_val+ans.dot(nv);
      tets_[*tet_itr].jpts_.push_back(jpt);
      ++judge_pts_cnt;
    }
    assert(judge_pts_cnt>=jpts.size());
    ++tet_itr;
  }

  // add new edges
  std::set<size_t> arround_v;
  for(const Eigen::Matrix<size_t,3,1> &b_tr : bound_tris) {
    arround_v.insert(b_tr[0]); arround_v.insert(b_tr[1]);
    arround_v.insert(b_tr[2]);
  }
  assert(inv_edge_ids.size() > arround_v.size());
  auto e_itr=inv_edge_ids.begin();
  for(const size_t v_id : arround_v) {
    REUSE_EDGE(e.vid_[0], v_id, *e_itr); ++e_itr;
  }

#if 0
  for(size_t ti=0; ti<tets_.size(); ++ti) {
    if(!tets_[ti].valid_) { continue; }
    for(size_t vi=0; vi<4; ++vi) {
      if(tets_[ti].vid_[vi]==e.vid_[0]) {
        std::cerr << "having e.vid[0] with :" << ti << std::endl;
        std::cerr << tets_[ti].vid_[0] << " " << tets_[ti].vid_[1]
                  << " " << tets_[ti].vid_[2] << " " << tets_[ti].vid_[3] << std::endl;
      }
    }
  }
#endif

}

void zsw::Triangulation::writeTet(const std::string &filepath, const size_t tet_id) const
{
  size_t tets_data[4] = {0,1,2,3};
  Eigen::Matrix<zsw::Scalar,3,4> pts_data;
  pts_data.block<3,1>(0,0)=vertices_[tets_[tet_id].vid_[0]].pt_;
  pts_data.block<3,1>(0,1)=vertices_[tets_[tet_id].vid_[1]].pt_;
  pts_data.block<3,1>(0,2)=vertices_[tets_[tet_id].vid_[2]].pt_;
  pts_data.block<3,1>(0,3)=vertices_[tets_[tet_id].vid_[3]].pt_;
  std::ofstream ofs;
  OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
  tet2vtk(ofs, pts_data.data(), 4, tets_data, 1);
}

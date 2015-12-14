#include "triangulation2.h"
#include <unordered_set>

#include <zswlib/mesh/vtk.h>
#include <zswlib/error_ctrl.h>
#include <zswlib/const_val.h>
#include <zswlib/zsw_log.h>

#include "sampling.h"
#include "basic_op.h"
#include "triangulation_fix.h"
#include "debug.h"
#include "config.h"

#define ADD_VERTEX(pt_type, pt) do{                                     \
    vertices_.push_back({true, pt_type, pt, Eigen::Matrix<zsw::Scalar,4,4>::Zero(), {}, {}}); \
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

#define FIND_ZERO_VERTEX_ID(va, vb, vid) do{            \
    auto it = ev_map.find(std::pair<size_t,size_t>(va,vb));     \
    ASSURE(it!=ev_map.end(), abort());                          \
    vid = it->second;                                           \
  }while(0)

#define REUSE_TET(v0, v1, v2, v3, t_id) do {                            \
    vertices_[v0].tet_ids_.push_back(t_id);   vertices_[v1].tet_ids_.push_back(t_id); \
    vertices_[v2].tet_ids_.push_back(t_id);   vertices_[v3].tet_ids_.push_back(t_id); \
    tets_[t_id] = {true, v0, v1, v2, v3};                          \
  } while(0)

#define CHECK_ADD_EDGE(v0, v1) do{                                      \
    bool need_add = true;                                               \
    for(size_t e_id : vertices_[v0].edge_ids_) {                        \
      if(edges_[e_id].vid_[0]==v1 || edges_[e_id].vid_[1]==v1) {        \
        need_add = false; break;                                        \
      }                                                                 \
    }                                                                   \
    if(need_add) { ADD_EDGE(v0, v1); }                                  \
  }while(0)

#define INTET_CONSTRAINT(int_judge, v0, v1, v2, v3) do{                 \
    int_judge.addConstraint(vertices_[v1].pt_, vertices_[v2].pt_, vertices_[v3].pt_, vertices_[v0].pt_); \
    int_judge.addConstraint(vertices_[v0].pt_, vertices_[v2].pt_, vertices_[v3].pt_, vertices_[v1].pt_); \
    int_judge.addConstraint(vertices_[v0].pt_, vertices_[v1].pt_, vertices_[v3].pt_, vertices_[v2].pt_); \
    int_judge.addConstraint(vertices_[v0].pt_, vertices_[v1].pt_, vertices_[v2].pt_, vertices_[v3].pt_); \
  }while(0)

#define UPDATE_MIN_JPT_DIS2FACE(vid0, vid1, vid2, v0, v1, v2) do{       \
    zsw::Scalar tmp_s_dis = calcPoint2TriSquaredDis(jpt_itr->pt_, v0, v1, v2); \
    if(squared_dis > tmp_s_dis) { target_face<<vid0, vid1, vid2; squared_dis=tmp_s_dis; target_bt_i=bt_i; } \
  }while(0)

size_t zsw::Triangulation::construct(const zsw::Scalar r, std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bo_pts,
                                  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bi_pts)
{
  assert(bo_pts.size()!=0 && bi_pts.size()!=0);
  tet_sample_r_ = r;

  Eigen::Matrix<zsw::Scalar,3,2> bbox;
  calcBBOX(bo_pts, bbox);
  zsw::Scalar scale = 0.5*(bbox.block<3,1>(0,0) - bbox.block<3,1>(0,1)).norm();
  Eigen::Matrix<zsw::Scalar,3,1> transform = 0.5*(bbox.block<3,1>(0,0)+bbox.block<3,1>(0,1));
  zsw::BoundSphere bs(BOUND_SPHERE_MODEL_FILE, scale, transform);

  size_t pt_id=0;
  std::vector<std::pair<Point, size_t>> tet_points;
  for(const Eigen::Matrix<zsw::Scalar,3,1> &tmp : bo_pts) {
    tet_points.push_back({Point(tmp[0], tmp[1], tmp[2]), pt_id++});
    vertices_.push_back({true, OUTER_POINT, tmp, Eigen::Matrix<zsw::Scalar,4,4>::Zero(), {}, {}});
  }

  for(const Eigen::Matrix<zsw::Scalar,3,1> &tmp : bi_pts) {
    tet_points.push_back({Point(tmp[0],tmp[1],tmp[2]), pt_id++});
    vertices_.push_back({true, INNER_POINT, tmp, Eigen::Matrix<zsw::Scalar,4,4>::Zero(), {}, {}});
  }

  // add bound sphere points
  const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bs_vertices = bs.getVertices();
  for(const Eigen::Matrix<zsw::Scalar,3,1> &v : bs_vertices) {
    tet_points.push_back({Point(v[0], v[1], v[2]), pt_id++});
    vertices_.push_back({true, BBOX_POINT, v, Eigen::Matrix<zsw::Scalar,4,4>::Zero(), {}, {}});
  }

  Delaunay delaunay(tet_points.begin(), tet_points.end());
  removeSliverTet(delaunay, vertices_);
  haveSliverTet(delaunay, vertices_);
  init(r, delaunay);
  CALL_FUNC(isGood(), return __LINE__);
  return 0;
}

void zsw::Triangulation::init(const zsw::Scalar r, Delaunay &delaunay)
{
  size_t t_id=0;
  for(Delaunay::Finite_cells_iterator cit=delaunay.finite_cells_begin();
      cit!=delaunay.finite_cells_end(); ++cit) {
    tets_.push_back({true, {cit->vertex(0)->info(), cit->vertex(1)->info(),
            cit->vertex(2)->info(), cit->vertex(3)->info()}});
    // judge points in tets
    size_t bo_cnt=0, bi_cnt=0;
    Eigen::Matrix<size_t,4,1> bo_vids;
    Eigen::Matrix<size_t,4,1> bi_vids;
    Eigen::Matrix<zsw::Scalar,3,4> bo_points;
    Eigen::Matrix<zsw::Scalar,3,4> bi_points;
    Tet &tmp_tet=tets_.back();
    for(size_t tvid : tmp_tet.vid_) {
      if(vertices_[tvid].pt_type_==OUTER_POINT) {
        bo_vids[bo_cnt] = tvid;
        bo_points.block<3,1>(0,bo_cnt++)=vertices_[tvid].pt_;
      } else if(vertices_[tvid].pt_type_==INNER_POINT) {
        bi_vids[bi_cnt] = tvid;
        bi_points.block<3,1>(0,bi_cnt++)=vertices_[tvid].pt_;
      } else {        bo_cnt=0; bi_cnt=0;        break;      }
    }

    if(bo_cnt==3) { sampleTriangle(bo_points.block<3,3>(0,0), r, 1, tmp_tet.jpts_); initQem(bo_vids[0], bo_vids[1], bo_vids[2]); }
    else if(bi_cnt==3) { sampleTriangle(bi_points.block<3,3>(0,0), r, -1, tmp_tet.jpts_); initQem(bi_vids[0], bi_vids[1], bi_vids[2]); }
    // vertex info
    vertices_[tmp_tet.vid_[0]].tet_ids_.push_back(t_id);  vertices_[tmp_tet.vid_[1]].tet_ids_.push_back(t_id);
    vertices_[tmp_tet.vid_[2]].tet_ids_.push_back(t_id);  vertices_[tmp_tet.vid_[3]].tet_ids_.push_back(t_id++);
  }

  // jpts_ptr_.reset(new zsw::Flann2<zsw::Scalar, 2>(all_jpts_[0].pt_.data(), all_jpts_.size()));

  // edge info
  size_t e_id=0;
  for(Delaunay::Finite_edges_iterator eit=delaunay.finite_edges_begin();
      eit!=delaunay.finite_edges_end(); ++eit) {
    edges_.push_back({true, {eit->first->vertex(eit->second)->info(), eit->first->vertex(eit->third)->info()}});
    Edge &tmp_edge=edges_.back();
    vertices_[tmp_edge.vid_[0]].edge_ids_.push_back(e_id);
    vertices_[tmp_edge.vid_[1]].edge_ids_.push_back(e_id++);
  }
  std::cerr << "tol edge size:" << e_id << std::endl;
}

void zsw::Triangulation::initQem(const size_t vid0, const size_t vid1, const size_t vid2)
{
  Eigen::Matrix<zsw::Scalar,3,1> v0 = vertices_[vid0].pt_;
  Eigen::Matrix<zsw::Scalar,3,1> v1 = vertices_[vid1].pt_;
  Eigen::Matrix<zsw::Scalar,3,1> v2 = vertices_[vid2].pt_;
  Eigen::Matrix<zsw::Scalar,3,1> vn = (v1-v0).cross(v2-v0); vn.normalize();
  Eigen::Matrix<zsw::Scalar,4,1> p;
  p.block<3,1>(0,0) = vn; p[3] = -vn.dot(v0);
  Eigen::Matrix<zsw::Scalar, 4,4> q = p * p.transpose();
  vertices_[vid0].qem_+=q;
  vertices_[vid1].qem_+=q;
  vertices_[vid2].qem_+=q;
}

bool zsw::Triangulation::linkCondition(const Edge &e) const
{
  std::set<size_t> adj_v[2]; // vertex link e.vid_[0] and e.vid_[1]
  for(size_t vind=0; vind<2; ++vind) {
    for(size_t t_id : vertices_[e.vid_[vind]].tet_ids_) {
      assert(tets_[t_id].valid_);
      adj_v[vind].insert(tets_[t_id].vid_[0]); adj_v[vind].insert(tets_[t_id].vid_[1]);
      adj_v[vind].insert(tets_[t_id].vid_[2]); adj_v[vind].insert(tets_[t_id].vid_[3]);
    }
  }
  size_t cv_cnt=0;
  {
    std::set<size_t>::iterator it0=adj_v[0].begin();
    std::set<size_t>::iterator it1=adj_v[1].begin();
    while(it0!=adj_v[0].end() && it1!=adj_v[1].end()) {
      if(*it0 == *it1) { ++cv_cnt; ++it0; ++it1; }
      else if(*it0>*it1) { ++it1; }
      else { ++it0; }
    }
    if(adj_v[0].find(e.vid_[0])!=adj_v[0].end() && adj_v[1].find(e.vid_[0])!=adj_v[1].end()) { --cv_cnt; }
    if(adj_v[0].find(e.vid_[1])!=adj_v[0].end() && adj_v[1].find(e.vid_[1])!=adj_v[1].end()) { --cv_cnt; }
  }

  std::unordered_set<size_t> fv; // vertex construct a face with edge e
  for(size_t t_id : vertices_[e.vid_[0]].tet_ids_) {
    for(size_t vid : tets_[t_id].vid_) {
      if(vid == e.vid_[1]) {
        fv.insert(tets_[t_id].vid_[0]); fv.insert(tets_[t_id].vid_[1]);
        fv.insert(tets_[t_id].vid_[2]); fv.insert(tets_[t_id].vid_[3]);
        break;
      }
    }
  }
  size_t fv_cnt=fv.size();
  if(fv.find(e.vid_[0]) != fv.end()) { --fv_cnt; }
  if(fv.find(e.vid_[1]) != fv.end()) { --fv_cnt; }

  assert(cv_cnt >=fv_cnt);
  return cv_cnt==fv_cnt;
}

void zsw::Triangulation::simpTolerance()
{
  std::queue<size_t> eids;
  std::set<size_t> eids_set;
  for(size_t i=0; i<edges_.size(); ++i) {
    const Edge &e=edges_[i];
    if(!e.valid_ || vertices_[e.vid_[0]].pt_type_!=vertices_[e.vid_[1]].pt_type_ ||
       (vertices_[e.vid_[0]].pt_type_!=OUTER_POINT && vertices_[e.vid_[0]].pt_type_!=INNER_POINT) )
      { continue; }
    eids.push(i); eids_set.insert(i);
  }
  while(!eids.empty()) {
    size_t cur_eid=eids.front(); eids_set.erase(cur_eid); eids.pop();
    if(!edges_[cur_eid].valid_) { continue; }
    tryCollapseBoundaryEdge(cur_eid, eids, eids_set);
  }
}

void zsw::Triangulation::testCollapseDebug(const size_t vid0, const size_t vid1)
{
  size_t e_id = -1;
  for(size_t i=0; i<edges_.size(); ++i) {
    size_t vcnt=0;
    if(edges_[i].vid_[0]==vid0 || edges_[i].vid_[1]==vid0) ++vcnt;
    if(edges_[i].vid_[0]==vid1 || edges_[i].vid_[1]==vid1) ++vcnt;
    if(vcnt == 2) { e_id=i; break; }
  }

  std::queue<size_t> eids;
  std::set<size_t> eids_set;
  tryCollapseBoundaryEdge(e_id, eids, eids_set);
}

void zsw::Triangulation::simpZeroSurface()
{
  std::queue<size_t> eids;
  std::set<size_t> eids_set;
  for(size_t i=0; i<edges_.size(); ++i) {
    const Edge &e=edges_[i];
    if(!e.valid_ || vertices_[e.vid_[0]].pt_type_!=vertices_[e.vid_[1]].pt_type_ ||
       (vertices_[e.vid_[0]].pt_type_!=ZERO_POINT && vertices_[e.vid_[0]].pt_type_!=ZERO_POINT) )
      { continue; }
    eids.push(i); eids_set.insert(i);
  }

  while(!eids.empty()) {
    size_t cur_eid = eids.front(); eids_set.erase(cur_eid); eids.pop();
    if(!edges_[cur_eid].valid_) { continue; }
    tryCollapseZeroEdge(cur_eid, eids, eids_set);
  }
  std::cout << __FILE__ << __LINE__ << std::endl;
}

void zsw::Triangulation::addZeroPoints(std::map<std::pair<size_t,size_t>, size_t, PairCompFunc> &ev_map)
{
  const size_t edges_size = edges_.size();
  for(size_t e_id=0; e_id<edges_size; ++e_id) {
    Edge &e = edges_[e_id];
    if(!e.valid_) { continue; }
    if(e.valid_ && vertices_[e.vid_[0]].pt_type_==OUTER_POINT && vertices_[e.vid_[1]].pt_type_==INNER_POINT ||
       vertices_[e.vid_[0]].pt_type_==INNER_POINT && vertices_[e.vid_[1]].pt_type_==OUTER_POINT) {
      size_t vid = vertices_.size();
      ADD_VERTEX(ZERO_POINT, (vertices_[e.vid_[0]].pt_+vertices_[e.vid_[1]].pt_)/2.0);
      invalidEdge(e_id);
      ADD_EDGE(vid, e.vid_[0]); ADD_EDGE(vid, e.vid_[1]);
      if(e.vid_[0] > e.vid_[1]) { std::swap(e.vid_[0], e.vid_[1]); }
      std::pair<size_t,size_t> ep = std::pair<size_t,size_t>(e.vid_[0], e.vid_[1]);
      ev_map[ep]=vid;
    }
  }
}

bool pairComp(const std::pair<size_t, size_t> &a, const std::pair<size_t, size_t> &b)
{
  if(a.first < b.first) { return true; }
  if(a.first==b.first && a.second<b.second) { return true; }
  return false;
}

void zsw::Triangulation::mutualTessellation()
{
  // add zero point with edge
  std::map<std::pair<size_t,size_t>, size_t, PairCompFunc> ev_map(pairComp); // bi bo edge to vertex index map
  addZeroPoints(ev_map);
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
      //if(isGood()!=0) { std::cout << ""_FILE__ << __LINE__ << std::endl; abort(); }
    }
    // bo : bi = 2 : 2
    else if(vo_cnt==2 && vi_cnt==2) {
      tessellation2v2(vo[0],vo[1],vi[0],vi[1],tet, ev_map);
      //if(isGood()!=0) { std::cout << __FILE__ << __LINE__ << std::endl; abort(); }
    }
    // bo : bi = 1 : 3
    else if(vo_cnt==1 && vi_cnt==3) {
      tessellation1v3(vo[0],vi[0],vi[1],vi[2],tet, ev_map);
      //if(isGood()!=0) { std::cout << __FILE__ << __LINE__ << std::endl; abort(); }
    }
  }

#ifdef ZSW_DEBUG
  size_t ntet_size = tets_.size();
  std::cerr << "size_old:" << tet_size << std::endl;
  for(size_t t_id=0; t_id<ntet_size; ++t_id) {
    Tet &tet=tets_[t_id];
    if(!tet.valid_) { continue; }
    size_t vo_cnt=0, vi_cnt=0;
    for(size_t vid : tet.vid_) {
      if(vertices_[vid].pt_type_ == OUTER_POINT) { vo_cnt++; }
      else if(vertices_[vid].pt_type_ == INNER_POINT){ vi_cnt++;  }
    }
    if(vo_cnt!=0 && vi_cnt!=0) {
      std::cerr << "error " << t_id << std::endl;
      abort();
    }
  }
#endif
}

void zsw::Triangulation::tessellation3v1(const size_t vo_0,
                                         const size_t vo_1, const size_t vo_2, const size_t vi_0,
                                         Tet &tet, std::map<std::pair<size_t,size_t>, size_t, PairCompFunc> &ev_map)
{
  // find the zero point
  size_t nv0, nv1, nv2;
  FIND_ZERO_VERTEX_ID(vo_0,vi_0,nv0);
  FIND_ZERO_VERTEX_ID(vo_1,vi_0,nv1);
  FIND_ZERO_VERTEX_ID(vo_2,vi_0,nv2);

  // invalid old tet and edges
  invalidTet(tet);
  // add tets
  Eigen::Matrix<size_t,4,1> tets_info[4];
  std::vector<size_t> tids;
  tids.push_back(tets_.size()); tids.push_back(tets_.size()+1);
  tids.push_back(tets_.size()+2);  tids.push_back(tets_.size()+3);
  tets_info[0] << vi_0, nv0, nv1, nv2;
  tets_info[1] << nv0, vo_0, vo_2, vo_1;
  tets_info[2] << nv1, nv0, vo_2, vo_1;
  tets_info[3] << nv2, nv1, nv0, vo_2;
  for(size_t ti=0; ti<4; ++ti) {
    ADD_TET(tets_info[ti][0], tets_info[ti][1], tets_info[ti][2], tets_info[ti][3]);
  }
  splitJudgePoints(tet.jpts_, tids);

  // add edges
  CHECK_ADD_EDGE(nv0,nv1); CHECK_ADD_EDGE(nv0, nv2);
  CHECK_ADD_EDGE(nv1,nv2); CHECK_ADD_EDGE(nv0, vo_1);
  CHECK_ADD_EDGE(nv0, vo_2); CHECK_ADD_EDGE(nv1,vo_2);

#ifdef ZSW_DEBUG
  checkTetEdgeExist(vi_0, nv0, nv1, nv2); checkTetEdgeExist(nv0, vo_0, vo_2, vo_1);
  checkTetEdgeExist(nv1, nv0, vo_2, vo_1); checkTetEdgeExist(nv2, nv1, nv0, vo_2);
#endif
}

void zsw::Triangulation::tessellation2v2(const size_t vo_0, const size_t vo_1,
                                         const size_t vi_0, const size_t vi_1,
                                         Tet &tet,
                                         std::map<std::pair<size_t,size_t>, size_t, PairCompFunc> &ev_map)
{
  size_t nv0, nv1, nv2, nv3;
  FIND_ZERO_VERTEX_ID(vo_0,vi_0,nv0);
  FIND_ZERO_VERTEX_ID(vo_0, vi_1,nv1);
  FIND_ZERO_VERTEX_ID(vo_1,vi_0,nv2);
  FIND_ZERO_VERTEX_ID(vo_1,vi_1,nv3);

  // invalid old tet and  add new tet
  invalidTet(tet);
  Eigen::Matrix<size_t,4,1> tets_info[6];
  std::vector<size_t> tids;
  tids.push_back(tets_.size());   tids.push_back(tets_.size()+1);
  tids.push_back(tets_.size()+2);   tids.push_back(tets_.size()+3);
  tids.push_back(tets_.size()+4);   tids.push_back(tets_.size()+5);
  tets_info[0] <<nv0, nv1, nv2, vo_0;
  tets_info[1] <<nv0, nv1, nv2, vi_0;
  tets_info[2] <<nv1, nv2, vi_0, vi_1;
  tets_info[3] <<nv1, nv2, nv3, vi_1;
  tets_info[4] <<nv2, nv1, vo_0, vo_1;
  tets_info[5] <<nv2, nv1, nv3, vo_1;
  for(size_t ti=0; ti<6; ++ti) {
    ADD_TET(tets_info[ti][0], tets_info[ti][1], tets_info[ti][2], tets_info[ti][3]);
  }
  splitJudgePoints(tet.jpts_, tids);
  // add edges
  CHECK_ADD_EDGE(nv0,nv2);  CHECK_ADD_EDGE(nv0,nv1);
  CHECK_ADD_EDGE(nv1,vi_0);  CHECK_ADD_EDGE(nv1,nv3);
  CHECK_ADD_EDGE(nv1,nv2);  CHECK_ADD_EDGE(nv1,vo_1);
  CHECK_ADD_EDGE(nv2,nv3); CHECK_ADD_EDGE(nv2,vo_0);
  CHECK_ADD_EDGE(nv2,vi_1);

#ifdef ZSW_DEBUG
  checkTetEdgeExist(nv0, nv1, nv2, vo_0);  checkTetEdgeExist(nv0, nv1, nv2, vi_0);
  checkTetEdgeExist(nv1, nv2, vi_0, vi_1);  checkTetEdgeExist(nv1, nv2, nv3, vi_1);
  checkTetEdgeExist(nv2, nv1, vo_0, vo_1);  checkTetEdgeExist(nv2, nv1, nv3, vo_1);
#endif
}

void zsw::Triangulation::tessellation1v3(const size_t vo_0, const size_t vi_0,
                                         const size_t vi_1, const size_t vi_2,
                                         Tet &tet,
                                         std::map<std::pair<size_t,size_t>, size_t, PairCompFunc> &ev_map)
{
  // find zero points
  size_t nv0, nv1, nv2;
  FIND_ZERO_VERTEX_ID(vo_0,vi_0,nv0);
  FIND_ZERO_VERTEX_ID(vo_0,vi_1,nv1);
  FIND_ZERO_VERTEX_ID(vo_0,vi_2,nv2);

  // invalid old tet and add new tets
  invalidTet(tet);
  Eigen::Matrix<size_t,4,1> tets_info[4];
  std::vector<size_t> tids;
  tids.push_back(tets_.size());  tids.push_back(tets_.size()+1);
  tids.push_back(tets_.size()+2);  tids.push_back(tets_.size()+3);
  tets_info[0]<<vo_0, nv0, nv1, nv2;
  tets_info[1]<<nv0, vi_0, vi_1, vi_2;
  tets_info[2]<<nv1, nv0, vi_1, vi_2;
  tets_info[3]<<nv2, nv0, nv1, vi_2;
  for(size_t ti=0; ti<4; ++ti) {
    ADD_TET(tets_info[ti][0], tets_info[ti][1], tets_info[ti][2], tets_info[ti][3]);
  }
  splitJudgePoints(tet.jpts_, tids);
  // ADD_TET(vo_0, nv0, nv1, nv2);  ADD_TET(nv0, vi_0, vi_1, vi_2);
  // ADD_TET(nv1, nv0, vi_1, vi_2);  ADD_TET(nv2, nv0, nv1, vi_2);

  // add edge
  CHECK_ADD_EDGE(nv0,nv1);  CHECK_ADD_EDGE(nv0,vi_1);
  CHECK_ADD_EDGE(nv0,vi_2);  CHECK_ADD_EDGE(nv1,vi_2);
  CHECK_ADD_EDGE(nv1,nv2);  CHECK_ADD_EDGE(nv0,nv2);

#ifdef ZSW_DEBUG
  checkTetEdgeExist(vo_0, nv0, nv1, nv2);  checkTetEdgeExist(nv0, vi_0, vi_1, vi_2);
  checkTetEdgeExist(nv1, nv0, vi_1, vi_2);  checkTetEdgeExist(nv2, nv0, nv1, vi_2);
#endif
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

void zsw::Triangulation::writeTetMesh(const std::string &filepath,
                                      std::vector<std::function<bool(const Tet &tet)>> ignore_tet_funcs) const
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
    for(std::function<bool(const Tet &tet)> &ig_func : ignore_tet_funcs) {
      if(ig_func(tet)) { ignore=true; break; }
    }
    if(ignore) { continue; }
    tets_data.push_back(tet.vid_[0]);
    tets_data.push_back(tet.vid_[1]);
    tets_data.push_back(tet.vid_[2]);
    tets_data.push_back(tet.vid_[3]);
  }
  size_t n_pts = pts_data.size()/3;
  size_t n_tets = tets_data.size()/4;
  NZSWLOG("zsw_info")  << "point size:" << n_pts << std::endl;
  NZSWLOG("zsw_info")  << "tet size:" << n_tets << std::endl;
  tet2vtk(ofs, &pts_data[0], n_pts, &tets_data[0], n_tets);
}

void zsw::Triangulation::writeSurface(const std::string &filepath, PointType pt_type) const
{
  std::ofstream ofs(filepath);
  std::set<size_t> valid_vid;
  std::vector<Eigen::Matrix<size_t,3,1>> faces;
  Eigen::Matrix<size_t,4,1> zv_id;
  for(const Tet &tet : tets_) {
    if(!tet.valid_) { continue; }
    size_t id_cnt=0;
    for(size_t v_id : tet.vid_) {
      if(vertices_[v_id].pt_type_ == pt_type) {
        zv_id[id_cnt++]=v_id;
      }
    }
    if(id_cnt==3) {
      valid_vid.insert(zv_id[0]); valid_vid.insert(zv_id[1]); valid_vid.insert(zv_id[2]);
      faces.push_back(zv_id.block<3,1>(0,0));
    }
  }
  std::map<size_t,size_t> vv_map;
  size_t vid_cnt=0;
  for(size_t v_id : valid_vid) {
    vv_map[v_id]=++vid_cnt;
    ofs << "v " << vertices_[v_id].pt_.transpose() << std::endl;
  }
  for(const Eigen::Matrix<size_t,3,1> &f : faces) {
    ofs << "f " << vv_map[f[0]] << " " << vv_map[f[1]] << " " << vv_map[f[2]] << std::endl;
  }
}

bool zsw::Triangulation::isKeepJpts(const zsw::Scalar pt_val, const Eigen::Matrix<zsw::Scalar,3,1> &pt,
                                    const std::list<Eigen::Matrix<size_t,3,1>> &bound_tris, const std::list<JudgePoint> &all_jpts,
                                    std::vector<std::pair<size_t, zsw::Scalar>> &jpts_update) const
{
  int bt_i=-1;
  size_t update_cnt=0;
  jpts_update.assign(all_jpts.size(), std::pair<size_t,zsw::Scalar>(-1, 0.0));
  for(const Eigen::Matrix<size_t,3,1> &b_tr  : bound_tris) {
    ++bt_i;
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

    auto up_itr = jpts_update.begin();
    auto jpt_itr = all_jpts.begin();
    for(; up_itr!=jpts_update.end(); ++up_itr, ++jpt_itr) {
      if(up_itr->first != -1) { continue; } // allready updated
      const JudgePoint &jpt = *jpt_itr;
      Eigen::Matrix<zsw::Scalar,3,1> ans = pplu.solve(jpt.pt_-pt);
      assert((A*ans-(jpt.pt_-pt)).norm()<zsw::const_val::eps);
      if((A*ans-(jpt.pt_-pt)).norm()>zsw::const_val::eps) { continue; }
      if(ans[0]<0 || ans[1]<0 || ans[2]<0 || ans[0]+ans[1]+ans[2]>1) { continue; } // not in tet
      zsw::Scalar jpt_val_cur=pt_val+ans.dot(nv);
      if(fabs(jpt_val_cur-jpt.val_exp_) > 1.0+zsw::const_val::eps) { return false; }
      else {  up_itr->first=bt_i; up_itr->second=jpt_val_cur;  ++update_cnt; }
    }
  }
  //std::cerr << "left jpts size:" << all_jpts.size()-update_cnt << std::endl;
  return (update_cnt==all_jpts.size()) || isKeepJptsLeft(pt_val, pt, bound_tris, all_jpts, jpts_update);
}

bool zsw::Triangulation::isKeepJptsLeft(const zsw::Scalar pt_val, const Eigen::Matrix<zsw::Scalar,3,1> &pt,
                                        const std::list<Eigen::Matrix<size_t,3,1>> &bound_tris,
                                        const std::list<JudgePoint> &all_jpts,
                                        std::vector<std::pair<size_t, zsw::Scalar>> &jpts_update) const
{
  auto up_itr = jpts_update.begin();
  auto jpt_itr = all_jpts.begin();
  for(; up_itr!=jpts_update.end(); ++up_itr, ++jpt_itr) {
    if(up_itr->first != -1) { continue; } // already have update info
    // find the nearest triangle -> target face, and its tet
    size_t bt_i=0;
    size_t target_bt_i=0;
    zsw::Scalar squared_dis = (jpt_itr->pt_-pt).squaredNorm();
    Eigen::Matrix<size_t,3,1> target_face = bound_tris.front();
    for(const Eigen::Matrix<size_t,3,1> &b_tr : bound_tris) {
      UPDATE_MIN_JPT_DIS2FACE(b_tr[0], b_tr[1], b_tr[2], vertices_[b_tr[0]].pt_, vertices_[b_tr[1]].pt_, vertices_[b_tr[2]].pt_);
      UPDATE_MIN_JPT_DIS2FACE(-1, b_tr[0], b_tr[1], pt, vertices_[b_tr[0]].pt_, vertices_[b_tr[1]].pt_);
      UPDATE_MIN_JPT_DIS2FACE(-1, b_tr[0], b_tr[2], pt, vertices_[b_tr[0]].pt_, vertices_[b_tr[2]].pt_);
      UPDATE_MIN_JPT_DIS2FACE(-1, b_tr[1], b_tr[2], pt, vertices_[b_tr[1]].pt_, vertices_[b_tr[2]].pt_);
      ++bt_i;
    }

    zsw::Scalar val0;
    Eigen::Matrix<zsw::Scalar,3,1> tri_v[3];
    if(target_face[0] == -1) { tri_v[0]=pt; val0=pt_val; }
    else {
      tri_v[0]=vertices_[target_face[0]].pt_;
      if(vertices_[target_face[0]].pt_type_==zsw::INNER_POINT) { val0=-1.0; }
      else if(vertices_[target_face[0]].pt_type_==zsw::ZERO_POINT) { val0=0.0; }
      else { val0=1.0; }
    }
    tri_v[1] = vertices_[target_face[1]].pt_;    tri_v[2] = vertices_[target_face[2]].pt_;
    Eigen::Matrix<zsw::Scalar,2,1> nv;
    if(vertices_[target_face[1]].pt_type_==zsw::INNER_POINT) { nv[0]=-1.0-val0; }
    else if(vertices_[target_face[1]].pt_type_==zsw::ZERO_POINT) { nv[0]=-val0; }
    else { nv[0]=1.0-val0; }
    if(vertices_[target_face[2]].pt_type_==zsw::INNER_POINT) { nv[1]=-1.0-val0; }
    else if(vertices_[target_face[2]].pt_type_==zsw::ZERO_POINT) { nv[1]=-val0; }
    else { nv[1]=1.0-val0; }

    // calc jpt's update data
    Eigen::Matrix<zsw::Scalar,3,2> A;
    A.block<3,1>(0,0)=tri_v[1] - tri_v[0];
    A.block<3,1>(0,1)=tri_v[2] - tri_v[0];
    Eigen::Matrix<zsw::Scalar,2,3> AT = A.transpose();
    Eigen::PartialPivLU<Eigen::Matrix<zsw::Scalar,2,2>> pplu;
    pplu.compute(AT*A);
    Eigen::Matrix<zsw::Scalar,2,1> ans = pplu.solve(AT*(jpt_itr->pt_ - tri_v[0]));
    zsw::Scalar jpt_val_cur=ans.dot(nv)+val0;
    if(fabs(ans.dot(nv)+val0-jpt_itr->val_exp_) > 1+zsw::const_val::eps) { return false; }
    else { up_itr->first=target_bt_i; up_itr->second=jpt_val_cur;    }
  }
  return true;
}

void zsw::Triangulation::edgeCollapse(const std::vector<size_t> &tet_ids,
                                      const std::list<Eigen::Matrix<size_t,3,1>> &bound_tris,
                                      const Eigen::Matrix<zsw::Scalar,3,1> &pt,
                                      const zsw::PointType pt_type,
                                      const std::vector<std::pair<size_t, zsw::Scalar>> &jpts_update,
                                      Edge &e,
                                      std::list<JudgePoint> &all_jpts,
                                      std::function<void(const size_t e_id)> eb_func)
{
  // invalid old tets and edges and vertex
  std::unordered_set<size_t> inv_edge_ids;
  vertices_[e.vid_[0]].valid_=false;  vertices_[e.vid_[1]].valid_=false;
  for(size_t e_id : vertices_[e.vid_[0]].edge_ids_) { inv_edge_ids.insert(e_id); }
  for(size_t e_id : vertices_[e.vid_[1]].edge_ids_) { inv_edge_ids.insert(e_id); }

  for(size_t t_id : tet_ids) { invalidTet(tets_[t_id]); }
  for(size_t e_id : inv_edge_ids) { invalidEdge(e_id); }

  // new vertex
  REUSE_VERTEX(pt, pt_type, e.vid_[0]);
  // add new tets
  assert(tet_ids.size()>=bound_tris.size());
  auto tet_itr = tet_ids.begin();
  std::vector<size_t> reuse_tet_ids;
  for(const Eigen::Matrix<size_t,3,1> &bt_r : bound_tris) {
    REUSE_TET(e.vid_[0], bt_r[0], bt_r[1], bt_r[2], *tet_itr);
    reuse_tet_ids.push_back(*tet_itr); ++tet_itr;
  }

  assert(all_jpts.size() == jpts_update.size());
  auto jpt_itr=all_jpts.begin();
  for(const std::pair<size_t,zsw::Scalar> &jpt_update : jpts_update) {
    assert(jpt_update.first < reuse_tet_ids.size());
    std::list<JudgePoint> &cur_tet_jpts = tets_[reuse_tet_ids[jpt_update.first]].jpts_;
    jpt_itr->val_cur_=jpt_update.second;
    auto tmp_jpt_itr = jpt_itr; ++jpt_itr;
    cur_tet_jpts.splice(cur_tet_jpts.end(), all_jpts, tmp_jpt_itr);
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
    const size_t tmp_eid=*e_itr;
    REUSE_EDGE(e.vid_[0], v_id, tmp_eid); ++e_itr;
    eb_func(tmp_eid);
  }
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

  bool zsw::Triangulation::ignoreWithPtType(const Tet &tet, PointType pt_type)
  {
    for(size_t i=0; i<4; ++i) {
      if(vertices_[tet.vid_[i]].pt_type_ == pt_type) {
        return true;
      }
    }
    return false;
  }

  bool zsw::Triangulation::ignoreOnlyWithPtType(const Tet &tet, PointType pt_type)
  {
    return vertices_[tet.vid_[0]].pt_type_==pt_type && vertices_[tet.vid_[1]].pt_type_==pt_type &&
      vertices_[tet.vid_[2]].pt_type_==pt_type && vertices_[tet.vid_[3]].pt_type_==pt_type;
  }

  bool zsw::Triangulation::ignoreNotWithPtType(const Tet &tet, PointType pt_type)
  {
    return !(vertices_[tet.vid_[0]].pt_type_==pt_type || vertices_[tet.vid_[1]].pt_type_==pt_type ||
             vertices_[tet.vid_[2]].pt_type_==pt_type || vertices_[tet.vid_[3]].pt_type_==pt_type);
  }


void zsw::Triangulation::writeTetMeshAdjVs(const std::string &filepath, const std::vector<size_t> &vids) const
{
  std::set<size_t> tet_ids;
  for(size_t vid : vids) {
    for_each(vertices_[vid].tet_ids_.begin(), vertices_[vid].tet_ids_.end(),
             [&tet_ids](const size_t t_id) { tet_ids.insert(t_id); });
  }

  std::ofstream ofs;
  OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
  std::vector<zsw::Scalar> pts_data;
  std::vector<size_t> tets_data;
  for(const Vertex &v : vertices_) {
    pts_data.push_back(v.pt_[0]);
    pts_data.push_back(v.pt_[1]);
    pts_data.push_back(v.pt_[2]);
  }
  size_t n_tets = 0;
  for(size_t tid : tet_ids) {
    if(!tets_[tid].valid_) { continue; }
    ++n_tets;
    tets_data.push_back(tets_[tid].vid_[0]);
    tets_data.push_back(tets_[tid].vid_[1]);
    tets_data.push_back(tets_[tid].vid_[2]);
    tets_data.push_back(tets_[tid].vid_[3]);
  }
  tet2vtk(ofs, &pts_data[0], vertices_.size(), &tets_data[0], n_tets);
}

void zsw::Triangulation::tryCollapseBoundaryEdge(const size_t e_id,
                                                 std::queue<size_t> &eids,
                                                 std::set<size_t> eids_set)
{
  Edge &e = edges_[e_id];
  const zsw::Scalar pt_val = (vertices_[e.vid_[0]].pt_type_==zsw::INNER_POINT) ? -1 : 1;
  // if e is invalid or is not bi and bo edge, if is bi bo edge is checked when push into the queue
  if(!e.valid_)   { return; }
  if(!linkCondition(e)) { return; }
  std::vector<size_t> tet_ids;
  for(size_t tid : vertices_[e.vid_[0]].tet_ids_) { tet_ids.push_back(tid); }
  for(size_t tid : vertices_[e.vid_[1]].tet_ids_) { tet_ids.push_back(tid); }
  std::sort(tet_ids.begin(), tet_ids.end());
  auto unq_end = std::unique(tet_ids.begin(), tet_ids.end());
  tet_ids.resize(std::distance(tet_ids.begin(), unq_end));
  // normal condition
  NormalConditionJudger ncj(NORMAL_CONT_TOL);
  initNormalCond(ncj, e);
  KernelRegionJudger krj;
  std::list<Eigen::Matrix<size_t,3,1>> bound_tris;
  std::list<JudgePoint> all_jpts;
  for(size_t tid : tet_ids) {
    all_jpts.insert(all_jpts.end(), tets_[tid].jpts_.begin(), tets_[tid].jpts_.end());
    // add kernel region constraint
    size_t vcnt=0;
    Eigen::Matrix<size_t,4,1> tmp_bound_tri;
    for(size_t tvid : tets_[tid].vid_) {
      if(tvid == e.vid_[0] || tvid==e.vid_[1]) {          tmp_bound_tri[3]=tvid;        }
      else {          tmp_bound_tri[vcnt++]=tvid;        }
    }
    if(vcnt==3) {
      krj.addConstraint(vertices_[tmp_bound_tri[0]].pt_, vertices_[tmp_bound_tri[1]].pt_,
                        vertices_[tmp_bound_tri[2]].pt_, vertices_[tmp_bound_tri[3]].pt_);
      bound_tris.push_back(tmp_bound_tri.block<3,1>(0,0));
    }
  }
  // find the candicate points :  jpt in kernel region
  std::vector<JudgePoint> candicate_pts;
  for(const JudgePoint &jpt : all_jpts) {
    if(fabs(jpt.val_exp_-pt_val)<0.5 && krj.judge(jpt.pt_) && ncj.judge(jpt.pt_)) {
      candicate_pts.push_back(jpt);
    }
  }
  static size_t step_info = 0;
  ++step_info;
  if(step_info %100 == 0) {
    NZSWLOG("zsw_info") << "edge:" << e.vid_[0] << " : "  << e.vid_[1] << std::endl;
    NZSWLOG("zsw_info") << "candicate_pts:" << candicate_pts.size() << std::endl;
    NZSWLOG("zsw_info") << "all_pts:" << all_jpts.size() << std::endl;
  }
  // debug
#if 0
  writeJudgePoints("/home/wegatron/tmp/simp_tol/debug/candicate_jpts", candicate_pts);
  writeJudgePoints("/home/wegatron/tmp/simp_tol/debug/jpts", all_jpts);
#endif
  //-- find the best point in the candicate_pts--
#if 0
  // sort by fabs(val_exp-val_cur) by decrease order
  std::sort(candicate_pts.begin(), candicate_pts.end(),
            [](const JudgePoint &a, const JudgePoint &b){ return fabs(a.val_cur_-a.val_exp_)>fabs(b.val_cur_-b.val_exp_); });
  const JudgePoint *merge_pt_ptr=nullptr;
  std::vector<std::pair<size_t,zsw::Scalar>> jpts_update;
  int step_info=0;
  for(const JudgePoint &jpt : candicate_pts) {
    if((++step_info)%100==0) { NZSWLOG("zsw_info") << step_info << " - "; }
    if(isKeepJpts(pt_val, jpt.pt_, bound_tris, all_jpts, jpts_update)) {
      merge_pt_ptr= &jpt; break;
    }
  }
  if(merge_pt_ptr != nullptr) {
    if(step_info %100 == 0) {
      NZSWLOG("zsw_info")  << "bc: collapse edge!!!" << std::endl;
    }
    edgeCollapse(tet_ids, bound_tris, merge_pt_ptr->pt_, vertices_[e.vid_[0]].pt_type_, jpts_update, e, all_jpts,
                 [&eids, &eids_set,this](const size_t e_id){
                   if(vertices_[edges_[e_id].vid_[0]].pt_type_==vertices_[edges_[e_id].vid_[1]].pt_type_
                      && (vertices_[edges_[e_id].vid_[0]].pt_type_==OUTER_POINT || vertices_[edges_[e_id].vid_[0]].pt_type_==INNER_POINT)
                      && eids_set.find(e_id)!=eids_set.end()) {
                     eids.push(e_id); eids_set.insert(e_id); }
                 });
  } else if(step_info %100 == 0) {      NZSWLOG("zsw_info")  << "bc: no poper merge point!"<< std::endl;    }
#else
  // using qem in accending order
  Eigen::Matrix<zsw::Scalar,4,4> cur_qem = vertices_[e.vid_[0]].qem_ + vertices_[e.vid_[1]].qem_;
  std::vector<std::pair<zsw::Scalar, size_t>> q_error; q_error.resize(candicate_pts.size());
  for(size_t q_i=0; q_i<q_error.size(); ++q_i) {
    Eigen::Matrix<zsw::Scalar,4,1> tmp_pt;
    tmp_pt.block<3,1>(0,0) = candicate_pts[q_i].pt_; tmp_pt[3]=1.0;
    q_error[q_i].first = tmp_pt.transpose() * cur_qem * tmp_pt;
    q_error[q_i].second = q_i;
  }
  std::sort(q_error.begin(), q_error.end(),
            [](const std::pair<zsw::Scalar, size_t> &a, const std::pair<zsw::Scalar, size_t> &b){ return a.first<b.first; });
  const JudgePoint *merge_pt_ptr=nullptr;
  std::vector<std::pair<size_t,zsw::Scalar>> jpts_update;
  //int step_info=0;
  for(const std::pair<zsw::Scalar, size_t> &qe : q_error) {
    const JudgePoint &jpt = candicate_pts[qe.second];
    //    if((++step_info)%100==0) { NZSWLOG("zsw_info") << step_info << " - "; }
    if(isKeepJpts(pt_val, jpt.pt_, bound_tris, all_jpts, jpts_update)) {
      merge_pt_ptr= &jpt; break;
    }
  }
  if(merge_pt_ptr != nullptr) {
    if(step_info %100 == 0) {
      NZSWLOG("zsw_info")  << "bc: collapse edge!!!" << std::endl;
    }
    edgeCollapse(tet_ids, bound_tris, merge_pt_ptr->pt_, vertices_[e.vid_[0]].pt_type_, jpts_update, e, all_jpts,
                 [&eids, &eids_set, this](const size_t e_id){
                   if(vertices_[edges_[e_id].vid_[0]].pt_type_==vertices_[edges_[e_id].vid_[1]].pt_type_
                      && (vertices_[edges_[e_id].vid_[0]].pt_type_==OUTER_POINT || vertices_[edges_[e_id].vid_[0]].pt_type_==INNER_POINT)
                      && eids_set.find(e_id)!=eids_set.end()) {
                     eids.push(e_id); eids_set.insert(e_id); }
                 });
    vertices_[e.vid_[0]].qem_ = cur_qem;
  } else if(step_info %100 == 0) {      NZSWLOG("zsw_info")  << "bc: no poper merge point!"<< std::endl;    }
#endif
}

void zsw::Triangulation::tryCollapseZeroEdge(const size_t e_id,
                                             std::queue<size_t> &eids,
                                             std::set<size_t> eids_set)
{
  Edge &e = edges_[e_id];
  std::vector<size_t> tet_ids;
  for(size_t tid : vertices_[e.vid_[0]].tet_ids_) { tet_ids.push_back(tid); }
  for(size_t tid : vertices_[e.vid_[1]].tet_ids_) { tet_ids.push_back(tid); }
  sort(tet_ids.begin(), tet_ids.end());
  auto unq_end = unique(tet_ids.begin(), tet_ids.end());
  tet_ids.resize(std::distance(tet_ids.begin(), unq_end));

  std::list<JudgePoint> all_jpts;
  std::list<Eigen::Matrix<size_t,3,1>> bound_tris;
  KernelRegionJudger krj;
  for(size_t tid : tet_ids) {
    all_jpts.insert(all_jpts.end(), tets_[tid].jpts_.begin(), tets_[tid].jpts_.end());
    // bound tri
    size_t vcnt=0;
    Eigen::Matrix<size_t,4,1> tmp_bound_tri;
    for(size_t tvid : tets_[tid].vid_) {
      if(tvid == e.vid_[0] || tvid==e.vid_[1]) {          tmp_bound_tri[3]=tvid;        }
      else {          tmp_bound_tri[vcnt++]=tvid;        }
    }
    if(vcnt==3) {
      krj.addConstraint(vertices_[tmp_bound_tri[0]].pt_, vertices_[tmp_bound_tri[1]].pt_,
                        vertices_[tmp_bound_tri[2]].pt_, vertices_[tmp_bound_tri[3]].pt_);
      bound_tris.push_back(tmp_bound_tri.block<3,1>(0,0));
    }
  }
  std::list<Eigen::Matrix<zsw::Scalar,3,1>> candicate_pts;
  // omp_lock_t data_lock;
  // omp_init_lock(&data_lock);
  //#pragma omp parallel for
  for(size_t i=0; i<tet_ids.size(); ++i) {
    const size_t *vid = tets_[tet_ids[i]].vid_;
    std::list<Eigen::Matrix<zsw::Scalar,3,1>> tmp_candicate_pts;
    sampleTet(vertices_[vid[0]].pt_, vertices_[vid[1]].pt_,
              vertices_[vid[2]].pt_, vertices_[vid[3]].pt_, tet_sample_r_, tmp_candicate_pts);
    tmp_candicate_pts.remove_if(std::bind(&KernelRegionJudger::judge, &krj, std::placeholders::_1));

    // omp_set_lock(&data_lock);
    candicate_pts.splice(candicate_pts.end(), tmp_candicate_pts);
    // omp_unset_lock(&data_lock);
  }
  // omp_destroy_lock(&data_lock);
  // check if keepJpts
  std::vector<std::pair<size_t,zsw::Scalar>> jpts_update;
  const Eigen::Matrix<zsw::Scalar,3,1> *merge_pt_ptr = nullptr;
  for(const Eigen::Matrix<zsw::Scalar,3,1> &pt : candicate_pts) {
    if(isKeepJpts(0, pt, bound_tris, all_jpts, jpts_update)) {  merge_pt_ptr = &pt; break;  }
  }
  if(merge_pt_ptr != nullptr) {
    NZSWLOG("zsw_info")  << "zc: collapse edge " << e.vid_[0] << " " << e.vid_[1]
                         << "to " << merge_pt_ptr->transpose() << std::endl;
    edgeCollapse(tet_ids, bound_tris, *merge_pt_ptr, zsw::ZERO_POINT, jpts_update,
                 e, all_jpts, [&eids, &eids_set,this](const size_t e_id) {
                   if(vertices_[edges_[e_id].vid_[0]].pt_type_==vertices_[edges_[e_id].vid_[1]].pt_type_
                     && vertices_[edges_[e_id].vid_[0]].pt_type_==ZERO_POINT
                      && eids_set.find(e_id)!=eids_set.end()){
                     eids.push(e_id); eids_set.insert(e_id);} });
    if(isGood()!=0) { abort(); }
  } else { NZSWLOG("zsw_info")  << "zc: no poper merge point!"<< std::endl; }
}


void zsw::Triangulation::writeBoundTris(const std::string &filepath, const size_t vid0, const size_t vid1)
{
  std::unordered_set<size_t> tet_ids;
  for(size_t tid : vertices_[vid0].tet_ids_) { tet_ids.insert(tid); }
  for(size_t tid : vertices_[vid1].tet_ids_) { tet_ids.insert(tid); }
  std::list<Eigen::Matrix<size_t,3,1>> bound_tris;
  Eigen::Matrix<size_t,4,1> tmp_bound_tri;
  for(size_t tid : tet_ids) {
    size_t vcnt = 0;
    for(size_t tvid : tets_[tid].vid_) {
      if(tvid == vid0 || tvid==vid1) {          tmp_bound_tri[3]=tvid;        }
      else {          tmp_bound_tri[vcnt++]=tvid+1;        }
    }
    if(vcnt==3) {
      bound_tris.push_back(tmp_bound_tri.block<3,1>(0,0));
    }
  }
  std::ofstream ofs;
  OPEN_STREAM(filepath, ofs, std::ofstream::out, return);
  for(const Vertex &vertex : vertices_) {
    ofs << "v " << vertex.pt_.transpose() << std::endl;
  }
  for(auto b_tr : bound_tris) {
    ofs << "f " << b_tr.transpose() << std::endl;
  }
  ofs.close();
}

#define ADD_NORMAL_CONSTRAINT(cur_vid, oth_vid) do{                     \
    for(size_t tid : vertices_[cur_vid].tet_ids_) {                     \
      size_t e_type_cnt=0, opposite_type_cnt=0;                         \
      size_t tmp_vid[3];                                                \
      for(size_t vid : tets_[tid].vid_) {                               \
        if(vid==cur_vid || vid==oth_vid) { continue; }                  \
        if(vertices_[vid].pt_type_ == e_pt_type) { tmp_vid[e_type_cnt]=vid; ++e_type_cnt; } \
        else if(vertices_[vid].pt_type_ == opposite_pt_type) { ++opposite_type_cnt; } \
      }                                                                 \
      if(e_type_cnt==2 && opposite_type_cnt==1) {                       \
        Eigen::Matrix<zsw::Scalar,3,1> tmp_normal =                     \
          (vertices_[tmp_vid[0]].pt_-vertices_[cur_vid].pt_).cross(vertices_[tmp_vid[1]].pt_-vertices_[cur_vid].pt_); \
        tmp_normal.normalize();                                         \
        ncj.addConstraint(vertices_[tmp_vid[0]].pt_, vertices_[tmp_vid[1]].pt_, tmp_normal); \
      }                                                                 \
    }                                                                   \
  }while(0)

void zsw::Triangulation::initNormalCond(NormalConditionJudger &ncj, const Edge &e) const
{
  // bound faces
  zsw::PointType e_pt_type = vertices_[e.vid_[0]].pt_type_;
  zsw::PointType opposite_pt_type = (e_pt_type==zsw::OUTER_POINT) ? zsw::INNER_POINT : zsw::OUTER_POINT;
  ADD_NORMAL_CONSTRAINT(e.vid_[0], e.vid_[1]);
  ADD_NORMAL_CONSTRAINT(e.vid_[1], e.vid_[0]);
}

size_t zsw::Triangulation::isGood() const
{
  // tet is good
  for(const Tet &tet : tets_) {
    if(!tet.valid_) { continue; }
    size_t vo_cnt=0, vi_cnt=0;
    for(size_t vid : tet.vid_) {
      if(!vertices_[vid].valid_) {
        NZSWLOG("zsw_err") << "There is valid tet include invalid vertex!!!" << std::endl;
        return __LINE__;
      }
      if(vertices_[vid].pt_type_ == OUTER_POINT) { vo_cnt++; }
      else if(vertices_[vid].pt_type_ == INNER_POINT){ vi_cnt++;  }
    }
    if(vo_cnt!=0 && vi_cnt!=0 && vo_cnt+vi_cnt!=4) {
      NZSWLOG("zsw_error")  << "There is tet include inner, outer, bbox vertex!!!" << std::endl;
      return __LINE__;
    }
  }

  // edge is good
  std::set<std::pair<size_t,size_t>,zsw::PairCompFunc> e_set_tr(pairComp), e_set_check(pairComp);
  size_t ev_cnt=0;
  for(const zsw::Edge &e : edges_) {
    if(!e.valid_) { continue; }
    ++ev_cnt;
    std::pair<size_t,size_t> ep;
    ep.first=e.vid_[0]; ep.second=e.vid_[1];
    if(ep.first>ep.second) { std::swap(ep.first, ep.second); }
    e_set_tr.insert(ep);
  }
  if(ev_cnt != e_set_tr.size()) {
    NZSWLOG("zsw_err") << "duplicate edge in edges_!!!" << std::endl;
    return __LINE__;
  }
  for(const zsw::Tet &tet : tets_) {
    if(!tet.valid_) { continue; }
    size_t vids[4]; std::copy(tet.vid_, tet.vid_+4, vids); std::sort(vids, vids+4);
    e_set_check.insert(std::pair<size_t,size_t>(vids[0], vids[1]));
    e_set_check.insert(std::pair<size_t,size_t>(vids[0], vids[2]));
    e_set_check.insert(std::pair<size_t,size_t>(vids[0], vids[3]));
    e_set_check.insert(std::pair<size_t,size_t>(vids[1], vids[2]));
    e_set_check.insert(std::pair<size_t,size_t>(vids[1], vids[3]));
    e_set_check.insert(std::pair<size_t,size_t>(vids[2], vids[3]));
  }

  {
    auto it_tr = e_set_tr.begin();
    auto it_check = e_set_check.begin();
    while(it_tr!=e_set_tr.end() && it_check!=e_set_check.end()) {
      if(it_tr->first != it_check->first || it_tr->second!=it_check->second) {
        writeTetMeshAdjVs("/home/wegatron/tmp/debug_adj.vtk", {it_tr->first});
        std::cerr << "it_tr:" << it_tr->first << " " << it_tr->second << std::endl;
        std::cerr << "it_check:" << it_check->first << " " << it_check->second << std::endl;
        std::cerr << "it_tr_type:" << vertices_[it_tr->first].pt_type_ << " " << vertices_[it_tr->second].pt_type_ << std::endl;
        std::cerr << "it_check_type:" << vertices_[it_check->first].pt_type_ << " " << vertices_[it_check->second].pt_type_ << std::endl;
        std::abort();
      }
      ++it_tr; ++it_check;
    }
  }

  if(e_set_check.size() != ev_cnt) {
    NZSWLOG("zsw_err") << "e_set_size error:" << "e_set_tr=" << e_set_tr.size() << ", e_set_check=" << e_set_check.size() << std::endl;
    return __LINE__;
  }

  // vertex(adjcent info) is good
  std::vector<size_t> ve_cnts(vertices_.size(),0);  // vertex's adjcent edge count
  std::vector<size_t> vt_cnts(vertices_.size(),0); // vertex's adjcent tet count

  for(size_t e_id=0; e_id<edges_.size(); ++e_id) {
    const zsw::Edge &edge = edges_[e_id];
    if(!edge.valid_) { continue; }
    ++ve_cnts[edge.vid_[0]];
    ++ve_cnts[edge.vid_[1]];
    for(size_t eei=0; eei<2; ++eei) {
      auto it = std::find(vertices_[edge.vid_[eei]].edge_ids_.begin(), vertices_[edge.vid_[eei]].edge_ids_.end(), e_id);
      if(it==vertices_[edge.vid_[eei]].edge_ids_.end()) {
        NZSWLOG("zsw_err")<<"vertex " << edge.vid_[eei]
                          << "has no adjcent info of edge"<< edge.vid_[0] << ":" << edge.vid_[1] << std::endl;
        return __LINE__;
      }
    }
  }
  for(size_t t_id=0; t_id<tets_.size(); ++t_id) {
    const zsw::Tet &tet = tets_[t_id];
    if(!tet.valid_) { continue; }
    for(size_t vid : tet.vid_) {
      ++vt_cnts[vid];
      auto it = std::find(vertices_[vid].tet_ids_.begin(), vertices_[vid].tet_ids_.end(), t_id);
      if(it==vertices_[vid].tet_ids_.end()) {
        NZSWLOG("zsw_err") << "vertex" << vid << "has no adjcent info of tet" << t_id << std::endl;
        return __LINE__;
      }
    }
  }

  for(size_t i=0; i<vertices_.size(); ++i) {
    if(!vertices_[i].valid_) { continue; }
    if(vertices_[i].edge_ids_.size() > ve_cnts[i]) {
      NZSWLOG("zsw_err") << "vertex" << i << "has extra info of adj edges!!!" << std::endl;
      return __LINE__;
    }
    if(vertices_[i].tet_ids_.size() > vt_cnts[i]) {
      NZSWLOG("zsw_err") << "vertex" << i << "has extra info of adj tets!!!" << std::endl;
      return __LINE__;
    }
  }
  return 0;
}

void zsw::Triangulation::checkTetEdgeExist(const size_t n0, const size_t n1, const size_t n2, const size_t n3)
{
  size_t e_cnt=0;
  for(size_t e_id : vertices_[n0].edge_ids_) {
    if(edges_[e_id].vid_[0]==n1 || edges_[e_id].vid_[0]==n2 ||  edges_[e_id].vid_[0]==n3) { ++e_cnt; }
    if(edges_[e_id].vid_[1]==n1 || edges_[e_id].vid_[1]==n2 ||  edges_[e_id].vid_[1]==n3) { ++e_cnt; }
  }

  for(size_t e_id : vertices_[n1].edge_ids_) {
    if(edges_[e_id].vid_[0]==n2 ||  edges_[e_id].vid_[0]==n3) { ++e_cnt; }
    if(edges_[e_id].vid_[1]==n2 ||  edges_[e_id].vid_[1]==n3) { ++e_cnt; }
  }

  for(size_t e_id : vertices_[n2].edge_ids_) {
    if(edges_[e_id].vid_[0]==n3) { ++e_cnt; }
    if(edges_[e_id].vid_[1]==n3) { ++e_cnt; }
  }

  if(e_cnt!=6) { std::cout << __FILE__ << __LINE__ << std::endl; std::abort(); }
}

void zsw::Triangulation::splitJudgePoints(std::list<zsw::JudgePoint> &jpts, const std::vector<size_t> &tids)
{
  // judge points splice
  const size_t tsize=tids.size();
  for(size_t ti=0; ti<tsize; ++ti) {
    zsw::KernelRegionJudger krj;
    const size_t * tmp_v = tets_[tids[ti]].vid_;
    krj.addConstraint(vertices_[tmp_v[0]].pt_, vertices_[tmp_v[1]].pt_,
                      vertices_[tmp_v[2]].pt_, vertices_[tmp_v[3]].pt_);
    krj.addConstraint(vertices_[tmp_v[0]].pt_, vertices_[tmp_v[1]].pt_,
                      vertices_[tmp_v[3]].pt_, vertices_[tmp_v[2]].pt_);
    krj.addConstraint(vertices_[tmp_v[0]].pt_, vertices_[tmp_v[3]].pt_,
                      vertices_[tmp_v[2]].pt_, vertices_[tmp_v[1]].pt_);
    krj.addConstraint(vertices_[tmp_v[3]].pt_, vertices_[tmp_v[1]].pt_,
                      vertices_[tmp_v[2]].pt_, vertices_[tmp_v[0]].pt_);
    for(auto iter=jpts.begin(); iter!=jpts.end(); ) {
      auto cur_iter=iter++;
      if(krj.judge(cur_iter->pt_)) { tets_[tids[ti]].jpts_.splice(tets_[tids[ti]].jpts_.end(), jpts, cur_iter); }
    }
  }

  for(auto iter=jpts.begin(); iter!=jpts.end();) {
    auto cur_iter=iter++;
    size_t target_ti=0;
    zsw::Scalar min_sq_dis=(cur_iter->pt_ - vertices_[tets_[tids[0]].vid_[0]].pt_).squaredNorm();
    for(size_t ti=0; ti<tsize; ++ti) {
      const size_t * tmp_v = tets_[tids[ti]].vid_;
      zsw::Scalar tmp_s_dis = calcPoint2TriSquaredDis(cur_iter->pt_, vertices_[tmp_v[1]].pt_,
                                                      vertices_[tmp_v[2]].pt_, vertices_[tmp_v[3]].pt_);
      tmp_s_dis = std::min(tmp_s_dis, calcPoint2TriSquaredDis(cur_iter->pt_, vertices_[tmp_v[0]].pt_,
                                                              vertices_[tmp_v[2]].pt_, vertices_[tmp_v[3]].pt_) );
      tmp_s_dis = std::min(tmp_s_dis, calcPoint2TriSquaredDis(cur_iter->pt_, vertices_[tmp_v[0]].pt_,
                                                              vertices_[tmp_v[1]].pt_, vertices_[tmp_v[3]].pt_) );
      tmp_s_dis = std::min(tmp_s_dis, calcPoint2TriSquaredDis(cur_iter->pt_, vertices_[tmp_v[0]].pt_,
                                                              vertices_[tmp_v[1]].pt_, vertices_[tmp_v[2]].pt_) );
      if(min_sq_dis > tmp_s_dis) { target_ti=ti; min_sq_dis = tmp_s_dis; }
    }
    tets_[target_ti].jpts_.splice(tets_[target_ti].jpts_.end(), jpts, cur_iter);
  }
}

#include "triangulation2.h"
#include <unordered_set>

#include <zswlib/mesh/vtk.h>
#include <zswlib/error_ctrl.h>

#include "sampling.h"

void zsw::KernelRegionJudger::addConstraint(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> &v1,
                       const Eigen::Matrix<zsw::Scalar,3,1> &v2, const Eigen::Matrix<zsw::Scalar,3,1> &vr)
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
}

bool zsw::KernelRegionJudger::judge(const Eigen::Matrix<zsw::Scalar,3,1> &pt)
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  return false;
}

zsw::Triangulation::Triangulation(const zsw::Scalar r, std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bo_pts,
                                  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &bi_pts)
{
  assert(bo_pts.size()!=0 && bi_pts.size()!=0);

  Eigen::Matrix<zsw::Scalar,3,2> bbox;
  bbox.block<3,1>(0,0)=bo_pts[0];
  bbox.block<3,1>(0,1)=bo_pts[0];

  size_t pt_id=0;
  std::vector<std::pair<Point, size_t>> tet_points;
  for(const Eigen::Matrix<zsw::Scalar,3,1> &tmp : bo_pts) {
    tet_points.push_back({Point(tmp[0],tmp[1],tmp[2]), pt_id++});
    vertices_.push_back({OUTER_POINT, tmp, {}, {}});
    for(size_t c_i=0; c_i<3; ++c_i) {
      if(tmp[c_i]<bbox(c_i,0)) { bbox(c_i,0)=tmp[c_i];}
      else if(tmp[c_i]>bbox(c_i,1)) { bbox(c_i,1)=tmp[c_i]; }
    }
  }
  for(const Eigen::Matrix<zsw::Scalar,3,1> &tmp : bi_pts) {
    tet_points.push_back({Point(tmp[0],tmp[1],tmp[2]), pt_id++});
    vertices_.push_back({INNER_POINT, tmp, {}, {}});
  }

  // add 8 bbox points
  for(size_t i=0; i<2; ++i) {
    for(size_t j=0; j<2; ++j) {
      for(size_t k=0; k<2; ++k) {
        tet_points.push_back({Point(bbox(0,i), bbox(1,j), bbox(2,k)), pt_id++});
        Eigen::Matrix<zsw::Scalar,3,1> tmp; tmp<<bbox(0,i), bbox(1,j), bbox(2,k);
        vertices_.push_back({BBOX_POINT, tmp, {}, {}});
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
    edges_.push_back({eit->first->vertex(eit->second)->info(), eit->first->vertex(eit->third)->info()});
    Edge &tmp_edge=edges_.back();
    vertices_[tmp_edge.vid_[0]].edge_ids_.push_back(e_id);
    vertices_[tmp_edge.vid_[1]].edge_ids_.push_back(e_id++);
  }
}

void zsw::Triangulation::simpTolerance()
{
  for(Edge e : edges_) {
    // if is bi and bo edge
    if(vertices_[e.vid_[0]].pt_type_!=vertices_[e.vid_[1]].pt_type_ ||
       (vertices_[e.vid_[0]].pt_type_!=OUTER_POINT && vertices_[e.vid_[0]].pt_type_!=INNER_POINT) ) {
      continue; }

    std::unordered_set<size_t> tet_ids;
    for(size_t tid : vertices_[e.vid_[0]].tet_ids_) { tet_ids.insert(tid); }
    for(size_t tid : vertices_[e.vid_[1]].tet_ids_) { tet_ids.insert(tid); }

    KernelRegionJudger krj;
    std::list<JudgePoint> all_jpts;
    for(size_t tid : tet_ids) {
      all_jpts.splice(all_jpts.begin(), tets_[tid].jpts_);

      // add kernel region constraint
      Eigen::Matrix<zsw::Scalar,3,1> tmp_v[4];
      size_t vcnt=0;
      for(size_t tvid : tets_[tid].vid_) {
        if(tvid == e.vid_[0] || tvid==e.vid_[1]) { tmp_v[3]<<vertices_[tvid].pt_[0], vertices_[tvid].pt_[1], vertices_[tvid].pt_[2]; }
        else { tmp_v[vcnt++] << vertices_[tvid].pt_[0], vertices_[tvid].pt_[1], vertices_[tvid].pt_[2]; }
      }
      if(vcnt==3) { krj.addConstraint(tmp_v[0], tmp_v[1], tmp_v[2], tmp_v[3]); }
    }

    // find the candicate points
    std::list<JudgePoint> candicate_pts;
    for(const JudgePoint &jpt : all_jpts) {
      // judge if jpt in kernel region
      if(krj.judge(jpt.pt_)) { candicate_pts.push_back(jpt); }
    }

    // find the best point in the candicate_pts
    zsw::Scalar max_error=-1;
    const JudgePoint *merge_point_ptr=nullptr;
    for(const JudgePoint &jpt : candicate_pts) {
      zsw::Scalar cur_error=fabs(jpt.val_cur_-jpt.val_exp_);
      if(cur_error>max_error && testCollapse(e, jpt.pt_, all_jpts)) {
        max_error=cur_error;
        merge_point_ptr=&jpt;
      }
    }
    edgeCollapse(e, merge_point_ptr->pt_, all_jpts);
  }
}

void zsw::Triangulation::mutualTessellation()
{
  for(Tet &tet : tets_) {
    size_t vo_cnt=0, vi_cnt=0;
    size_t vo[4], vi[4];
    for(size_t vid : tet.vid_) {
      if(vertices_[vid].pt_type_ == OUTER_POINT) { vo[vo_cnt++]=vid; }
      else { vi[vi_cnt++]=vid; }
    }
    // bo : bi = 3 : 1
    if(vo_cnt==3 && vi_cnt==1) {
      tessllelation3v1(vo[0],vo[1],vo[2],vi[0],tet);
    }
    // bo : bi = 2 : 2
    else if(vo_cnt==2 && vo_cnt==vi_cnt) {
      tessllelation2v2(vo[0],vo[1],vi[0],vi[1],tet);
    }
    // bo : bi = 1 : 3
    else if(vo_cnt==1 && vi_cnt==3) {
      tessllelation1v3(vo[0],vi[0],vi[1],vi[2],tet);
    }
  }
}

void zsw::Triangulation::tessllelation3v1(const size_t vo_0, const size_t vo_1,
                                          const size_t vo_2, const size_t vi_0,
                                          Tet &tet)
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  tet.valid_=false;
  size_t nv0=vertices_.size();
  size_t nv1=nv0+1;
  size_t nv2=nv1+1;
  vertices_.push_back({ZERO_POINT,(vertices_[vo_0].pt_+vertices_[vi_0].pt_)/2.0,{},{}}); //nv0
  vertices_.push_back({ZERO_POINT,(vertices_[vo_1].pt_+vertices_[vi_0].pt_)/2.0,{},{}}); //nv1
  vertices_.push_back({ZERO_POINT,(vertices_[vo_2].pt_+vertices_[vi_0].pt_)/2.0,{},{}}); //nv2

  // add tet
  tets_.push_back({true, {vi_0,nv0,nv1,nv2}});
  tets_.push_back({true, {nv0,vo_0,vo_2,vo_1}});
  tets_.push_back({true, {nv1,nv0,vo_2,vo_1}});
  tets_.push_back({true, {nv2,nv1,nv0,vo_2}});
  // add edge
  std::cerr << "Haven't add edge!!!" << std::endl;
  // add additional info of vertex
  std::cerr << "Haven't add additional info of vertex!!!'" << std::endl;
}

void zsw::Triangulation::tessllelation2v2(const size_t vo_0, const size_t vo_1,
                                          const size_t vi_0, const size_t vi_1,
                                          Tet &tet)
{
  tet.valid_=false;
  std::cerr << " invalid the old tet's edge!!!" << std::endl;
  size_t nv0=vertices_.size();
  size_t nv1=nv0+1;
  size_t nv2=nv1+1;
  size_t nv3=nv2+1;
  vertices_.push_back({ZERO_POINT, (vertices_[vo_0].pt_+vertices_[vi_0].pt_)/2.0, {}, {}}); // nv0
  vertices_.push_back({ZERO_POINT, (vertices_[vo_0].pt_+vertices_[vi_1].pt_)/2.0, {}, {}}); // nv1
  vertices_.push_back({ZERO_POINT, (vertices_[vo_1].pt_+vertices_[vi_0].pt_)/2.0, {}, {}}); // nv2
  vertices_.push_back({ZERO_POINT, (vertices_[vo_1].pt_+vertices_[vi_1].pt_)/2.0, {}, {}}); // nv3

  // add tet
  tets_.push_back({true,{nv0,nv1,nv2,vo_0},{}}); //0
  tets_.push_back({true,{nv0,nv1,nv2,vi_0},{}}); //1
  tets_.push_back({true,{nv1,nv2,vi_0,vi_1},{}}); //2
  tets_.push_back({true,{nv1,nv2,nv3,vi_1},{}}); //3
  tets_.push_back({true,{nv2,nv1,vo_0,vo_1},{}}); //4
  tets_.push_back({true,{nv2,nv1,nv3,vo_1},{}}); //5
  // add edge
  edges_.push_back({nv0, vo_0});  edges_.push_back({nv0, nv1});
  edges_.push_back({nv0, nv2});  edges_.push_back({nv0, vi_0});

  edges_.push_back({nv1, vo_0});  edges_.push_back({nv1, vi_0});
  edges_.push_back({nv1, vi_1});  edges_.push_back({nv1, vo_1});
  edges_.push_back({nv1, nv2});  edges_.push_back({nv1, nv3});

  edges_.push_back({nv2, vo_0});  edges_.push_back({nv2, vo_1});
  edges_.push_back({nv2, vi_0});  edges_.push_back({nv2, vi_1});
  edges_.push_back({nv2, nv3});

  edges_.push_back({nv3, vo_1});  edges_.push_back({nv3, vi_1});
  // other vertex info-egde and tet
  std::cerr << "other vertex info-egde and tet missing!!!! " << std::endl;
}

void zsw::Triangulation::tessllelation1v3(const size_t vo_0, const size_t vi_0,
                                          const size_t vi_1, const size_t vi_2,
                                          Tet &tet)
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  tet.valid_=false;
  size_t nv0=vertices_.size();
  size_t nv1=nv0+1;
  size_t nv2=nv1+1;
  vertices_.push_back({ZERO_POINT,(vertices_[vo_0].pt_+vertices_[vi_0].pt_)/2.0,{},{}}); //nv0
  vertices_.push_back({ZERO_POINT,(vertices_[vo_0].pt_+vertices_[vi_1].pt_)/2.0,{},{}}); //nv1
  vertices_.push_back({ZERO_POINT,(vertices_[vo_0].pt_+vertices_[vi_2].pt_)/2.0,{},{}}); //nv2

  // add tets
  tets_.push_back({true, {vo_0,nv0,nv1,nv2},{}});
  tets_.push_back({true, {nv0,vi_0,vi_1,vi_2},{}});
  tets_.push_back({true, {nv1,nv0,vi_1,vi_2},{}});
  tets_.push_back({true, {nv2,nv0,nv1,vi_2},{}});

  // add edge
  std::cerr << "Haven't add edge!!!" << std::endl;
  // add additional vertex info
  std::cerr << "Haven't add  additional vertex info!!!" << std::endl;
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

bool zsw::Triangulation::testCollapse(Edge &e, const Eigen::Matrix<zsw::Scalar,3,1> &pt, std::list<JudgePoint> jpts) const
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
  return false;
}

void zsw::Triangulation::edgeCollapse(Edge &e, const Eigen::Matrix<zsw::Scalar,3,1> &pt, std::list<JudgePoint> jpts)
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
}

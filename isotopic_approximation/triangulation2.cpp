#include "triangulation2.h"
#include <unordered_set>

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
  initTets(delaunay);
}

void zsw::Triangulation::initTets(Delaunay &delaunay)
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
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
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
}

void zsw::Triangulation::writeTetMesh(const std::string &filepath, size_t mask) const
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
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

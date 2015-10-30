#include "triangulation2.h"

zsw::Triangulation::Triangulation(const zsw::Scalar r, std::vector<Point> &bo_points, std::vector<Point> &bi_points)
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
}
void zsw::Triangulation::simpTolerance()
{
  for(Edge e : edges_) {
    if(vertices_[e.vid_[0]].pt_type==OUTER_POINT && vertices_[e.vid_[1]].pt_type==OUTER_POINT) {
      // find the candicate point and do edge collapse
      std::list<JudgePoint> all_jpts;
      unorded_set<size_t> tet_ids;
      for(size_t tid : vertices_[e.vid_[0]].tet_ids_) { tet_ids.insert(tid); }
      for(size_t tid : vertices_[e.vid_[1]].tet_ids_) { tet_ids.insert(tid); }
      KernelRegionJudger krj;
      for(size_t tid : tet_ids) {
        all_jpts.splice(all_jpts.begin(), tets_[tid].jpts_);
        // tet to v0 v1 v2 vr
        Eigen::Matrix<zsw::Scalar,3,1> tmp_v[4];
        size_t vcnt=0;
        for(size_t tvid : tets_[tid].vid_) {
          if(tvid == e.vid_[0] || tvid==e.vid_[1]) { tmp_v[3]<<vertices_[tvid].pt_[0], vertices_[tvid].pt_[1], vertices_[tvid].pt_[2]; }
          else { tmp_v[vcnt++] << vertices_[tvid].pt_[0], vertices_[tvid].pt_[1], vertices_[tvid].pt_[2]; }
        }
        if(vcnt==3) { krj.addConstraint(tmp_v[0], tmp_v[1], tmp_v[2], tmp_v[3]); }
      }
      std::list<JudgePoint> candicate_pts;
      for(const JudgePoint &jpt : all_jpts) {
        // judge if jpt in kernel region
      }
    } else if(vertices_[e.vid_[0]].pt_type==INNER_POINT && vertices_[e.vid_[1]].pt_type==INNER_POINT) {

    }
  }
}

void zsw::Triangulation::mutualTessellation()
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
}

void zsw::Triangulation::writeTetMesh(const string &filepath, size_t mask)
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
}

void zsw::Triangulation::writeSurface(const string &filepath, PointType pt_tyte)
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
}

void zsw::Triangulation::edgeCollapse(size_t eid, const Point &pt)
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
}

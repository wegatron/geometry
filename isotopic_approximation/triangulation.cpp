#include "triangulation.h"

//#include <algorithm>
#include <fstream>
#include <algorithm>
#include <zswlib/mesh/vtk.h>
#include <zswlib/error_ctrl.h>
#include <zswlib/zsw_clock_c11.h>
#include <zswlib/zsw_stl_ext.h>
#include "optimizer.h"
#include "sampling.h"

using namespace std;

namespace zsw
{
  TetMesh::TetMesh(const vector<Point> bz_points, const vector<Point> &bo_points, const vector<Point> &bi_points)
  {
    // add points
    for(const Point &pt : bz_points) { vertices_.push_back({0, pt, -1, {},{}}); }
    for(const Point &pt : bo_points) { vertices_.push_back({1, pt, -1, {},{}}); }
    for(const Point &pt : bi_points) {  vertices_.push_back({-1, pt, -1, {},{}}); }

    // 3d triangulation
    std::vector<pair<Point, size_t>> tmo; // outer tetmesh and inner tetmesh
    size_t cnt = -1;
    for(const Point &pt : bz_points) {  tmo.push_back(pair<Point, size_t>(pt, ++cnt)); }
    vector<pair<Point, size_t>> tmi(tmo);

    for(const Point &pt : bo_points) { tmo.push_back(pair<Point, size_t>(pt, ++cnt)); }
    for(const Point &pt : bi_points) { tmi.push_back(pair<Point, size_t>(pt, ++cnt)); }

    Delaunay to(tmo.begin(), tmo.end());
    Delaunay ti(tmi.begin(), tmi.end());

    ASSURE(to.is_valid(), exit(__LINE__));
    ASSURE(ti.is_valid(), exit(__LINE__));

    // add tets
    size_t tet_id=0;
    addTets(to, tet_id);
    std::cerr << "[INFO] stage 1 finished! tet size:" << tet_id << std::endl;
    addTets(ti, tet_id);
    std::cerr << "[INFO] stage 2 finished! tet size:" << tet_id << std::endl;
    addEdges(ti,to);
    std::cerr << "[INFO] stage 3 finish add edges!" << std::endl;
  }

  void TetMesh::addEdges(const Delaunay &ti, const Delaunay &to)
  {
    // add edges
    size_t edge_id=0;
    for(Delaunay::Finite_edges_iterator eit=to.finite_edges_begin();
        eit!=to.finite_edges_end(); ++eit) {
      const size_t vind0=eit->first->vertex(eit->second)->info();
      const size_t vind1=eit->first->vertex(eit->third)->info();
      // std::cerr << "add edge:" << eit->second << ":" << eit->third << std::endl;
      edges_.push_back({vind0, vind1, 0, true, {}});
      vertices_[vind0].edge_ids_.push_back(edge_id);
      vertices_[vind1].edge_ids_.push_back(edge_id);
      ++edge_id;
    }

    for(Delaunay::Finite_edges_iterator eit=ti.finite_edges_begin();
        eit!=ti.finite_edges_end(); ++eit) {
      if(vertices_[eit->second].pt_type_==0 && vertices_[eit->third].pt_type_==0) { continue; }
      edges_.push_back({eit->second, eit->third, 0, true, {}});
      vertices_[eit->second].edge_ids_.push_back(edge_id);
      vertices_[eit->third].edge_ids_.push_back(edge_id);
      ++edge_id;
    }

    for(Edge & te : edges_) {
      updateFv(te);
    }
  }

  void TetMesh::updateFv(Edge &edge)
  {
    set<size_t> te_fv;
    for(size_t tet_id : vertices_[edge.vind0_].tet_ids_) {
      if(!tets_[tet_id].valid_) { continue; }
      if(tets_[tet_id].vind0_==edge.vind1_ || tets_[tet_id].vind1_==edge.vind1_ || tets_[tet_id].vind2_==edge.vind1_
         || tets_[tet_id].vind3_==edge.vind1_) { // a tet with v0 and v1
        Tet &tt=tets_[tet_id];
        if(tt.vind0_ != edge.vind0_ && tt.vind0_ !=edge.vind1_) { te_fv.insert(tt.vind0_); }
        if(tt.vind1_ != edge.vind0_ && tt.vind1_ !=edge.vind1_) { te_fv.insert(tt.vind1_); }
        if(tt.vind2_ != edge.vind0_ && tt.vind2_ !=edge.vind1_) { te_fv.insert(tt.vind2_); }
        if(tt.vind3_ != edge.vind0_ && tt.vind3_ !=edge.vind1_) { te_fv.insert(tt.vind3_); }
      }
    }
    edge.fv_cnt_=te_fv.size();
    edge.fv_.resize(edge.fv_cnt_);
    copy(te_fv.begin(), te_fv.end(), edge.fv_.begin());
  }

  void TetMesh::addTets(const Delaunay &td, size_t &tet_id)
  {
    for(Delaunay::Finite_cells_iterator cit=td.finite_cells_begin();
        cit!=td.finite_cells_end(); ++cit) {
      if(vertices_[cit->vertex(0)->info()].pt_type_ ==vertices_[cit->vertex(1)->info()].pt_type_ &&
         vertices_[cit->vertex(1)->info()].pt_type_ ==vertices_[cit->vertex(2)->info()].pt_type_ &&
         vertices_[cit->vertex(2)->info()].pt_type_ ==vertices_[cit->vertex(3)->info()].pt_type_ )  {
        continue;
      }
      Tet tmp_tet = {cit->vertex(0)->info(), cit->vertex(1)->info(), cit->vertex(2)->info(), cit->vertex(3)->info(), true};
      tets_.push_back(tmp_tet);
      vertices_[tmp_tet.vind0_].tet_ids_.push_back(tet_id);
      vertices_[tmp_tet.vind1_].tet_ids_.push_back(tet_id);
      vertices_[tmp_tet.vind2_].tet_ids_.push_back(tet_id);
      vertices_[tmp_tet.vind3_].tet_ids_.push_back(tet_id);
      ++tet_id;
    }
  }


  void TetMesh::simplify(const zsw::Scalar tol)
  {
    // collapse ZEdges
    bool collapsable = true;
    size_t collapse_cnt=0;
    while(collapsable) {
      collapsable = false;
      // here the edges_'s vector size should not change
      for(Edge &edge : edges_) {
        if(!edge.valid_) { continue; }
        if(vertices_[edge.vind0_].pt_type_==0 && vertices_[edge.vind1_].pt_type_==0
           && collapseEdge(edge, tol)) { collapsable=true; ++collapse_cnt; }
      }

      std::cout << "Collapsed times:" << collapse_cnt << std::endl;
#ifdef COLLAPSE_ONCE
      collapsable=false;
#endif
    }
  }

  bool TetMesh::collapseEdge(Edge &edge, const zsw::Scalar tol)
  {
    /**
     *check if edge satisfy the link condition: the linked points of v0 and v1 is a valid triangle
                      **/
    // linked vids
    set<size_t> vid_set;
    set<size_t> linked_vids;
    for(size_t tet_id : vertices_[edge.vind0_].tet_ids_) {
      if(!tets_[tet_id].valid_) { continue; }
      vid_set.insert(tets_[tet_id].vind0_);
      vid_set.insert(tets_[tet_id].vind1_);
      vid_set.insert(tets_[tet_id].vind2_);
      vid_set.insert(tets_[tet_id].vind3_);
    }
    for(size_t tet_id : vertices_[edge.vind1_].tet_ids_) {
      if(!tets_[tet_id].valid_) { continue; }
      if(vid_set.find(tets_[tet_id].vind0_)!=vid_set.end()) { linked_vids.insert(tets_[tet_id].vind0_); }
      if(vid_set.find(tets_[tet_id].vind1_)!=vid_set.end()) { linked_vids.insert(tets_[tet_id].vind1_); }
      if(vid_set.find(tets_[tet_id].vind2_)!=vid_set.end()) { linked_vids.insert(tets_[tet_id].vind2_); }
      if(vid_set.find(tets_[tet_id].vind3_)!=vid_set.end()) { linked_vids.insert(tets_[tet_id].vind3_); }
    }

    //linked_vids.erase(edge.vind0_); linked_vids.erase(edge.vind1_);
    if(linked_vids.size() -2!= edge.fv_cnt_) {
      assert(linked_vids.size()-2 > edge.fv_cnt_);
      return false;
    }

    // kernel region sampling point
    Eigen::Matrix<zsw::Scalar,3,1> res_pt;
    if(!findKernelRegionPoint(edge, res_pt, tol)) { return false; }

    // do collapse edge
    Point pt(res_pt[0], res_pt[1], res_pt[2]);
    collapseEdge(edge, pt);
    return true;
  }

  bool TetMesh::findKernelRegionPoint(const Edge &edge, Eigen::Matrix<zsw::Scalar,3,1> &pt, const zsw::Scalar tol) const
  {
#ifdef FAKE_KERNEL_REGION_POINT
    std::cerr << "FAKE_KERNEL_REGION_POINT" << std::endl;
    pt= Point((vertices_[edge.vind0_].pt_[0] + vertices_[edge.vind1_].pt_[0])/2,
              (vertices_[edge.vind0_].pt_[1] + vertices_[edge.vind1_].pt_[1])/2,
              (vertices_[edge.vind0_].pt_[2] + vertices_[edge.vind1_].pt_[2])/2);
    return true;
#else
    // find the one ring faces, and add constraint
    Eigen::Matrix<Ipopt::Number,3,1> cx;
    cx << (vertices_[edge.vind0_].pt_[0]+vertices_[edge.vind1_].pt_[0])/2.0,
      (vertices_[edge.vind0_].pt_[1]+vertices_[edge.vind1_].pt_[1])/2.0,
      (vertices_[edge.vind0_].pt_[2]+vertices_[edge.vind1_].pt_[2])/2.0;
    Ipopt::SmartPtr<zsw::Optimizer> opt= new zsw::Optimizer(cx);
    addOneRingConstraint(edge.vind0_, edge.vind1_, *opt);
    addOneRingConstraint(edge.vind1_, edge.vind0_, *opt);

    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    app->Options()->SetNumericValue("tol", tol);
    app->Options()->SetIntegerValue("max_iter", 100);
    app->Options()->SetIntegerValue("print_level", 1);
    app->Options()->SetIntegerValue("print_frequency_iter", 100);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("output_file", "ipopt.out");

    if(app->Initialize() != Ipopt::Solve_Succeeded) {
      std::cerr << "Error during initialization!" << std::endl;
      return false;
    }

    // Ask Ipopt to solve the problem
    if(app->OptimizeTNLP(opt) == Ipopt::Solve_Succeeded) {
      std::cout << "Solved fine!!!" << std::endl;
      opt->getResult(pt);
      return true;
    }
#endif
    return false;
  }

  void TetMesh::addOneRingConstraint(const size_t vid, const size_t evid, zsw::Optimizer &opt) const
  {
    for(size_t tet_id : vertices_[vid].tet_ids_) {
      if(tets_[tet_id].vind0_ == evid || tets_[tet_id].vind1_ == evid
         || tets_[tet_id].vind2_ == evid || tets_[tet_id].vind3_ == evid) {
        continue;
      }
      size_t face_vids[3];
      size_t *ptr = face_vids;
      if(tets_[tet_id].vind0_ != vid) { *ptr=tets_[tet_id].vind0_; ++ptr; }
      if(tets_[tet_id].vind1_ != vid) { *ptr=tets_[tet_id].vind1_; ++ptr; }
      if(tets_[tet_id].vind2_ != vid) { *ptr=tets_[tet_id].vind2_; ++ptr; }
      if(tets_[tet_id].vind3_ != vid) { *ptr=tets_[tet_id].vind3_;}
      Eigen::Matrix<zsw::Scalar,3,1> v0, v1, v2, vr;
      // same side
      v0 << vertices_[face_vids[0]].pt_[0], vertices_[face_vids[0]].pt_[1], vertices_[face_vids[0]].pt_[2];
      v1 << vertices_[face_vids[1]].pt_[0], vertices_[face_vids[1]].pt_[1], vertices_[face_vids[1]].pt_[2];
      v2 << vertices_[face_vids[2]].pt_[0], vertices_[face_vids[2]].pt_[1], vertices_[face_vids[2]].pt_[2];
      vr << vertices_[vid].pt_[0], vertices_[vid].pt_[1], vertices_[vid].pt_[2];
      Eigen::Matrix<zsw::Scalar,3,1> n = (v1-v0).cross(v2-v0);
      Eigen::Matrix<zsw::Scalar,3,1> m = vr-v0;
      zsw::Scalar k=n.dot(m);
      Eigen::Matrix<zsw::Scalar,3,1> kn = k * n;
      opt.addConstraint(kn[0],kn[1],kn[2], kn.dot(v0)); // bug fix replace k*kn --> kn
    }
  }

  void TetMesh::collapseEdge(Edge &e, const Point &pt)
  {
    e.valid_=false;
    // update vertices_ and tets_
    vertices_[e.vind1_].father_=e.vind0_;
    vertices_[e.vind0_].pt_=pt;
    vector<size_t> ntet_ids;
    for(size_t tet_id : vertices_[e.vind0_].tet_ids_) {
      if(!tets_[tet_id].valid_) { continue; }
      if(tets_[tet_id].vind0_==e.vind1_ || tets_[tet_id].vind1_==e.vind1_ ||
         tets_[tet_id].vind2_==e.vind1_ || tets_[tet_id].vind3_==e.vind1_) { tets_[tet_id].valid_=false; continue; }
      ntet_ids.push_back(tet_id);
    }
    for(size_t tet_id : vertices_[e.vind1_].tet_ids_) {
      if(!tets_[tet_id].valid_) { continue; }
      // if(tets_[tet_id].vind0_==e.vind0_ || tets_[tet_id].vind1_==e.vind0_ ||
      //    tets_[tet_id].vind2_==e.vind0_ || tets_[tet_id].vind3_==e.vind0_) { tets_[tet_id].valid_=false; continue; }
      if(tets_[tet_id].vind0_ == e.vind1_) { tets_[tet_id].vind0_=e.vind0_; }
      else if(tets_[tet_id].vind1_==e.vind1_) { tets_[tet_id].vind1_=e.vind0_; }
      else if(tets_[tet_id].vind2_==e.vind1_) { tets_[tet_id].vind2_=e.vind0_; }
      else if(tets_[tet_id].vind3_==e.vind1_) { tets_[tet_id].vind3_=e.vind0_; }
      ntet_ids.push_back(tet_id);
    }
    vertices_[e.vind0_].tet_ids_=ntet_ids;

    // update edge's vind1_ -> vind0_ ,for te.vind0_ te.vind1_,
    for(size_t e_id : vertices_[e.vind1_].edge_ids_) {
      if(!edges_[e_id].valid_) { continue; }
      size_t tvid = -1;
      if(edges_[e_id].vind0_==e.vind1_)  {
        tvid = edges_[e_id].vind1_;
        edges_[e_id].vind0_=e.vind0_;
      } else {
        tvid = edges_[e_id].vind0_;
        edges_[e_id].vind1_=e.vind0_;
      }
      if(binary_search(e.fv_.begin(), e.fv_.end(), tvid)) {  edges_[e_id].valid_=false; continue;  }
      vertices_[e.vind0_].edge_ids_.push_back(e_id);
    }
    // updatye te.fv_ and update te.fv_cnt_ for e.vind1_ -> e.vind0_
    set<size_t> vid_set;
    for(size_t e_id : vertices_[e.vind1_].edge_ids_) {
      if(!edges_[e_id].valid_) { continue; }
      (edges_[e_id].vind0_!=e.vind1_) ? vid_set.insert(edges_[e_id].vind0_) : vid_set.insert(edges_[e_id].vind1_);
    }
    for(size_t vid : vid_set) {
      for(size_t e_id : vertices_[vid].edge_ids_) {
        if(binary_search(edges_[e_id].fv_.begin(), edges_[e_id].fv_.end(), e.vind1_)) {
          if(!binary_search(edges_[e_id].fv_.begin(), edges_[e_id].fv_.end(), e.vind0_)) {
            // insert e.vind0_ into fv_
            edges_[e_id].fv_.insert(lower_bound(edges_[e_id].fv_.begin(), edges_[e_id].fv_.end(), e.vind0_), e.vind0_);
          } else {  --edges_[e_id].fv_cnt_;  }
        }
      }
    }

    // update vind0' edge which can construct a face with vind1
    vid_set.clear();
    for(size_t e_id : vertices_[e.vind0_].edge_ids_) {
      if(!edges_[e_id].valid_) { continue; }
      (edges_[e_id].vind0_!=e.vind0_) ? vid_set.insert(edges_[e_id].vind0_) : vid_set.insert(edges_[e_id].vind1_);
    }
    for(size_t vid: vid_set) {
      for(size_t e_id : vertices_[vid].edge_ids_) {
        if(binary_search(edges_[e_id].fv_.begin(), edges_[e_id].fv_.end(), e.vind1_)) { updateFv(edges_[e_id]); }
      }
    }
  }

  void TetMesh::writeVtk(const std::string &filepath, const char filter_pt_type) const
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
      // useless tet the same type

      // assert(!(vertices_[tet.vind0_].pt_type_==vertices_[tet.vind1_].pt_type_ &&
      //          vertices_[tet.vind1_].pt_type_==vertices_[tet.vind2_].pt_type_ &&
      //          vertices_[tet.vind2_].pt_type_==vertices_[tet.vind3_].pt_type_));

      if(vertices_[tet.vind0_].pt_type_==filter_pt_type || vertices_[tet.vind1_].pt_type_==filter_pt_type
         || vertices_[tet.vind2_].pt_type_==filter_pt_type || vertices_[tet.vind3_].pt_type_==filter_pt_type) {
        continue;
      }

      tets_data.push_back(tet.vind0_);
      tets_data.push_back(tet.vind1_);
      tets_data.push_back(tet.vind2_);
      tets_data.push_back(tet.vind3_);
    }
    size_t n_pts = pts_data.size()/3;
    size_t n_tets = tets_data.size()/4;
    std::cout << "[INFO] point size:" << n_pts << std::endl;
    std::cout << "[INFO] tet size:" << n_tets << std::endl;
    tet2vtk(ofs, &pts_data[0], n_pts, &tets_data[0], n_tets);
  }

  void TetMesh::cleanVertices()
  {
    // set new index
    const std::vector<Vertex> old_vertices=vertices_;
    vertices_.clear();
    std::vector<size_t> index(old_vertices.size(), -1);
    size_t id = -1;
    for(size_t i=0; i<old_vertices.size(); ++i) {
      if(old_vertices[i].father_ == -1) {
        index[i] = ++id;
        vertices_.push_back(old_vertices[i]);
      }
    }

    std::vector<size_t> id_map(old_vertices.size(), -1);
    for(size_t i=0; i<old_vertices.size(); ++i) {
      // size_t tmp=i;
      // for(; index[tmp]==-1; tmp=old_vertices[tmp].father_);
      // id_map[i]=index[tmp];
      size_t tmp = i;
      while(old_vertices[tmp].father_ != -1) { tmp = old_vertices[tmp].father_ ; }
      id_map[i] = index[tmp];
    }

    // using id_map to update edges and tets
    for(Tet &tet : tets_) {
      if(!tet.valid_) { continue; }
      tet.vind0_ = id_map[tet.vind0_];
      tet.vind1_ = id_map[tet.vind1_];
      tet.vind2_ = id_map[tet.vind2_];
      tet.vind3_ = id_map[tet.vind3_];
    }

    for(Edge &edge : edges_) {
      edge.vind0_ = id_map[edge.vind0_];
      edge.vind1_ = id_map[edge.vind1_];
    }
  }

  void TetMesh::writeSurface(const std::string &filepath, const char pt_type) const
  {
    std::set<Eigen::Matrix<size_t,3,1>, bool (*)(Eigen::Matrix<size_t,3,1>, Eigen::Matrix<size_t,3,1>)>
      faces([](Eigen::Matrix<size_t,3,1> a, Eigen::Matrix<size_t,3,1> b){ return a[0]<b[0] || a[1]<b[1] || a[2]<b[2]; });
    vector<size_t> used_vind_map(vertices_.size(),0);
    size_t used_vid=0;
    for(const Tet &tet : tets_) {
      if(!tet.valid_) { continue; }
      Eigen::Matrix<size_t,4,1> face;
      size_t i=0;
      if(vertices_[tet.vind0_].pt_type_==pt_type) { face[i++]=tet.vind0_; }
      if(vertices_[tet.vind1_].pt_type_==pt_type) { face[i++]=tet.vind1_; }
      if(vertices_[tet.vind2_].pt_type_==pt_type) { face[i++]=tet.vind2_; }
      if(vertices_[tet.vind3_].pt_type_==pt_type) { face[i++]=tet.vind3_; }
      if(i==3) {
        sort(face.data(), face.data()+3);
        faces.insert(face.block<3,1>(0,0));
        used_vind_map[face[0]]=-1;
        used_vind_map[face[1]]=-1;
        used_vind_map[face[2]]=-1;
      }
    }
    ofstream ofs;
    OPEN_STREAM(filepath, ofs, std::ofstream::out, return);

    size_t index=0;
    for(size_t i=0; i<used_vind_map.size(); ++i) {
      if(used_vind_map[i] == -1) {
        ofs << "v " << vertices_[i].pt_[0] << " " << vertices_[i].pt_[1]
            << " " << vertices_[i].pt_[2] << std::endl;
        used_vind_map[i] = ++index;
      }
    }
    for(const Eigen::Matrix<size_t,3,1> &face : faces) {
      ofs << "f " << used_vind_map[face[0]] << " " << used_vind_map[face[1]] << " " << used_vind_map[face[2]] << std::endl;
    }
    ofs.close();
  }

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

  bool zsw::KernelRegionJudger::judge(const Point &pt)
  {
    Eigen::Matrix<zsw::Scalar,3,1> ept; ept << pt[0], pt[1], pt[2];
    for(size_t i=0; i<vec_v0.size(); ++i) {
      if(vec_vn[i].dot(ept-vec_v0[i]) < 0) {
        return false;
      }
    }
    return true;
  }

#ifdef ZSW_DEBUG
  bool TetMesh::testCollapseEdge(size_t vind0, size_t vind1)
  {
    for(Edge &edge : edges_) {
      // std::cerr << "edge:" << edge.vind0_ << " : " << edge.vind1_ << std::endl;
      if((edge.vind0_==vind0 && edge.vind1_==vind1) || (edge.vind0_==vind1 && edge.vind1_==vind0)) {
        collapseEdge(edge, 1e-3);
        std::cout << "edge: " << edge.vind0_ << ":" << edge.vind1_ << "collapsed fine!!!" << std::endl;
        return true;
      }
    }
    return false;
  }
#endif
}

#define _EXPORTING

#include "bilateral_normal_filter.h"

#include <iostream>
#include <fstream>
#include <list>
#include <queue>
#include <unordered_set>
#include <boost/foreach.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include "EdgeFlipMean.h"
#include "zsw_clock.h"
#include "util.h"
#include "const_val.h"

using namespace std;

ZSW_API zsw::BilateralNormalFilter::BilateralNormalFilter(std::shared_ptr<zsw::BNFParam> param)
{
  param_=param;
}

ZSW_API zsw::Scalar zsw::BilateralNormalFilter::calBc(const size_t fid)
{
  Eigen::Vector3d ci = fc_.block<3,1>(0,fid);
  zsw::Scalar dis = 0.0;
  assert(ref_face_set_[fid].size()!=0);
  BOOST_FOREACH(size_t t_fid, ref_face_set_[fid]) {
    Eigen::Vector3d cj=fc_.block<3,1>(0,t_fid);
    dis += (ci-cj).norm();
  }
  dis /= ref_face_set_[fid].size();
  return 2*dis*dis;
}

ZSW_API void zsw::BilateralNormalFilter::smooth()
{
  std::vector<zsw::FakeSet<size_t>> tmp_ring;
  rRingFacesByVertex(*tm_,1,tmp_ring);

  Eigen::Matrix<zsw::Scalar, 3, Eigen::Dynamic> vertices(3, tm_->n_vertices());
  for(size_t j=0; j<tm_->n_vertices(); ++j) {
    vertices.block<3,1>(0,j)=tm_->point(zsw::mesh::TriMesh::VertexHandle(j));
  }
#pragma omp parallel for
  for(size_t i=0; i<tm_->n_faces(); ++i) {
    const Eigen::Vector3d &ni=tm_->normal(zsw::mesh::TriMesh::FaceHandle(i));
    bool need_smooth=false;
    BOOST_FOREACH(size_t fid, tmp_ring[i]) {
      const Eigen::Vector3d &nj=tm_->normal(zsw::mesh::TriMesh::FaceHandle(fid));
      if(ni.dot(nj)<param_->smooth_threshold_) {
        need_smooth=true;
        break;
      }
    }
    if(need_smooth) {
      for(zsw::mesh::TriMesh::ConstFaceVertexIter fv_it = tm_->cfv_begin(zsw::mesh::TriMesh::FaceHandle(i)); fv_it.is_valid(); ++fv_it) {
        Eigen::Vector3d tmp_v = Eigen::Vector3d::Zero();
        size_t vn = 0;
        for(zsw::mesh::TriMesh::ConstVertexVertexIter vv_it = tm_->cvv_begin(fv_it); vv_it.is_valid(); ++vv_it) {
          tmp_v += vertices.block<3,1>(0,vv_it->idx());
          ++vn;
        }
        tm_->set_point(fv_it,0.45/vn*tmp_v+0.55*tm_->point(fv_it));
      }
    }
  }
//  if(param_->need_flipedge_) {
//    Sn3DGraphics::EdgeFlipMean ef; ef.run(*tm_);
//  }
  param_->need_smooth_=false;
  param_->bnf_type_=zsw::BNFParam::BASIC_FE_RING;
  param_->rn_=1;
  param_->st_=2;
  BilateralNormalFilter bnf(param_); bnf.run(tm_);
}

ZSW_API zsw::BilateralNormalFilter::~BilateralNormalFilter() {
}

ZSW_API void zsw::BilateralNormalFilter::preCompute()
{
  // flip edge if needed
  if(param_->need_flipedge_) {
    Sn3DGraphics::EdgeFlipMean ef; ef.run(*tm_);
    //param_->need_flip_edge_ = false;
  }
  // calc face center and face area adn face normal
  zsw::mesh::calcFaceArea(*tm_, fa_);
  zsw::mesh::calcFaceCenter(*tm_, fc_);
  tm_->request_face_normals(); tm_->update_face_normals();

  std::cerr << "compute ring:" << std::endl;
  zsw::common::Clock clock;
  // compute one_ring
  switch(param_->bnf_type_) {
  case zsw::BNFParam::BASIC_FE_RING:
    rRingFacesByEdge(*tm_,param_->rn_,ref_face_set_);
    break;
  case zsw::BNFParam::BASIC_FV_RING:
    rRingFacesByVertex(*tm_,param_->rn_,ref_face_set_);
    break;
  case zsw::BNFParam::EXTENDED_RADIUS_ERING:
    dynamicRringByEdge(*tm_,param_->normal_threshold_,param_->r_min_,param_->r_max_,fc_,ref_face_set_);
    break;
  case zsw::BNFParam::EXTENDED_RADIUS_VRING:
    dynamicRringByVertex(*tm_,param_->normal_threshold_,param_->r_min_,param_->r_max_,fc_,ref_face_set_);
    break;
  default:
    std::cerr << "Un regonize RING!" << std::endl;
    break;
  }
  std::cout << "compute ring using:" << clock.time() << " milliseconds\n" << std::endl;
  std::list<Eigen::Triplet<zsw::Scalar> > triplet_list;
#pragma omp parallel for
  for(size_t i=0; i<tm_->n_faces(); ++i) {
    const zsw::Scalar b_c=calBc(i);
    const Eigen::Vector3d ci=fc_.block<3,1>(0,i);
    const Eigen::Vector3d ni=tm_->normal(zsw::mesh::TriMesh::FaceHandle(i));
    BOOST_FOREACH(size_t fid, ref_face_set_[i]) {
      const Eigen::Vector3d cj=fc_.block<3,1>(0,fid);
      const Eigen::Vector3d nj=tm_->normal(zsw::mesh::TriMesh::FaceHandle(fid));

      const zsw::Scalar w_tmp = fa_(fid)*pow(zsw::const_val::e, -(cj-ci).dot(cj-ci)/b_c)*pow(zsw::const_val::e, -(nj-ni).dot(nj-ni)/param_->bs_);
#pragma omp critical
      triplet_list.push_back(Eigen::Triplet<zsw::Scalar>(i, fid, w_tmp));
    }
  }
  weight_.resize(tm_->n_faces(),tm_->n_faces());
  weight_.setFromTriplets(triplet_list.begin(), triplet_list.end());
#if 0
  std::cout << weight_ << std::endl;
  exit(0);
#endif
}

ZSW_API void zsw::BilateralNormalFilter::filterNormal() {
  Eigen::Matrix<zsw::Scalar, 3, Eigen::Dynamic> tmp_normal(3,tm_->n_faces());
#pragma omp parallel for
  for(size_t i=0; i<tm_->n_faces(); ++i) {
      tmp_normal.block<3,1>(0, i) = tm_->normal(zsw::mesh::TriMesh::FaceHandle(i));
    }
#pragma omp parallel for
  for(size_t i=0; i<tm_->n_faces(); ++i) {
      Eigen::Vector3d new_ni = Eigen::MatrixXd::Zero(3,1);
      zsw::Scalar wa = 0.0;
      BOOST_FOREACH(size_t fid, ref_face_set_[i]) {
        wa += weight_.coeff(i, fid);
        new_ni += tmp_normal.block<3,1>(0,fid)*weight_.coeff(i, fid);
      }
      new_ni/=wa;
      if(new_ni.norm()>1e-8) {
          tm_->set_normal(zsw::mesh::TriMesh::FaceHandle(i), new_ni/new_ni.norm());
        }
    }
}
ZSW_API void zsw::BilateralNormalFilter::updateVertex() {
#pragma omp parallel for
  for(zsw::mesh::TriMesh::VertexIter v_it = tm_->vertices_begin();
      v_it!=tm_->vertices_end(); ++v_it) {
      for(zsw::mesh::TriMesh::ConstVertexFaceIter vf_it = tm_->cvf_begin(v_it);
          vf_it.is_valid(); ++vf_it) {
          size_t vid[3] = { v_it->idx(), -1, -1};
          size_t ind = 0;
          for(zsw::mesh::TriMesh::ConstFaceVertexIter fv_it = tm_->cfv_begin(vf_it); fv_it.is_valid(); ++fv_it) {
              if(fv_it->idx() != vid[0]) {
                  vid[++ind] = fv_it->idx();
                }
            }
          Eigen::Vector3d edge[2] = { tm_->point(zsw::mesh::TriMesh::VertexHandle(vid[1]))-tm_->point(zsw::mesh::TriMesh::VertexHandle(vid[0])),
                                    tm_->point(zsw::mesh::TriMesh::VertexHandle(vid[2]))-tm_->point(zsw::mesh::TriMesh::VertexHandle(vid[0]))};
          Eigen::Vector3d fnormal = tm_->normal(vf_it);
          Eigen::Vector3d pt = tm_->point(v_it);
          Eigen::Vector3d npt = pt+(fnormal.dot(edge[0])+fnormal.dot(edge[1]))/18.0*fnormal;
          tm_->set_point(v_it, npt);
        }
    }
}

void zsw::dynamicRringByEdge(const zsw::mesh::TriMesh &trimesh,
                               const zsw::Scalar normal_threshold,
                               const size_t r_min, const size_t r_max,
                               const Eigen::Matrix<zsw::Scalar,3,Eigen::Dynamic> &fc,
                               std::vector<zsw::FakeSet<size_t>> &ring)
{
  ring.clear(); ring.resize(trimesh.n_faces());
#pragma omp parallel for
  for(size_t i=0; i<trimesh.n_faces(); ++i) {
    size_t piquet=i;
    size_t cnt=0;
    std::unordered_set<size_t> tmp_ring;
    std::queue<size_t> qu; qu.push(piquet);
    Eigen::Matrix<zsw::Scalar,3,1> cur_normal=trimesh.normal(zsw::mesh::TriMesh::FaceHandle(i));
    while(!qu.empty() && cnt<r_max) {
      size_t cur_fid=qu.front(); qu.pop();
      for(zsw::mesh::TriMesh::ConstFaceFaceIter ff_it=trimesh.cff_begin(zsw::mesh::TriMesh::FaceHandle(cur_fid));
          ff_it.is_valid(); ++ff_it) {
        if(tmp_ring.find(ff_it->idx())==tmp_ring.end()) {
          qu.push(ff_it->idx());
          if(cnt<r_min || cur_normal.dot(trimesh.normal(ff_it))>normal_threshold) {
            tmp_ring.insert(ff_it->idx());
          }
        }
      } // end for
      if(cur_fid==piquet) {
        if(!qu.empty()) {
          piquet=qu.back();
        }
        ++cnt;
      }
    } // end while
    tmp_ring.erase(i);
    ring[i].initFromSet(tmp_ring);
  } // end for(size_t i=0; i<trimesh_n_faces(); ++i) {
}

void zsw::dynamicRringByVertex(const zsw::mesh::TriMesh &trimesh,
                               const zsw::Scalar normal_threshold,
                               const size_t r_min, const size_t r_max,
                               const Eigen::Matrix<zsw::Scalar,3,Eigen::Dynamic> &fc,
                               std::vector<zsw::FakeSet<size_t>> &ring)
{
  ring.clear(); ring.resize(trimesh.n_faces());
#pragma omp parallel for
  for(size_t i=0; i<trimesh.n_faces(); ++i) {
    size_t piquet=i;
    size_t cnt=0;
    std::unordered_set<size_t> tmp_ring;
    std::queue<size_t> qu; qu.push(piquet);
    Eigen::Matrix<zsw::Scalar,3,1> cur_normal=trimesh.normal(zsw::mesh::TriMesh::FaceHandle(i));
    while(!qu.empty() && cnt<r_max) {
      size_t cur_fid=qu.front(); qu.pop();
      for(zsw::mesh::TriMesh::ConstFaceVertexIter fv_it=trimesh.cfv_begin(zsw::mesh::TriMesh::FaceHandle(cur_fid));
          fv_it.is_valid(); ++fv_it) {
        for(zsw::mesh::TriMesh::ConstVertexFaceIter vf_it=trimesh.cvf_begin(fv_it); vf_it.is_valid(); ++vf_it) {
          if(tmp_ring.find(vf_it->idx())==tmp_ring.end()) {
            qu.push(vf_it->idx());
            if(cnt<r_min || cur_normal.dot(trimesh.normal(vf_it))>normal_threshold) {
              tmp_ring.insert(vf_it->idx());
            }
          }
        }
      }
      if(cur_fid==piquet) {
        if(!qu.empty()) {
          piquet=qu.back();
        }
        ++cnt;
      }
    } // end while
    tmp_ring.erase(i);
    ring[i].initFromSet(tmp_ring);
  } // end for(size_t i=0; i<trimesh_n_faces(); ++i) {
}

ZSW_API void zsw::BilateralNormalFilter::run(std::shared_ptr<zsw::mesh::TriMesh> tm)
{
  tm_=tm;
  zsw::common::Clock clock;
  preCompute();
  std::cerr << "precompute(including flip edge)  using:" << clock.time() << " milliseconds\n" << std::endl;

  for(size_t i=0; i<param_->st_; ++i) {
    std::cout << "smooth normal step " << i << std::endl;
    filterNormal();
  }

  std::cerr << "filter normal using:" << clock.time() << "milliseconds\n" << std::endl;
  for(size_t i=0; i<param_->ut_; ++i) {
    std::cout << "update vertex step " << i << std::endl;
    updateVertex();
  }
  std::cerr << "update vertex using:" << clock.time() << "milliseconds\n" << std::endl;
  if(param_->need_smooth_)
  {
    smooth();
  }
  std::cerr << "smooth using:" << clock.time() << "milliseconds\n" << std::endl;
}

void zsw::rRingFacesByEdge(const zsw::mesh::TriMesh &trimesh, const size_t rn, std::vector<zsw::FakeSet<size_t>> &ring)
{
  ring.clear(); ring.resize(trimesh.n_faces());
#pragma omp parallel for
  for(size_t i=0; i<trimesh.n_faces(); ++i) {
    size_t piquet=i;
    size_t cnt=0;
    std::unordered_set<size_t> tmp_ring;
    std::queue<size_t> qu; qu.push(piquet);
    while(!qu.empty() && cnt<rn) {
      size_t cur_fid=qu.front(); qu.pop();
      for(zsw::mesh::TriMesh::ConstFaceFaceIter ff_it=trimesh.cff_begin(zsw::mesh::TriMesh::FaceHandle(cur_fid));
          ff_it.is_valid(); ++ff_it) {
        if(tmp_ring.find(ff_it->idx())==tmp_ring.end()) {
          qu.push(ff_it->idx());
          tmp_ring.insert(ff_it->idx());
        }
      } // end for
      if(cur_fid==piquet) {
        if(!qu.empty()) {
          piquet=qu.back();
        }
        ++cnt;
      }
    } // end while
    tmp_ring.erase(i);
    ring[i].initFromSet(tmp_ring);
  } // end for(size_t i=0; i<trimesh_n_faces(); ++i) {
}

void zsw::rRingFacesByVertex(const zsw::mesh::TriMesh &trimesh, const size_t rn, std::vector<zsw::FakeSet<size_t>> &ring)
{
  ring.clear(); ring.resize(trimesh.n_faces());
#pragma omp parallel for
  for(size_t i=0; i<trimesh.n_faces(); ++i) {
    size_t piquet=i;
    size_t cnt=0;
    std::unordered_set<size_t> tmp_ring;
    std::queue<size_t> qu; qu.push(piquet);
    while(!qu.empty() && cnt<rn) {
      size_t cur_fid=qu.front(); qu.pop();
      for(zsw::mesh::TriMesh::ConstFaceVertexIter fv_it=trimesh.cfv_begin(zsw::mesh::TriMesh::FaceHandle(cur_fid));
          fv_it.is_valid(); ++fv_it) {
        for(zsw::mesh::TriMesh::ConstVertexFaceIter vf_it=trimesh.cvf_begin(fv_it); vf_it.is_valid(); ++vf_it) {
          if(tmp_ring.find(vf_it->idx())==tmp_ring.end()) {
            qu.push(vf_it->idx());
            tmp_ring.insert(vf_it->idx());
          }
        }
      }
      if(cur_fid==piquet) {
        if(!qu.empty()) {
          piquet=qu.back();
        }
        ++cnt;
      }
    } // end while
    tmp_ring.erase(i);
    ring[i].initFromSet(tmp_ring);
  } // end for(size_t i=0; i<trimesh_n_faces(); ++i) {
}

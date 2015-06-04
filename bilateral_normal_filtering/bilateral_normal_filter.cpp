#include "bilateral_normal_filter.h"

#include <fstream>
#include <zjucad/matrix/io.h>
#include "vtk.h"
using namespace std;

#define E 2.7182818284590

typedef   zjucad::matrix::matrix<double> matrixd;
typedef zjucad::matrix::matrix<size_t> matrixst;

zsw::BilateralNormalFilter::BilateralNormalFilter()
{
  ut_ = 15; // suitable for most cases
  b_s_ = 0.1; // the most important parameter, 0.1 is suitable for preserving features
  st_ = 3; // usually 3, can be larger if noise is in high level, usually around 5
  ring_type_ = ONE_EDGE_RING;
  // b_c_ = 0.3; // should dynamiclly caculated as the average distance of 1-ring face center distance
}

void zsw::BilateralNormalFilter::filter(jtf::mesh::tri_mesh &trimesh)
{
  preCompute(trimesh);
  std::cerr << "precompute done!" << std::endl;
  for(size_t i=0; i<st_; ++i) {
    std::cout << "smooth normal step " << i << std::endl;
    filterNormal(trimesh);
  }

  for(size_t i=0; i<ut_; ++i) {
    std::cout << "update vertex step " << i << std::endl;
    updateVertex(trimesh);
  }
}

void zsw::BilateralNormalFilter::filterNormal(jtf::mesh::tri_mesh &trimesh)
{
  using namespace zjucad::matrix;
  const matrixd &node = trimesh.trimesh_.node_;
  const matrixst &mesh = trimesh.trimesh_.mesh_;
  const matrixd &face_area = trimesh.face_area_;
  matrixd &normal = trimesh.face_normal_;
  matrixd tmp_normal = normal;
  assert(node.size(1) == 3 && mesh.size(1) == 3 && normal.size(2)==mesh.size(2));
  #pragma omp parallel for
  for(size_t i=0; i<mesh.size(2); ++i) {
    double wa = 0.0;
    matrixd new_ni = zjucad::matrix::zeros(3,1);
    for(const size_t fid : one_ring_[i]) {
      const matrixd nj = normal(colon(), fid);
      wa += weight_.coeff(i, fid);
      new_ni += weight_.coeff(i, fid) * nj;
    }
    new_ni /= wa;
    if(zjucad::matrix::norm(new_ni) > 1e-5) {
      tmp_normal(colon(), i) = new_ni/norm(new_ni); // normalize
    }
#if DEBUG
    else {
      std::cerr << "normal so min" << std::endl;
    }
#endif
  }
  normal = tmp_normal;
}

void zsw::BilateralNormalFilter::preCompute(const jtf::mesh::tri_mesh &trimesh)
{
  // compute one_ring
  using namespace zjucad::matrix;
  if(ring_type_ == ONE_EDGE_RING) {
    processEdgeOneRing(trimesh, one_ring_);
  } else if(ring_type_ == ONE_VERTEX_RING) {
    processVertexOneRing(trimesh, one_ring_);
  }

  // compute face center
  const matrixst &mesh = trimesh.trimesh_.mesh_;
  fc_.resize(3, mesh.size(2));
  const matrixd &node = trimesh.trimesh_.node_;
#pragma omp parallel for
  for(size_t i=0;i<mesh.size(2); ++i) {
    fc_(colon(), i) = (node(colon(), mesh(0,i)) + node(colon(), mesh(1,i)) + node(colon(), mesh(2,i)))/3.0;
  }

  const zjucad::matrix::matrix<double> &face_area = trimesh.face_area_;
  const matrixd &normal = trimesh.face_normal_;
  std::list<Eigen::Triplet<double>> triplet_list;
  for(size_t i = 0; i<mesh.size(2); ++i) {
    double b_c = calBc(i, mesh, node);
    matrixd ci = fc_(colon(), i);
    matrixd ni = normal(colon(), i);
    for( size_t fid : one_ring_[i]) {
      const matrixd cj = fc_(colon(), fid);
      const matrixd nj = normal(colon(), fid);
      const double w_tmp = face_area[fid]*pow(E, -dot(cj-ci, cj-ci)/b_c)*pow(E, -dot(nj-ni, nj-ni)/b_s_);
      triplet_list.push_back(Eigen::Triplet<double>(i, fid, w_tmp));
    }
  }
  weight_.resize(mesh.size(2), mesh.size(2));
  weight_.setFromTriplets(triplet_list.begin(), triplet_list.end());
}

void zsw::BilateralNormalFilter::updateVertex(jtf::mesh::tri_mesh &trimesh)
{
  using namespace zjucad::matrix;

  matrixd &node = trimesh.trimesh_.node_;
  matrixst &mesh = trimesh.trimesh_.mesh_;
  matrixd &normal = trimesh.face_normal_;
  for(size_t i=0; i<mesh.size(2); ++i) {
    matrixd v0 = node(colon(), mesh(0,i));
    matrixd v1 = node(colon(), mesh(1,i));
    matrixd v2 = node(colon(), mesh(2,i));
    matrixd ni = normal(colon(), i);
    node(colon(), mesh(0,i)) = v0 + ni*dot(ni, v2+v1-2*v0)/18.0;
    node(colon(), mesh(1,i)) = v1 + ni*dot(ni, v2+v0-2*v1)/18.0;
    node(colon(), mesh(2,i)) = v2 + ni*dot(ni, v0+v1-2*v2)/18.0;
  }
}

double zsw::BilateralNormalFilter::calBc(const size_t fid,
                                         const zjucad::matrix::matrix<size_t> &mesh,
                                         const zjucad::matrix::matrix<double> &node)
{
  using namespace zjucad::matrix;
  matrixd ci = fc_(colon(), fid);
  double dis = 0.0;
  assert(one_ring_[fid].size()!=0);
  for(size_t t_fid : one_ring_[fid]) {
    matrixd cj = fc_(colon(), t_fid);
    dis += norm(ci-cj);
  }
  dis /= one_ring_[fid].size();
  return 2*dis*dis;
}

void zsw::writeTriMesh(const std::string &filename, const zjucad::matrix::matrix<size_t> &mesh,
                       const zjucad::matrix::matrix<double> &node,
                       const zjucad::matrix::matrix<double> &normal)
{
  assert(mesh.size(1) == 3 && mesh.size(2) !=0);
  assert(node.size(1) == 3 && node.size(2) !=0);
  assert(normal.size(1) == 3 && normal.size(2) !=0);
  assert(normal.size(2) == mesh.size(2));

  std::ofstream ofs(filename);
  if(!ofs) {
    cerr << "can open file: " << filename << " for write!" << endl;
    exit(1);
  }
  for(size_t i=0; i<node.size(2); ++i) {
    ofs << "v " << node(0,i) << " " << node(1,i) << " " << node(2,i) << std::endl;
  }

  for(size_t i=0; i<normal.size(2);++i) {
    ofs << "vn " << normal(0,i) << " " << normal(1,i) << " " << normal(2,i) << std::endl;
  }

  for(size_t i=0; i<mesh.size(2); ++i) {
    ofs << "f " << mesh(0,i)+1 << "//" << i  << " " << mesh(1,i)+1 << "//" << i  << " " << mesh(2,i)+1 << "//" << i << std::endl;
  }

  // for(size_t i=0; i<mesh.size(2); ++i) {
  //   ofs << "f " << mesh(0,i)+1 <<  " " << mesh(1,i)+1 << " " << mesh(2,i)+1 << std::endl;
  // }

  ofs.close();
}

void zsw::writeVtk(const std::string &filename, const zjucad::matrix::matrix<size_t> &mesh,
                   const zjucad::matrix::matrix<double> &node)
{
  size_t vn = node.size(2);

  std::ofstream outf(filename.c_str());

  outf.precision(15);

  if (outf.fail()) {
    std::cerr << "# [ ERROR: draw singularity to vtk ] cannot open file: "
              << filename << std::endl;
  }

  tri2vtk(outf, &node[0], vn, &mesh[0], mesh.size(2));
  outf.close();
}

void zsw::processEdgeOneRing(const jtf::mesh::tri_mesh &trimesh, std::vector<set<size_t>> &one_ring)
{
  using namespace zjucad::matrix;
  const matrixst &mesh = trimesh.trimesh_.mesh_;
  one_ring.resize(mesh.size(2));
  #pragma omp parallel for
  for(size_t fid=0; fid<mesh.size(2); ++fid) {
    zjucad::matrix::matrix<size_t> vid = mesh(colon(), fid);
    {
      std::pair<size_t, size_t> res = trimesh.ea_->query(vid[0], vid[1]);
      if(res.first != fid && res.first!=-1)   { one_ring[fid].insert(res.first); }
      else if(res.second!=fid && res.second!=-1) { one_ring[fid].insert(res.second); }
    }

    {
      std::pair<size_t, size_t> res = trimesh.ea_->query(vid[0], vid[2]);
      if(res.first != fid && res.first!=-1)   { one_ring[fid].insert(res.first); }
      else if(res.second!=fid && res.second!=-1) { one_ring[fid].insert(res.second); }
    }

    {
      std::pair<size_t, size_t> res = trimesh.ea_->query(vid[1], vid[2]);
      if(res.first != fid && res.first!=-1)   { one_ring[fid].insert(res.first); }
      else if(res.second!=fid && res.second!=-1) { one_ring[fid].insert(res.second); }
    }
  }
}

void zsw::processVertexOneRing(const jtf::mesh::tri_mesh &trimesh, std::vector<set<size_t>> &one_ring)
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
}

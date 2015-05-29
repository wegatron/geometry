#include "bilateral_normal_filter.h"

#include <fstream>

using namespace std;

#define E 2.7182818284590

typedef   zjucad::matrix::matrix<double> matrixd;
typedef zjucad::matrix::matrix<size_t> matrixst;

zsw::BilateralNormalFilter::BilateralNormalFilter()
{
  st_ = 3;
  ut_ = 16;
  b_c_ = 0.1;
  b_s_ = 0.3;
}
void zsw::BilateralNormalFilter::filter(jtf::mesh::tri_mesh &trimesh)
{
  for(size_t i=0; i<st_; ++i) {
    std::cout << "smooth normal step " << i << std::endl;
    filterNormal(trimesh);
  }
  writeTriMesh("/home/wegatron/tmp/tooth_debug0.obj", trimesh.trimesh_.mesh_, trimesh.trimesh_.node_, trimesh.face_normal_);
  for(size_t i=0; i<ut_; ++i) {
    std::cout << "update vertex step " << i << std::endl;
    updateVertex(trimesh);
  }
  jtf::mesh::save_obj("/home/wegatron/tmp/tooth_res.obj", trimesh.trimesh_.mesh_, trimesh.trimesh_.node_);
}

void zsw::BilateralNormalFilter::filterNormal(jtf::mesh::tri_mesh &trimesh)
{
  using namespace zjucad::matrix;
  matrixd &node = trimesh.trimesh_.node_;
  matrixst &mesh = trimesh.trimesh_.mesh_;
  matrixd &normal = trimesh.face_normal_;
  matrixd &face_area = trimesh.face_area_;
  assert(node.size(1) == 3 && mesh.size(1) == 3 && normal.size(2)==mesh.size(2));
  for(size_t i=0; i<mesh.size(2); ++i) {
    // caculate c_i, n_i
    matrixd ci = ( node(colon(), mesh(0,i)) + node(colon(), mesh(1,i)) + node(colon(), mesh(2,i)) )/3;
    matrixd ni = normal(colon(), i);
    vector<size_t> fid_one_ring;
    if(!queryFidOneRing(i, trimesh, fid_one_ring)) { continue; } // boundary cell
    double wa = 0.0;
    matrixd new_ni = zjucad::matrix::zeros(3,1);
    for(const size_t fid : fid_one_ring) {
      const matrixd cj = ( node(colon(), mesh(0,fid)) + node(colon(), mesh(1,fid)) + node(colon(), mesh(2,fid)) )/3;
      const matrixd nj = normal(colon(), fid);
      const double w_tmp = face_area[fid]*pow(E, -dot(cj-ci, cj-ci)/b_c_)*pow(E, -dot(nj-ni, nj-ni)/b_s_);
      wa += w_tmp;
      new_ni += w_tmp * nj;
    }
    normal(colon(), i) = new_ni/wa;
  }
}

bool zsw::BilateralNormalFilter::queryFidOneRing(const size_t fid, const jtf::mesh::tri_mesh &trimesh,  std::vector<size_t> &fid_one_ring)
{
  zjucad::matrix::matrix<size_t> vid = trimesh.trimesh_.mesh_(zjucad::matrix::colon(), fid);

  {
    std::pair<size_t, size_t> res = trimesh.ea_->query(vid[0], vid[1]);
    if(res.first == -1 || res.second == -1) {
      // std::cout << vid[0] << " " << vid[1] << std::endl;
      // std::cout << res.first << " " << res.second << std::endl;
      return false;
    }
    if(res.first != fid) { fid_one_ring.push_back(res.first); }
    else { fid_one_ring.push_back(res.second); }
  }
  {
    std::pair<size_t, size_t> res = trimesh.ea_->query(vid[0], vid[2]);
    if(res.first == -1 || res.second == -1) {
      // std::cout << vid[0] << " " << vid[2] << std::endl;
      // std::cout << res.first << " " << res.second << std::endl;
      return false;
    }
    if(res.first != fid) { fid_one_ring.push_back(res.first); }
    else { fid_one_ring.push_back(res.second); }
  }
  {
    std::pair<size_t, size_t> res = trimesh.ea_->query(vid[1], vid[2]);
    if(res.first == -1 || res.second == -1) {
      // std::cout << vid[1] << " " << vid[2] << std::endl;
      // std::cout << res.first << " " << res.second << std::endl;
      return false;
    }
    if(res.first != fid) { fid_one_ring.push_back(res.first); }
    else { fid_one_ring.push_back(res.second); }
  }
  return true;
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

void zsw::BilateralNormalFilter::postProcessing(jtf::mesh::tri_mesh &trimesh)
{
  std::cerr << "Function " << __FUNCTION__ << "in " << __FILE__ << __LINE__  << " haven't implement!!!" << std::endl;
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

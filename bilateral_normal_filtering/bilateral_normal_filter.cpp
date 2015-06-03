#include "bilateral_normal_filter.h"

#include <fstream>
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
  // b_c_ = 0.3; // should dynamiclly caculated as the average distance of 1-ring face center distance
}

void zsw::BilateralNormalFilter::filter(jtf::mesh::tri_mesh &trimesh)
{
  for(size_t i=0; i<st_; ++i) {
    std::cout << "smooth normal step " << i << std::endl;
    filterNormal(trimesh);
  }

  for(size_t i=0; i<ut_; ++i) {
    std::cout << "update vertex step " << i << std::endl;
    updateVertex(trimesh);
  }
  // jtf::mesh::save_obj("/home/wegatron/tmp/tooth_res.obj", trimesh.trimesh_.mesh_, trimesh.trimesh_.node_);
  // writeTriMesh("/home/wegatron/tmp/tooth_debug0.obj", trimesh.trimesh_.mesh_, trimesh.trimesh_.node_, trimesh.face_normal_);
}

void zsw::BilateralNormalFilter::filterNormal(jtf::mesh::tri_mesh &trimesh)
{
  using namespace zjucad::matrix;
  matrixd &node = trimesh.trimesh_.node_;
  matrixst &mesh = trimesh.trimesh_.mesh_;
  matrixd &normal = trimesh.face_normal_;
  matrixd &face_area = trimesh.face_area_;
  matrixd tmp_normal = normal;
  assert(node.size(1) == 3 && mesh.size(1) == 3 && normal.size(2)==mesh.size(2));
#if !ONE_RING_I
  preProcess(trimesh);
#endif
  for(size_t i=0; i<mesh.size(2); ++i) {
    // caculate c_i, n_i
    matrixd ci = ( node(colon(), mesh(0,i)) + node(colon(), mesh(1,i)) + node(colon(), mesh(2,i)) )/3;
    matrixd ni = normal(colon(), i);
    vector<size_t> fid_one_ring;
#if ONE_RING_I
    if(!queryFidOneRingI(i, trimesh, fid_one_ring)) { continue; } // boundary cell
#else
    queryFidOneRingII(i, mesh, fid_one_ring);
#endif
    if(fid_one_ring.size()<3) {
      std::cerr << "[INFO] fid" << i << " one ring: "<< fid_one_ring.size() << std::endl;
    }
    double wa = 0.0;
    matrixd new_ni = zjucad::matrix::zeros(3,1);
    b_c_ = calBc(i, fid_one_ring, mesh, node);
    for(const size_t fid : fid_one_ring) {
      const matrixd cj = ( node(colon(), mesh(0,fid)) + node(colon(), mesh(1,fid)) + node(colon(), mesh(2,fid)) )/3;
      const matrixd nj = normal(colon(), fid);
      const double w_tmp = face_area[fid]*pow(E, -dot(cj-ci, cj-ci)/b_c_)*pow(E, -dot(nj-ni, nj-ni)/b_s_);
      wa += w_tmp;
      new_ni += w_tmp * nj;
    }
    new_ni /= wa;
    if(zjucad::matrix::norm(new_ni) > 1e-8) {
      tmp_normal(colon(), i) = new_ni/norm(new_ni); // normalize
    } else {
      std::cerr << "normal so min" << std::endl;
    }
  }
  normal = tmp_normal;
}

bool zsw::BilateralNormalFilter::queryFidOneRingI(const size_t fid, const jtf::mesh::tri_mesh &trimesh,  std::vector<size_t> &fid_one_ring)
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

void zsw::BilateralNormalFilter::queryFidOneRingII(const size_t fid, const matrixst &mesh, std::vector<size_t> &fid_one_ring)
{
  typedef std::multimap<size_t, size_t>::iterator mit;
  assert(v2f_.size()!=0);
  set<size_t> tmp_one_ring_f;
  for(size_t i=0; i<3; ++i) {
    std::pair<mit, mit> ret = v2f_.equal_range(mesh(i, fid));
    for(mit tmp_it = ret.first; tmp_it!= ret.second; ++tmp_it) {
      tmp_one_ring_f.insert(tmp_it->second);
    }
  }
  tmp_one_ring_f.erase(fid);
  fid_one_ring.resize(tmp_one_ring_f.size());
  std::copy(tmp_one_ring_f.begin(), tmp_one_ring_f.end(), fid_one_ring.begin());
}

void zsw::BilateralNormalFilter::preProcess(const jtf::mesh::tri_mesh &trimesh)
{
  using namespace zjucad::matrix;
  const matrixst &mesh = trimesh.trimesh_.mesh_;
  for(size_t i=0; i<mesh.size(2); ++i) {
    v2f_.insert(std::pair<size_t, size_t>(mesh(0,i),i));
    v2f_.insert(std::pair<size_t, size_t>(mesh(1,i),i));
    v2f_.insert(std::pair<size_t, size_t>(mesh(2,i),i));
  }
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

double zsw::BilateralNormalFilter::calBc(const size_t fid, const std::vector<size_t> &fid_one_ring,
                                         const zjucad::matrix::matrix<size_t> &mesh,
                                         const zjucad::matrix::matrix<double> &node)
{
  using namespace zjucad::matrix;
  matrixd ci = node(colon(), mesh(0, fid)) + node(colon(), mesh(1, fid)) + node(colon(), mesh(2, fid));
  double dis = 0.0;
  for(size_t t_fid : fid_one_ring) {
    matrixd cj = node(colon(), mesh(0, t_fid)) + node(colon(), mesh(1, t_fid)) + node(colon(), mesh(2, t_fid));
    dis += norm(ci-cj);
  }
  dis /= fid_one_ring.size();
  return dis*dis;
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

#include <iostream>
#include <fstream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <zswlib/const_val.h>
#include <zswlib/error_ctrl.h>
#include <zswlib/mesh/vtk.h>
#include "../sampling.h"

using namespace std;

void test_calcLocalCoordinate()
{
  zsw::Sampler sampler;
  Eigen::Matrix<zsw::Scalar,3,3> tri_points;
  tri_points<< 1,0,0,
    0,1,0,
    0,0,1;
  Eigen::Matrix<zsw::Scalar,3,1> translate;
  Eigen::Matrix<zsw::Scalar,3,3> rotate;
  sampler.calcLocalCoordinate(tri_points, translate, rotate);

  Eigen::Matrix<zsw::Scalar,3,3> local_points = tri_points;
  local_points += translate * Eigen::Matrix<zsw::Scalar,1,3>::Ones();
  local_points = rotate * local_points;

  std::cerr << "local_points:\n" << local_points << std::endl;

  std::cerr << "points transform back:\n" << (rotate.inverse() * local_points) - translate * Eigen::Matrix<zsw::Scalar,1,3>::Ones() << std::endl;
}

void test_sampleTriangle()
{
  zsw::Sampler sampler;
  Eigen::Matrix<zsw::Scalar,3,3> tri_points;
  tri_points<< 1,0,0,
    0,1,0,
    0,0,1;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> samples;
  sampler.sampleTriangle(tri_points, 0.1, samples);
  // write obj file
  ofstream ofs("/home/wegatron/tmp.obj", std::ofstream::out);
  ofs << "f 1 2 3" << std::endl;
  for(int i=0; i<3; ++i) {
    ofs << "v " << tri_points.block<3,1>(0,i).transpose() << std::endl;
  }
  for(Eigen::Matrix<zsw::Scalar,3,1> &tmp_sample : samples) {
    ofs << "v " << tmp_sample.transpose() << std::endl;
  }
  ofs.close();
}

void test_sameSide()
{
  Eigen::Matrix<zsw::Scalar,3,1> v0, v1, vr, vt;
  v0 << -1, 0, 0; v1 << 5, 1, 0;
  vr << 0, 0, 0; vt << 0, 0.1667, 0;
  std::cerr << zsw::sameSide2D(v0,v1,vr,vt) << std::endl;
}

void test_projectToLine()
{
  Eigen::Matrix<zsw::Scalar,3,1> a = Eigen::Matrix<zsw::Scalar,3,1>::Random();
  Eigen::Matrix<zsw::Scalar,3,1> b = Eigen::Matrix<zsw::Scalar,3,1>::Random();
  Eigen::Matrix<zsw::Scalar,3,1> c0 = Eigen::Matrix<zsw::Scalar,3,1>::Random();
  Eigen::Matrix<zsw::Scalar,3,1> c1 = c0;
  zsw::Sampler sampler;
  sampler.projectToLine(a, b, c1);
  if(fabs((c1-c0).dot(b-a)) > zsw::const_val::eps) {
    std::cerr << "Error with projectToLine!" << std::endl;
  } else if(((c1-a).cross(b-a)).norm() > zsw::const_val::eps) {
    std::cerr << "Error with projectToLine!" << std::endl;
  }
  std::cerr << "a " << a.transpose() << std::endl;
  std::cerr << "b " << b.transpose() << std::endl;
  std::cerr << "c0 " << c0.transpose() << std::endl;
  std::cerr << "c1 " << c1.transpose() << std::endl;
}

void test_resolvePoint()
{
  {
    Eigen::Matrix<zsw::Scalar,3,3> tri_points;
    tri_points<<0,0,2,
      0,4,2,
      0,0,0;
    Eigen::Matrix<zsw::Scalar,3,1> sample_point;
    sample_point << 2,3,0;
    zsw::Sampler sampler;
    sampler.resolvePoint(tri_points, sample_point);
    std::cerr << sample_point.transpose() << std::endl;
  }

  {
    Eigen::Matrix<zsw::Scalar,3,3> tri_points;
    tri_points<<0,4,2,
      0,0,2,
      0,0,0;
    Eigen::Matrix<zsw::Scalar,3,1> sample_point;
    sample_point << 2,3,0;
    zsw::Sampler sampler;
    sampler.resolvePoint(tri_points, sample_point);
    std::cerr << sample_point.transpose() << std::endl;
  }
}

void test_sampleing()
{
  zsw::mesh::TriMesh tm;
  if(!OpenMesh::IO::read_mesh(tm, "/home/wegatron/workspace/geometry/data/dragon.obj")) {
    std::cerr << "[ERROR] can't read mesh: /home/wegatron/workspace/geometry/data/dragon.obj"<< std::endl;
    return;
  }
  zsw::Sampler sampler;
  std::vector<Eigen::Matrix<zsw::Scalar,3,1>> samples;
  sampler.sampleSigmaDense(tm, 0.03, samples);

  // write obj file
  ofstream ofs("/home/wegatron/tmp.obj", std::ofstream::out);
  for(zsw::mesh::TriMesh::ConstFaceIter f_it=tm.faces_begin(); f_it!=tm.faces_end(); ++f_it) {
    ofs << "f";
    for(zsw::mesh::TriMesh::FaceVertexIter fv_it=tm.fv_iter(*f_it); fv_it.is_valid(); ++fv_it) {
      ofs << " " << fv_it->idx()+1;
    }
    ofs << endl;
  }
  for(zsw::mesh::TriMesh::ConstVertexIter v_it=tm.vertices_begin(); v_it!=tm.vertices_end(); ++v_it) {
    ofs << "v " << (tm.point(*v_it)).transpose() << endl;
  }

  for(Eigen::Matrix<zsw::Scalar,3,1> &tmp_sample : samples) {
    ofs << "v " << tmp_sample.transpose() << std::endl;
  }
  ofs.close();
}

void test_sample_tet()
{
  Eigen::Matrix<zsw::Scalar,3,4> tet_points;
  tet_points << 2,0,0,1,
                          0,4,-1,1,
                          0,0,0,4;
  std::vector<zsw::Point> sample_points;
  zsw::sampleTet(0.4, tet_points, sample_points);
  std::ofstream ofs;
  OPEN_STREAM("/home/wegatron/test_sample_tet.vtk", ofs, std::ofstream::out, return);
  vector<zsw::Scalar> pt_data(12, 0);
  copy(tet_points.data(), tet_points.data()+12, &pt_data[0]);
  for(const zsw::Point &pt : sample_points) {
    pt_data.push_back(pt[0]);
    pt_data.push_back(pt[1]);
    pt_data.push_back(pt[2]);
  }
  size_t tets[4] = {0,1,2,3};
  tet2vtk(ofs, &pt_data[0], 4+sample_points.size(), &tets[0], 1);
}

int main(int argc, char *argv[])
{
  test_sample_tet();
  //test_sampleTriangle();
  //test_sampleing();
  //test_calcLocalCoordinate();
  // test_sameSide();
  // test_projectToLine();
  // test_resolvePoint();
  return 0;
}

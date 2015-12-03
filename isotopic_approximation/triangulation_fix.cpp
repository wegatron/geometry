#include "triangulation_fix.h"

namespace zsw{
  void removeSliverTet(Delaunay &delaunay, const std::vector<zsw::Vertex> &vertices)
  {
    bool flipable = true;
    while(flipable) {
      flipable = false;
      for(Delaunay::Finite_cells_iterator cit=delaunay.finite_cells_begin();
          cit!=delaunay.finite_cells_end(); ++cit) {
        std::pair<size_t,size_t> e[2];
        if(isSliverTet(vertices[cit->vertex(0)->info()].pt_, vertices[cit->vertex(1)->info()].pt_,
                       vertices[cit->vertex(2)->info()].pt_, vertices[cit->vertex(3)->info()].pt_, e[0], e[1])) {
          if(delaunay.flip(cit, e[0].first, e[0].second) || delaunay.flip(cit, e[1].first, e[1].second)) {
            flipable = true;
            break;
          }
        }
      }
    }
  }


  bool isSliverTet(const Eigen::Matrix<zsw::Scalar,3,1> &v0, const Eigen::Matrix<zsw::Scalar,3,1> &v1,
                   const Eigen::Matrix<zsw::Scalar,3,1> &v2, const Eigen::Matrix<zsw::Scalar,3,1> &v3,
                   std::pair<size_t,size_t> &e0, std::pair<size_t,size_t> &e1)
  {
    Eigen::Matrix<zsw::Scalar,3,1> va = v1 - v0;
    Eigen::Matrix<zsw::Scalar,3,1> vb = v2 - v0;
    Eigen::Matrix<zsw::Scalar,3,1> vc = v3 - v0;
    Eigen::Matrix<zsw::Scalar,3,1> vn = va.cross(vb); vn.normalize();
    if(fabs(vc.dot(vn)) > 10*zsw::const_val::eps) { return false; }
    if((va.cross(vc)).dot(va.cross(vb)) < 0) {
      e0 = std::pair<size_t,size_t>(0,1);
      e1 = std::pair<size_t,size_t>(2,3);
    } else if((vb.cross(vc)).dot(vb.cross(va)) < 0) {
      e0 = std::pair<size_t,size_t>(0,2);
      e1 = std::pair<size_t,size_t>(1,3);
    } else {
      e0 = std::pair<size_t,size_t>(0,3);
      e1 = std::pair<size_t,size_t>(1,2);
    }
    return true;
  }

  void haveSliverTet(Delaunay &delaunay, const std::vector<zsw::Vertex> &vertices)
  {
    for(Delaunay::Finite_cells_iterator cit=delaunay.finite_cells_begin();
        cit!=delaunay.finite_cells_end(); ++cit) {
      std::pair<size_t,size_t> e[2];
      if(isSliverTet(vertices[cit->vertex(0)->info()].pt_, vertices[cit->vertex(1)->info()].pt_,
                     vertices[cit->vertex(2)->info()].pt_, vertices[cit->vertex(3)->info()].pt_, e[0], e[1])) {
        std::cerr << "Have sliver Tet:" << cit->vertex(0)->info() << " "
                  << cit->vertex(1)->info() << " " << cit->vertex(2)->info() << " "
                  << cit->vertex(3)->info() << std::endl;
      }
    }
  }

  zsw::Scalar calcTetQuality(const Eigen::Matrix<zsw::Scalar,3,1> &v0,
                             const Eigen::Matrix<zsw::Scalar,3,1> &v1,
                             const Eigen::Matrix<zsw::Scalar,3,1> &v2,
                             const Eigen::Matrix<zsw::Scalar,3,1> &v3)
  {
    static const zsw::Scalar N = 9*sqrt(3)/8.0;
    Eigen::Matrix<zsw::Scalar,3,3> vol_tet_mat;
    vol_tet_mat.block<3,1>(0,0) = v1-v0;
    vol_tet_mat.block<3,1>(0,1) = v2-v0;
    vol_tet_mat.block<3,1>(0,2) = v3-v0;
    zsw::Scalar vol_tet = fabs(vol_tet_mat.determinant()/6.0);
    Eigen::PartialPivLU<Eigen::Matrix<zsw::Scalar,3,3>> pplu;
    pplu.compute(2*vol_tet_mat.transpose());
    zsw::Scalar v_v0 = v0.squaredNorm();
    zsw::Scalar v_v1 = v1.squaredNorm();
    zsw::Scalar v_v2 = v2.squaredNorm();
    zsw::Scalar v_v3 = v3.squaredNorm();
    Eigen::Matrix<zsw::Scalar,3,1> b; b<<v_v1 - v_v0, v_v2-v_v0, v_v3-v_v0;
    Eigen::Matrix<zsw::Scalar,3,1> c = pplu.solve(b);
    zsw::Scalar r = (c-v0).norm();
    return pow(N*vol_tet, 0.333)/r;
  }

}

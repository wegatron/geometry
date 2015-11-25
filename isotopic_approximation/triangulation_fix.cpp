#include "triangulation_fix.h"


namespace zsw{
  void removeSliverTet(Delaunay &delaunay, const std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &vertices)
  {
    for(Delaunay::Finite_cells_iterator cit=delaunay.finite_cells_begin();
        cit!=delaunay.finite_cells_end(); ++cit) {
      std::pair<size_t,size_t> e[2];
      if(isSliverTet(vertices[cit->vertex(0)->info()], vertices[cit->vertex(1)->info()],
                     vertices[cit->vertex(2)->info()], vertices[cit->vertex(3)->info()], e[0], e[1])) {
        if(delaunay.flip(cit, e[0].first, e[0].second)) {
        } else if(delaunay.flip(cit, e[1].first, e[1].second)) {
        } else {
          std::cerr << "cant remove sliver tet:" << cit->vertex(0)->info()
                    << cit->vertex(1)->info() << cit->vertex(2)->info()
                    << cit->vertex(3)->info() << std::endl;
        }
      }
    }
  }
}

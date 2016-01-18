#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_set>
#include <zswlib/error_ctrl.h>
#include "isotopic_approximation.h"
#include "basic_op.h"
#include "bound_sphere.h"
#include "sampling.h"
#include "smoother.h"

using namespace std;

namespace zsw{
  void Approximation::mutuallTessellation()
  {
    const TTds &tds = tw_->getTds();
    std::vector<Vhd> vhds;
    for(auto vit=tds.vertices_begin(); vit!=tds.vertices_end(); ++vit) {
      if(vit->info().pt_type_==zsw::INNER_POINT) { vhds.push_back(vit); }
    }

    for(auto tmp_vhd : vhds) {
      assert(tmp_vhd->info().pt_type_==zsw::INNER_POINT);
      do {
        std::vector<TTds::Edge> edges;
        tds.incident_edges(tmp_vhd, std::back_inserter(edges));
        auto it=find_if(edges.begin(), edges.end(), [](const TTds::Edge &e){
            return e.first->vertex(e.second)->info().pt_type_==zsw::OUTER_POINT ||
            e.first->vertex(e.third)->info().pt_type_==zsw::OUTER_POINT;});
        if(it==edges.end()) { break; }
        TTds::Edge &e = *it;
        Point &pa = e.first->vertex(e.second)->point();
        Point &pb = e.first->vertex(e.third)->point();
        Point pt((pa[0]+pb[0])/2.0, (pa[1]+pb[1])/2.0, (pa[2]+pb[2])/2.0);
        tw_->insertInEdge(e, pt, zsw::ZERO_POINT);
      } while(1);
    }
  }
}

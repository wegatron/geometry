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

  void Approximation::simpBZEdges()
  {
    std::unordered_map<std::string, TTds::Edge> bz_map;
    const TTds &tds=tw_->getTds();
    for(auto eit=tds.edges_begin(); eit!=tds.edges_end(); ++eit) {
      if(!tw_->isBZEdge(*eit)) { continue; }
      std::string key_str=edge2key(*eit);
      bz_map.insert(std::make_pair(key_str, *eit));
    }
    size_t zb_c_step=0;
    size_t zb_total=0;
    while(!bz_map.empty()) {
      TTds::Edge e=bz_map.begin()->second; bz_map.erase(bz_map.begin());
      if(!tw_->isBZEdge(e) && !tw_->isZeroEdge(e)) { continue; }
      if(zb_total++ % 100==0) { std::cout  << "[INFO] all edge tried collapsed " << zb_total << std::endl; }
      if(tryCollapseBZEdge(e, bz_map, true)) {
        if(++zb_c_step%50==0) { std::cout << "[INFO] all edge collapsed " << zb_c_step << std::endl; }
      }
    }
    NZSWLOG("zsw_log") << "zb tried collapse:" << zb_total << std::endl;
    NZSWLOG("zsw_log") << "zb collapsed total:" << zb_c_step << std::endl;
    NZSWLOG("zsw_log") << "zb collapsed suc:" << zb_c_step*1.0/zb_total << std::endl;
  }

  void Approximation::bzEdgeBack(Vhd vhd, std::unordered_map<std::string, TTds::Edge> &edge_map) const
  {
    const TTds &tds=tw_->getTds();
    std::vector<TTds::Edge> edges;
    tds.incident_edges(vhd, std::back_inserter(edges));
    for(const TTds::Edge &e : edges) {
      if(!tw_->isBZEdge(e) && !tw_->isZeroEdge(e)) { continue; }
      std::string key_str=edge2key(e);
      edge_map.insert(std::make_pair(key_str, e));
    }
  }

}

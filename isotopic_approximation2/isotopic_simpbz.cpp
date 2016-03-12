#include <iostream>
#include <fstream>
#include <algorithm>
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
    std::queue<std::pair<TTds::Edge,size_t>> bz_q;
    const TTds &tds=tw_->getTds();
    for(auto eit=tds.edges_begin(); eit!=tds.edges_end(); ++eit) {
      if(tw_->isBZEdge(*eit)) { bz_q.push(std::make_pair(*eit, 0)); }
    }
    tw_->resetVertexLastUpdate();
    size_t zb_step = 0;
    size_t zb_step_suc = 0;
    print_sp_size_= true;
    while(!bz_q.empty()) {
      TTds::Edge e = bz_q.front().first;
      size_t last_update = bz_q.front().second;
      bz_q.pop();
      if(!tw_->isBZEdge(e) && !tw_->isZeroEdge(e)) { continue; }
      if(e.first->vertex(e.second)->info().last_update_ > last_update ||
         e.first->vertex(e.third)->info().last_update_ > last_update) { continue; }
      if(++zb_step % 100==0) { std::cout  << "[INFO] all edge tried collapsed " << zb_step << std::endl; print_sp_size_= true; }
      if(tryCollapseBZEdge(e, bz_q, zb_step_suc, true)) {
        if(++zb_step_suc %50==0) {
          std::cout << "[INFO] all edge collapsed " << zb_step_suc << std::endl;
          writeTetMesh(tmp_outdir_+"simp_bz"+to_string(zb_step_suc)+".vtk", {zsw::ignore_out, zsw::ignore_bbox});
        }
      }
    }
    NZSWLOG("zsw_log") << "zb tried collapse:" << zb_step << std::endl;
    NZSWLOG("zsw_log") << "zb collapsed total:" << zb_step_suc << std::endl;
    NZSWLOG("zsw_log") << "zb collapsed suc:" << zb_step_suc*1.0/zb_step << std::endl;
  }

  void Approximation::bzEdgeBack(Vhd vhd, std::queue<std::pair<TTds::Edge, size_t>> &zb_q, size_t cur_update) const
  {
    const TTds &tds=tw_->getTds();
    std::vector<TTds::Edge> edges;
    tds.incident_edges(vhd, std::back_inserter(edges));
    for(const TTds::Edge &e : edges) {
      if(tw_->isBZEdge(e) || tw_->isZeroEdge(e)) {        zb_q.push(std::make_pair(e, cur_update));      }
    }
  }

}

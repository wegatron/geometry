#define _EXPORTING

#include <queue>
#include <boost/assert.hpp>
#include <zswlib/mesh/mesh_op.h>

void zsw::mesh::rRingVertex(const zsw::mesh::TriMesh &tm, const size_t r, std::vector<zsw::FakeSet<size_t>> &ring)
{
  assert(r>0);
  //BOOST_VERIFY_MSG(r>0, "rRing vertex r should > 0!");
  ring.clear();  ring.resize(tm.n_vertices());
  for(size_t i=0; i<ring.size(); ++i) {
    std::queue<size_t> q; q.push(i);
    size_t r_cnt = 0;
    size_t mark = i;
    std::set<size_t> tmp_ring; tmp_ring.insert(i);
    while(!q.empty() && r_cnt<r) {
      size_t vid = q.front(); q.pop();
      for(zsw::mesh::TriMesh::ConstVertexVertexIter vvit = tm.cvv_begin(zsw::mesh::TriMesh::VertexHandle(vid)); vvit.is_valid(); ++vvit) {
        if(tmp_ring.find(vvit->idx()) == tmp_ring.end()) {
          tmp_ring.insert(vvit->idx());
          q.push(vvit->idx());
        }
      }
      if(vid == mark) {
        ++r_cnt;
        if(!q.empty()) mark = q.back();
      }
    }
    tmp_ring.erase(i);
    ring[i].initFromSet(tmp_ring);
  }
}

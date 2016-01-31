#include "kdtree_warp.h"

namespace zsw
{
  KdTreeWarper::KdTreeWarper() {}

  void KdTreeWarper::buildTree(zsw::Scalar *data, size_t n)
  {
    zsw::Scalar *ptr=data;
    for(size_t i=0; i<n; ++i) {
      kdtree_.insert(KdTreeNode(*ptr, *(ptr+1), *(ptr+2), i));
      ptr+=3;
    }
  }

  std::pair<size_t,zsw::Scalar> KdTreeWarper::queryNearest(zsw::Scalar *pt) const
  {
    KdTreeNode tmp_node(*pt, *(pt+1), *(pt+2), -1);
    std::pair<KdTreeType::const_iterator,double> res = kdtree_.find_nearest(tmp_node);
    std::pair<size_t, zsw::Scalar> ret;
    ret.first=res.first->index_;
    ret.second=res.second;
    return ret;
  }
}

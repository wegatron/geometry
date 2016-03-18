#ifndef KDTREE_WARP_H
#define KDTREE_WARP_H

#include <Eigen/Dense>
#include <kdtree++/kdtree.hpp>
#include <zswlib/config.h>

namespace zsw{

  struct KdTreeNode
  {
    zsw::Scalar xyz_[3];
    size_t index_;
    typedef zsw::Scalar value_type;
    KdTreeNode(zsw::Scalar x, zsw::Scalar y, zsw::Scalar z, size_t index) {
      xyz_[0]=x;
      xyz_[1]=y;
      xyz_[2]=z;
      index_=index;
    }

    zsw::Scalar operator[](size_t n) const
    {
      return xyz_[n];
    }

    zsw::Scalar distance( const KdTreeNode &node)
    {
      zsw::Scalar x = xyz_[0] - node.xyz_[0];
      zsw::Scalar y = xyz_[1] - node.xyz_[1];
      zsw::Scalar z = xyz_[2] - node.xyz_[2];
      return x*x+y*y+z*z;
    }
  };

  typedef KDTree::KDTree<3,KdTreeNode> KdTreeType;

  class KdTreeWarper
  {
  public:
    KdTreeWarper();
    void buildTree(zsw::Scalar *data, size_t n);
    std::pair<size_t,zsw::Scalar> queryNearest(zsw::Scalar *nearest) const;

    template<int col_n>
      void queryNearest(const Eigen::Matrix<zsw::Scalar,3,col_n> &pts,
                        std::vector<size_t> &indices,
                        std::vector<zsw::Scalar> &dist) const
      {
        for(size_t i=0; i<col_n; ++i) {
          KdTreeNode tmp_node(pts(0,i), pts(1,i), pts(2,i), -1);
          std::pair<KdTreeType::const_iterator,double> res = kdtree_.find_nearest(tmp_node);
          indices.push_back(res.first->index_);
          dist.push_back(res.second);
        }
      }
    template<typename _OutPutIterator>
      _OutPutIterator findWithinR(const Eigen::Matrix<zsw::Scalar,3,1> pt, const zsw::Scalar r, _OutPutIterator out) const
      {
        KdTreeNode tmp_node(pt(0), pt(1), pt(2), -1);
        return kdtree_.find_within_range(tmp_node, r, out);
      }
  private:
    KdTreeType kdtree_;
  };
}

#endif /* KDTREE_WARP_H */

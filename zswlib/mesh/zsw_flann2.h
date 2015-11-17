#ifndef ZSW_FLANN2_H
#define ZSW_FLANN2_H

#include <vector>
#include <flann/flann.hpp>
#include <Eigen/Dense>
#include <memory>

namespace zsw {
  template<typename ScalarType, size_t ext_info_size=0>
    class Flann2 {
  public:
  Flann2(ScalarType *data, const size_t points_number) : points(data, points_number, 3+ext_info_size) {
      index.reset(new flann::Index<flann::L2_3D<ScalarType> >(points, flann::KDTreeIndexParams(4)));
      index->buildIndex();
    }

  void queryNearest(Eigen::Matrix<ScalarType, 3, Eigen::Dynamic> &q_points, std::vector<size_t> &indices,
                    std::vector<ScalarType> &dist) {
    std::vector<std::vector<int> > tmp_indices;
    std::vector<std::vector<ScalarType> > tmp_dist;
    flann::Matrix<ScalarType> tmp_q_points(q_points.data(), q_points.cols(), 3);
    index->knnSearch(tmp_q_points, tmp_indices, tmp_dist, 1, flann::SearchParams(128));
    indices.resize(q_points.cols());
    for(size_t i=0; i<tmp_indices.size(); ++i) {
      indices[i] = tmp_indices[i].front();
    }
    dist.resize(tmp_dist.size());
    for(size_t i=0; i<tmp_dist.size(); ++i) {
      dist[i] = tmp_dist[i].front();
    }
  }

  void queryInR(Eigen::Matrix<ScalarType,3+ext_info_size,1> &q_point,
                const ScalarType radius,
                std::vector<int> &indices,
                std::vector<ScalarType> &dist)
  {
    flann::Matrix<ScalarType> tmp_q_point(q_point.data(),1,5);
    std::vector<std::vector<int>> tmp_indices;
    std::vector<std::vector<ScalarType>> tmp_dist;
    index->radiusSearch(tmp_q_point,tmp_indices,tmp_dist,radius,flann::SearchParams(128));
    indices=tmp_indices[0];
    dist=tmp_dist[0];
  }

  private:
  std::shared_ptr<flann::Index<flann::L2_3D<ScalarType> > > index;
  flann::Matrix<ScalarType> points; // row-major order
  };
}

#endif /* ZSW_FLANN2_H */

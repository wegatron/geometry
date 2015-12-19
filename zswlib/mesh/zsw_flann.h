#ifndef ZSW_FLANN_H
#define ZSW_FLANN_H

#include <vector>
#include <memory>
#include <flann/flann.hpp>
#include <Eigen/Dense>

namespace zsw {
  template<typename ScalarType>
    class Flann {
  public:
    Flann(ScalarType *data, const size_t points_number) : points(data, points_number, 3) {
      index.reset(new flann::Index<flann::L2<ScalarType> >(points, flann::KDTreeIndexParams(4)));
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
  private:
    std::shared_ptr<flann::Index<flann::L2<ScalarType> > > index;
    flann::Matrix<ScalarType> points; // row-major order
  };
}

#endif /* ZSW_FLANN_H */

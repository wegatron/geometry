#ifndef ZSW_FLANN_H
#define ZSW_FLANN_H

#include <vector>
#include <memory>
#include <flann/flann.hpp>
#include <Eigen/Dense>

#define MAX_ITER 10000

namespace zsw {
  template<typename ScalarType>
    class Flann {
  public:
  Flann(ScalarType *data, const size_t points_number) : points(data, points_number, 3) {
      index.reset(new flann::Index<flann::L2<ScalarType> >(points, flann::KDTreeIndexParams(1)));
      index->buildIndex();
    }

    template<int col_n>
      void queryNearest(const Eigen::Matrix<ScalarType, 3, col_n> &q_points, std::vector<size_t> &indices,
                        std::vector<ScalarType> &dist) {
      std::vector<std::vector<int> > tmp_indices;
      std::vector<std::vector<ScalarType> > tmp_dist;
      const flann::Matrix<ScalarType> tmp_q_points(const_cast<ScalarType*>(q_points.data()), q_points.cols(), 3);
      index->knnSearch(tmp_q_points, tmp_indices, tmp_dist, 1, flann::SearchParams(MAX_ITER));
      indices.resize(q_points.cols());
      for(size_t i=0; i<tmp_indices.size(); ++i) {
        indices[i] = tmp_indices[i].front();
      }
      dist.resize(tmp_dist.size());
      for(size_t i=0; i<tmp_dist.size(); ++i) {
        dist[i] = tmp_dist[i].front();
      }
    }

    void queryNearest(const std::vector<Eigen::Matrix<ScalarType, 3, 1>> &q_points, std::vector<size_t> &indices,
                      std::vector<ScalarType> &dist) {
      std::vector<std::vector<int> > tmp_indices;
      std::vector<std::vector<ScalarType> > tmp_dist;
      const flann::Matrix<ScalarType> tmp_q_points(const_cast<ScalarType*>(q_points[0].data()), q_points.size(), 3);
      index->knnSearch(tmp_q_points, tmp_indices, tmp_dist, 1, flann::SearchParams(MAX_ITER));
      indices.resize(q_points.size());
      for(size_t i=0; i<tmp_indices.size(); ++i) {
        indices[i] = tmp_indices[i].front();
      }
      dist.resize(tmp_dist.size());
      for(size_t i=0; i<tmp_dist.size(); ++i) {
        dist[i] = tmp_dist[i].front();
      }
    }

    void queryKnn(const std::vector<Eigen::Matrix<ScalarType,3,1>> &q_points,
                  std::vector<std::vector<size_t>> &indices,
                  std::vector<std::vector<zsw::Scalar>> &dists,
                  size_t count)
    {
      const flann::Matrix<ScalarType> tmp_q_points(const_cast<ScalarType*>(q_points[0].data()), q_points.size(), 3);
      //specifies the maximum leafs to visit, CHECKS UNLIMITED
      index->knnSearch(tmp_q_points, indices, dists, count, flann::SearchParams(MAX_ITER));
    }

    void addPoints(ScalarType *data, const size_t points_number) {
      flann::Matrix<ScalarType> n_pts(data, points_number,3);
      index->addPoints(n_pts);
    }
  private:
    std::shared_ptr<flann::Index<flann::L2<ScalarType> > > index;
    flann::Matrix<ScalarType> points; // row-major order
  };
}

#endif /* ZSW_FLANN_H */

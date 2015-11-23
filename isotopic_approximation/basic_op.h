#ifndef BASIC_OP_H
#define BASIC_OP_H

#include <Eigen/Dense>
#include <zswlib/config.h>


namespace zsw{

  template<typename SCALAR, int DIMENSION, template<typename, typename...> class Container>
    void calcBBOX(const Container<Eigen::Matrix<SCALAR, DIMENSION, 1>> &vertices,
                  Eigen::Matrix<SCALAR,DIMENSION,2> &bbox)
    {
      assert(vertices.size()>0);
      bbox.block(0,0,3,1)=vertices[0];
      bbox.block(0,1,3,1)=vertices[0];
      for(const Eigen::Matrix<SCALAR,3,1> &vertex : vertices) {
        for(size_t di=0; di<DIMENSION; ++di) {
          if(vertex[di] < bbox(di,0)) { bbox(di,0) = vertex[di]; }
          else if(vertex[di] > bbox(di,1)) { bbox(di,1) = vertex[di]; }
        }
      }
    }

  zsw::Scalar calcPoint2TriSquaredDis(const Eigen::Matrix<zsw::Scalar,3,1> &pt,
                               const Eigen::Matrix<zsw::Scalar,3,1> &v0,
                               const Eigen::Matrix<zsw::Scalar,3,1> &v1,
                               const Eigen::Matrix<zsw::Scalar,3,1> &v2);

}

#endif /* BASIC_OP_H */

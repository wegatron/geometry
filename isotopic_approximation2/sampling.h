#ifndef SAMPLING_H
#define SAMPLING_H

#include <vector>
#include <zswlib/mesh/mesh_type.h>
#include <zswlib/const_val.h>
#include "basic_data_structure.h"

namespace zsw
{
  void sampleTriangle(const Eigen::Matrix<zsw::Scalar, 3, 3> &tri_points, const zsw::Scalar r,
                      std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &samples);


  /// \brief uniform sample a ref triangle with r dense, and mapping it into another triangle.
  ///
  /// uniform sample a ref triangle with r dense, maping it into sample_tri_pts.
  ///
  /// \param
  ///
  void sampleTriangleRefTriangle(const Eigen::Matrix<zsw::Scalar, 3, 3> &sample_tri_pts,
                                 const Eigen::Matrix<zsw::Scalar, 3, 3> &ref_tri_pts, const zsw::Scalar r,
                                 std::vector<Eigen::Matrix<zsw::Scalar,3,1>> &samples);

  template< template<typename, typename...> class Container>
    void sampleTet(const Eigen::Matrix<zsw::Scalar,3,1> &v0,
                   const Eigen::Matrix<zsw::Scalar,3,1> &v1,
                   const Eigen::Matrix<zsw::Scalar,3,1> &v2,
                   const Eigen::Matrix<zsw::Scalar,3,1> &v3,
                   const zsw::Scalar r,
                   Container<Eigen::Matrix<zsw::Scalar,3,1>> &samples)
    {
      Eigen::Matrix<zsw::Scalar,3,1> n=v1-v0;
      Eigen::Matrix<zsw::Scalar,3,1> m=v2-v0;
      Eigen::Matrix<zsw::Scalar,3,1> p=v3-v0;
      zsw::Scalar step_n=r/n.norm();
      zsw::Scalar step_m=r/m.norm();
      zsw::Scalar step_p=r/p.norm();

      bool sn_flg=true;
      for(zsw::Scalar sn=0; sn_flg; sn+=step_n) {
        if(sn>1) { sn=1; sn_flg=false; }
        const Eigen::Matrix<zsw::Scalar,3,1> cur_n=v0+sn*n;
        const zsw::Scalar max_sm=1-sn-zsw::const_val::eps;

        bool sm_flg=true;
        for(zsw::Scalar sm=0; sm_flg; sm+=step_m) {
          if(sm>max_sm) { sm=max_sm; sm_flg=false; }
          const Eigen::Matrix<zsw::Scalar,3,1> cur_m=cur_n+sm*m;
          const zsw::Scalar max_sp=max_sm-sm;

          bool sp_flg=true;
          for(zsw::Scalar sp=0; sp_flg; sp+=step_p) {
            if(sp>max_sp) { sp=max_sp; sp_flg=false; }
            samples.push_back(cur_m+sp*p);
          }
        }
      }
    }
}

#endif /* SAMPLING_H */

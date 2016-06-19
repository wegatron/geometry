#ifndef BILATER_NORMAL_FILTER_H
#define BILATER_NORMAL_FILTER_H

#include <vector>
#include <set>
#include <list>
#include <Eigen/Sparse>
#include "config.h"
#include "mesh_type.h"
#include "data_type.h"

namespace zsw
{
  struct BNFParam{
    enum BNF_TYPE {
      BASIC_FE_RING, // 标准双边滤波 边邻域
      BASIC_FV_RING, // 标准双边滤波 顶点邻域
      EXTENDED_RADIUS_VRING, // 扩展的双边滤波 顶点邻域
      EXTENDED_RADIUS_ERING // 扩展的双边滤波 边邻域
    };

    BNF_TYPE bnf_type_;
    size_t rn_; // r ring
    size_t st_; // filter normal iter times
    size_t ut_; // update vertices iter times
    zsw::Scalar bs_; // normal的高斯函数参数 e^(-bs*(n1.dot(n2))^2)
    zsw::Scalar normal_threshold_; // (r_min,r_max)之间判断面片是否需要作为参考的阈值
    zsw::Scalar smooth_threshold_; // 需要smooth锐利的边的法向夹角阈值
    size_t r_min_; // 选取的最小拓扑半径
    size_t r_max_; // 选取的最大拓扑半径
    bool need_flipedge_; // 是否需要翻边
    bool need_smooth_; // 是否需要laplace smooth锐利的边
  };

  class ZSW_API BilateralNormalFilter //final
  {
  public:
    BilateralNormalFilter(std::shared_ptr<zsw::BNFParam> param);
    void setParam(std::shared_ptr<zsw::BNFParam> param) {
      param_=param;
    }
    void run(std::shared_ptr<zsw::mesh::TriMesh> tm);
    ~BilateralNormalFilter();
  private:
    /***
     * \brief 预计算，计算每一个面片的参考邻域的面片，以及参考邻域中的每一个面片的weight值.
     */
    void preCompute();

    /***
     * \brief 双边滤波更新法向
     **/
    void filterNormal();

    /***
     * \breif 根据更新完的法向更新顶点位置
     **/
    void updateVertex();

    /***
     * \brief 计算一个面片的 高斯函数参数(该面片的所有参考面片到该面片的平均距离的平方) e^(-(1.0/bc)*(ci-cj)^2)
     * \param fid 面片的id
     * \return bc
     **/
    zsw::Scalar calBc(const size_t fid);

    /***
     * \brief 对锐利的边（法向点乘积小于smooth_threshold)的面片上的点做一个局部的laplace smooth
     **/
    void smooth();
    std::shared_ptr<zsw::BNFParam> param_;
    std::shared_ptr<zsw::mesh::TriMesh> tm_;
    // data from precompute
    Eigen::Matrix<zsw::Scalar,Eigen::Dynamic,1> fa_; // face area
    Eigen::Matrix<zsw::Scalar,3,Eigen::Dynamic> fc_; // face center
    Eigen::SparseMatrix<zsw::Scalar> weight_; // weight sparse matrix
    std::vector<zsw::FakeSet<size_t>> ref_face_set_; // 每一个面片的参考面片id
  };

  /***
   * \brief 计算一个面片的参考面片（以面片边的一环邻域）使用如下策略:
   * r_min 最小滤波拓扑半径, if r_{ij}<r_min then j \in Ni
   * r_max 最大滤波拓扑半径, if r_{ij}>r_max then j \notin Ni
   * if r_max>r_{ij}>r_{min} then 如果normal(i).dot(normal(j))>normal_threshold 则 j \in Ni
   * \param trimesh 输入mesh
   * \param normal_threshold 法向阈值
   * \param r_min 最小半径
   * \param r_max 最大半径
   * \param fc 每个面片的中心点位置
   * \return ring 输出每个面片的参考面片
   **/
  void dynamicRringByEdge(const zsw::mesh::TriMesh &trimesh,
                    const zsw::Scalar normal_threshold,
                    const size_t r_min, const size_t r_max,
                    const Eigen::Matrix<zsw::Scalar,3,Eigen::Dynamic> &fc,
                    std::vector<zsw::FakeSet<size_t>> &ring);

  /***
   * \brief 计算一个面片的参考面片（以面片顶点的一环邻域）使用如下策略:
   * r_min 最小滤波拓扑半径, if r_{ij}<r_min then j \in Ni
   * r_max 最大滤波拓扑半径, if r_{ij}>r_max then j \notin Ni
   * if r_max>r_{ij}>r_{min} then 如果normal(i).dot(normal(j))>normal_threshold 则 j \in Ni
   * \param trimesh 输入mesh
   * \param normal_threshold 法向阈值
   * \param r_min 最小半径
   * \param r_max 最大半径
   * \param fc 每个面片的中心点位置
   * \return ring 输出每个面片的参考面片
   **/
  void dynamicRringByVertex(const zsw::mesh::TriMesh &trimesh,
                    const zsw::Scalar normal_threshold,
                    const size_t r_min, const size_t r_max,
                    const Eigen::Matrix<zsw::Scalar,3,Eigen::Dynamic> &fc,
                    std::vector<zsw::FakeSet<size_t>> &ring);


  /***
   * \brief 计算一个mesh的每个面片的r环邻域的面片。邻域的定义为共享边的面片。
   * \param trimesh 输入mesh
   * \param rn rn环邻域
   * \param ring 输出
   **/
  void rRingFacesByEdge(const zsw::mesh::TriMesh &trimesh, const size_t rn, std::vector<zsw::FakeSet<size_t>> &ring);

  /***
   * \brief 计算一个mesh的每个面片的r环邻域的面片。邻域的定义为共享顶点的面片。
   * \param trimesh 输入mesh
   * \param rn rn环邻域
   * \param ring 输出
   **/
  void rRingFacesByVertex(const zsw::mesh::TriMesh &trimesh, const size_t rn, std::vector<zsw::FakeSet<size_t>> &ring);
}
#endif /* BILATER_NORMAL_FILTER_H */

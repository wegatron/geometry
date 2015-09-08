/***************************************************
 * 使用翻边的方式，使特征边变得更加尖锐
 * \note 实现论文Optimizing 3D Triangulations Using Discrete Curvature Analysis。并添加了一些优化
 * \author 张健
 ***************************************************/

#ifndef SN3D_GRAPHICS_EDGE_FLIP_MEAN_H
#define SN3D_GRAPHICS_EDGE_FLIP_MEAN_H

#include <zswlib/config.h>
#include <zswlib/mesh/mesh_type.h>

namespace Sn3DGraphics
{
  /** \brief  使用翻边的方式，使特征边变得更加尖锐
   *  \note  实现论文Optimizing 3D Triangulations Using Discrete Curvature Analysis。并添加了一些优化
   */
  class ZSW_API EdgeFlipMean
  {
  public:
    EdgeFlipMean();
    void run(zsw::mesh::TriMesh& mesh);

    ///是否使用曲率来指导翻边
    void set_is_use_curvature(bool f) { _isUseCurv = f; }
    ///得到是否使用曲率来指导翻边
    bool get_is_use_curvature() const { return _isUseCurv; }

  protected:
    /** \brief 计算vHandle的主曲率方向和最大，最小曲率值
     *  \note  maxk为最大主曲率，mink为最小主曲率。direction表示(最大主曲率方向，最小主曲率方向，法相方向)
     */
    void compute_principal_curvature_direction(zsw::mesh::TriMesh& mesh, zsw::mesh::TriMesh::VertexHandle vHandle,
                                               Eigen::Matrix3d& direction, float& maxk, float& mink);
    float compute_face_area(zsw::mesh::TriMesh& mesh, zsw::mesh::TriMesh::FaceHandle fHandle);
    ///计算hHandle受曲率影响后的向量
    void compute_halfedge_vec_changeby_curvature(zsw::mesh::TriMesh& mesh, zsw::mesh::TriMesh::HalfedgeHandle hHandle,
                                                 zsw::mesh::TriMesh::Normal& vec);
    ///计算hHandle受曲率影响后的对应面片的法向量
    void compute_face_normal_changeby_curvature(zsw::mesh::TriMesh& mesh, zsw::mesh::TriMesh::HalfedgeHandle hHandle,
                                                zsw::mesh::TriMesh::Normal& normal);
    ///eHandle翻边后，更新网格属性
    void update_property_for_flip(zsw::mesh::TriMesh& mesh, zsw::mesh::TriMesh::EdgeHandle eHandle);
    ///ehandle翻边后，返回受影响的其他边
    void get_affect_edge(zsw::mesh::TriMesh& mesh, zsw::mesh::TriMesh::EdgeHandle ehandle, std::vector<zsw::mesh::TriMesh::EdgeHandle>& affectEdge);
    ///计算ehandle翻边后的代价。代价越大越好
    double compute_flip_cost(zsw::mesh::TriMesh& mesh, zsw::mesh::TriMesh::EdgeHandle  ehandle);
    ///计算每个半边的cost
    float compute_cost_about_halfedge(zsw::mesh::TriMesh& mesh, zsw::mesh::TriMesh::HalfedgeHandle hHandle);


    ///使用曲率主方向优化翻边
    float _midCurv;
    float _maxRatio;
    bool _isUseCurv;

    ///记录最大主曲率的方向
    OpenMesh::VPropHandleT<zsw::mesh::TriMesh::Normal> VPropMaxVec;
    ///记录沿最大主曲率放大的比例
    OpenMesh::VPropHandleT<float> VPropRatio;

    ///根据点的曲率，半边改变后的向量
    OpenMesh::HPropHandleT<zsw::mesh::TriMesh::Normal> HPropVec;
    ///记录半边对应face的法向量
    OpenMesh::HPropHandleT<zsw::mesh::TriMesh::Normal> HPropNormal;

  };
}

#endif //SN3D_GRAPHICS_EDGE_FLIP_MEAN_H

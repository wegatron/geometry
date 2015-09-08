#include <zswlib/mesh/EdgeFlipMean.h>

#include <set>

using namespace std;
using namespace Eigen;
using namespace OpenMesh;

#define _EXPORTING

namespace Sn3DGraphics
{
  struct EdgeFlipMeanCost
  {
    int id;
    double cost;
    EdgeFlipMeanCost(int i, double c) : id(i), cost(c) {}
    bool operator <(const EdgeFlipMeanCost& e) const
    {
      if(cost == e.cost) return id > e.id;
      return cost > e.cost;
    }
    bool operator >(const EdgeFlipMeanCost& e) const
    {
      if(cost == e.cost) return id < e.id;
      return cost < e.cost;
    }
  };

  Sn3DGraphics::EdgeFlipMean::EdgeFlipMean() : _midCurv(0.25), _maxRatio(5), _isUseCurv(true) {}
  void EdgeFlipMean::run(zsw::mesh::TriMesh& mesh)
  {
    //计算每个点的主曲率方向和放大比例
    if(_isUseCurv)
      {
        mesh.add_property(VPropMaxVec);
        mesh.add_property(VPropRatio);
        Eigen::Matrix3d m;
        float maxk, mink;
        for(zsw::mesh::TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
          {
            compute_principal_curvature_direction(mesh, v_it, m, maxk, mink);
            mesh.property(VPropMaxVec, v_it) = m.block(0, 0, 3, 1);
            mesh.property(VPropMaxVec, v_it).normalize();

            float curvLen = fabs(maxk) - fabs(mink); 
            float ratio = curvLen / _midCurv;
            if(ratio > _maxRatio) ratio = _maxRatio;
            mesh.property(VPropRatio, v_it) = ratio;
          }
      }

    //根据每个点的曲率信息，形变由它发出的半边向量
    mesh.add_property(HPropVec);
    for(zsw::mesh::TriMesh::HalfedgeIter h_it = mesh.halfedges_begin(); h_it != mesh.halfedges_end(); ++h_it)
      {
        compute_halfedge_vec_changeby_curvature(mesh, h_it, mesh.property(HPropVec, h_it));
      }
    //计算每个面片对应的法向量
    mesh.add_property(HPropNormal);
    for(zsw::mesh::TriMesh::HalfedgeIter h_it = mesh.halfedges_begin(); h_it != mesh.halfedges_end(); ++h_it)
      {
        compute_face_normal_changeby_curvature(mesh, h_it, mesh.property(HPropNormal, h_it));
      }


    //初始化flipCost，计算所有边的翻转代价
    int num = mesh.n_edges();
    vector<double> cost(mesh.n_edges(), -1);
    set<EdgeFlipMeanCost> flipCost;
    for(zsw::mesh::TriMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
      {
        if(mesh.is_flip_ok(e_it))
          {
            double t = compute_flip_cost(mesh, e_it);
            if(t > 0) {
              flipCost.insert(EdgeFlipMeanCost(e_it.handle().idx(), t));
              cost[e_it.handle().idx()] = t;
            }
          }
      }

    //迭代，依次翻转所有满足条件的边
    while(!flipCost.empty())
      {
        //flip edge
        EdgeFlipMeanCost temp = *(flipCost.begin());
        flipCost.erase(flipCost.begin());
        cost[temp.id] = -1;

        mesh.flip(zsw::mesh::TriMesh::EdgeHandle(temp.id));
        //跟新HPropVec和HPropNormal属性
        update_property_for_flip(mesh, zsw::mesh::TriMesh::EdgeHandle(temp.id));
        //update affect edge flip cost
        vector<zsw::mesh::TriMesh::EdgeHandle> affectEdge;
        get_affect_edge(mesh, zsw::mesh::TriMesh::EdgeHandle(temp.id), affectEdge);
        for(size_t i = 0; i < affectEdge.size(); ++i)
          {
            int id = affectEdge[i].idx();
            if(cost[id] != -1)
              {
                flipCost.erase(EdgeFlipMeanCost(id, cost[id]));
                cost[id] = -1;
              }
            if(mesh.is_flip_ok(affectEdge[i]))
              {
                double t = compute_flip_cost(mesh, affectEdge[i]);
                if(t > 0)
                  {
                    flipCost.insert(EdgeFlipMeanCost(id, t));
                    cost[id] = t;
                  }
              }
          }
      }

    mesh.remove_property(HPropVec);
    mesh.remove_property(HPropNormal);
    if(_isUseCurv)
      {
        mesh.remove_property(VPropMaxVec);
        mesh.remove_property(VPropRatio);
      }
  }

  void EdgeFlipMean::compute_principal_curvature_direction(zsw::mesh::TriMesh& mesh, zsw::mesh::TriMesh::VertexHandle vHandle,
                                                           Eigen::Matrix3d& direction, float& maxk, float& mink)
  {

    float area = 0;
    Eigen::Matrix3d tensor;
    tensor.setZero();
    zsw::mesh::TriMesh::Point center = mesh.point(vHandle);
    for(zsw::mesh::TriMesh::VertexOHalfedgeIter vh_it = mesh.voh_iter(vHandle); vh_it; ++vh_it)
      {
        //compute tensor
        zsw::mesh::TriMesh::Point v = mesh.point(mesh.to_vertex_handle(vh_it));
        v = v - center;
        float len = v.norm();
        if(len > 1e-12)
          {
            v /= len;
            tensor += mesh.calc_dihedral_angle(vh_it)*len*v*v.transpose();
          }
        //compute area
        zsw::mesh::TriMesh::FaceHandle fHandle = mesh.face_handle(vh_it.handle());
        if(fHandle.is_valid())
          {
            float t = compute_face_area(mesh, fHandle);
            area += t;
          }
      }
    //compute eigenVector;
    if(area < 1e-12)
      {
        direction.setIdentity();
        maxk = 1;
        mink = 1;
        return;
      }
    tensor /= area;
    Eigen::SelfAdjointEigenSolver<Matrix3d> eigensolver(tensor);
    if(eigensolver.info() != Eigen::Success)
      {
        direction.setIdentity();
        maxk = 1;
        mink = 1;
        return;
      }
    Eigen::Vector3d eigenValue = eigensolver.eigenvalues();
    Eigen::Matrix3d eigenVector = eigensolver.eigenvectors();
    //choose maxVector minVector normal
    int nid;
    double min;
    vector<int> id(3);
    for(int i = 0; i < 3; ++i) id[i] = i;
    for(int i = 0; i < 2; ++i)
      {
        min = eigenValue[i];
        nid = i;
        for(int j = i+1; j < 3; ++j)
          {
            if(fabs(eigenValue[j]) > fabs(min))
              {
                min = eigenValue[j];
                nid = j;
              }
            if(i != nid)
              {
                float t = eigenValue[i];
                eigenValue[i] = eigenValue[nid];
                eigenValue[nid] = t;
                int tt = id[i];
                id[i] = id[nid];
                id[nid] = tt;
              }
          }
      }
    direction.block(0, 0, 3, 1) = eigenVector.block(0, id[1], 3, 1);
    direction.block(0, 1, 3, 1) = eigenVector.block(0, id[0], 3, 1);
    direction.block(0, 2, 3, 1) = eigenVector.block(0, id[2], 3, 1);
    maxk = eigenValue[0];
    mink = eigenValue[1];
  }


  float EdgeFlipMean::compute_face_area(zsw::mesh::TriMesh& mesh, zsw::mesh::TriMesh::FaceHandle fHandle)
  {
    zsw::mesh::TriMesh::FaceVertexIter fv_it = mesh.fv_iter(fHandle);
    zsw::mesh::TriMesh::FaceVertexIter fv_it1 = fv_it++;
    zsw::mesh::TriMesh::FaceVertexIter fv_it2 = fv_it++;
    return (mesh.point(fv_it1) - mesh.point(fv_it)).cross(mesh.point(fv_it2) - mesh.point(fv_it)).norm() / 2;
  }

  void EdgeFlipMean::compute_halfedge_vec_changeby_curvature(zsw::mesh::TriMesh& mesh, zsw::mesh::TriMesh::HalfedgeHandle hHandle,
                                                             zsw::mesh::TriMesh::Normal& vec)
  {
    zsw::mesh::TriMesh::VertexHandle center = mesh.from_vertex_handle(hHandle);
    zsw::mesh::TriMesh::VertexHandle v = mesh.to_vertex_handle(hHandle);
    vec = mesh.point(v) - mesh.point(center);
    if(!_isUseCurv) return;
    float len = vec.dot(mesh.property(VPropMaxVec, center));
    vec += len*mesh.property(VPropRatio, center) * mesh.property(VPropMaxVec, center);
  }

  void EdgeFlipMean::compute_face_normal_changeby_curvature(zsw::mesh::TriMesh& mesh, zsw::mesh::TriMesh::HalfedgeHandle hHandle,
                                                            zsw::mesh::TriMesh::Normal& normal)
  {
    if(mesh.face_handle(hHandle).is_valid() == false) return;
    zsw::mesh::TriMesh::HalfedgeHandle hHandle1 = mesh.opposite_halfedge_handle(mesh.prev_halfedge_handle(hHandle));
    normal = mesh.property(HPropVec, hHandle).cross(mesh.property(HPropVec, hHandle1));
    float len = normal.norm();
    if(len > 1e-12) normal /= len;
    else normal.setZero();
  }

  void EdgeFlipMean::update_property_for_flip(zsw::mesh::TriMesh& mesh, zsw::mesh::TriMesh::EdgeHandle eHandle)
  {
    zsw::mesh::TriMesh::HalfedgeHandle hHandle1 = mesh.halfedge_handle(eHandle, 0);
    zsw::mesh::TriMesh::HalfedgeHandle hHandle2 = mesh.halfedge_handle(eHandle, 1);
    compute_halfedge_vec_changeby_curvature(mesh, hHandle1, mesh.property(HPropVec, hHandle1));
    compute_halfedge_vec_changeby_curvature(mesh, hHandle2, mesh.property(HPropVec, hHandle2));

    zsw::mesh::TriMesh::HalfedgeHandle hHandles[6];
    hHandles[0] = hHandle1;
    hHandles[1] = hHandle2;
    hHandles[2] = mesh.next_halfedge_handle(hHandle1);
    hHandles[3] = mesh.prev_halfedge_handle(hHandle1);
    hHandles[4] = mesh.next_halfedge_handle(hHandle2);
    hHandles[5] = mesh.prev_halfedge_handle(hHandle2);

    for(int i = 0; i < 6; ++i)
      {
        compute_face_normal_changeby_curvature(mesh, hHandles[i], mesh.property(HPropNormal, hHandles[i]));
      }
  }


  void EdgeFlipMean::get_affect_edge(zsw::mesh::TriMesh& mesh, zsw::mesh::TriMesh::EdgeHandle ehandle, vector<zsw::mesh::TriMesh::EdgeHandle>& affectEdge)
  {
    vector<zsw::mesh::TriMesh::HalfedgeHandle> temp;
    temp.reserve(12);
    //寻找被影响的一环邻域
    zsw::mesh::TriMesh::HalfedgeHandle hHandle1 = mesh.halfedge_handle(ehandle, 0);
    if(mesh.face_handle(hHandle1).is_valid())
      {
        temp.push_back(mesh.next_halfedge_handle(hHandle1));
        temp.push_back(mesh.prev_halfedge_handle(hHandle1));
      }
    zsw::mesh::TriMesh::HalfedgeHandle hHandle2 = mesh.halfedge_handle(ehandle, 1);
    if(mesh.face_handle(hHandle2).is_valid())
      {
        temp.push_back(mesh.next_halfedge_handle(hHandle2));
        temp.push_back(mesh.prev_halfedge_handle(hHandle2));
      }
    //寻找被影响的两环邻域
    size_t  num = temp.size();
    for(size_t i = 0; i < num; ++i)
      {
        zsw::mesh::TriMesh::HalfedgeHandle hHandle = mesh.opposite_halfedge_handle(temp[i]);
        if(mesh.face_handle(hHandle).is_valid())
          {
            temp.push_back(mesh.next_halfedge_handle(hHandle));
            temp.push_back(mesh.prev_halfedge_handle(hHandle));
          }
      }
    //得到最终结果
    affectEdge.resize(temp.size());
    for(size_t i = 0; i < temp.size(); ++i)
      {
        affectEdge[i] = mesh.edge_handle(temp[i]);
      }
  }

  double EdgeFlipMean::compute_flip_cost(zsw::mesh::TriMesh& mesh, zsw::mesh::TriMesh::EdgeHandle ehandle)
  {
    //寻找需要计算的半边
    vector<zsw::mesh::TriMesh::HalfedgeHandle> hHandles;
    hHandles.reserve(10);
    hHandles.push_back(mesh.halfedge_handle(ehandle, 0));
    hHandles.push_back(mesh.halfedge_handle(ehandle, 1));
    for(size_t i = 0; i < 2; ++i)
      {
        zsw::mesh::TriMesh::HalfedgeHandle hHandle = mesh.next_halfedge_handle(hHandles[i]);
        hHandles.push_back(hHandle);
        hHandles.push_back(mesh.opposite_halfedge_handle(hHandle));
        hHandle = mesh.next_halfedge_handle(hHandle);
        hHandles.push_back(hHandle);
        hHandles.push_back(mesh.opposite_halfedge_handle(hHandle));
      }
    //计算翻边前cost
    float beforeCost = 0;
    for(size_t i = 0; i < hHandles.size(); ++i)
      {
        beforeCost += compute_cost_about_halfedge(mesh, hHandles[i]);
      }

    mesh.flip(ehandle);
    update_property_for_flip(mesh, ehandle);

    //计算翻边后的Cost
    float afterCost = 0;
    for(size_t i = 0; i < hHandles.size(); ++i)
      {
        afterCost += compute_cost_about_halfedge(mesh, hHandles[i]);
      }

    mesh.flip(ehandle);
    update_property_for_flip(mesh, ehandle);

    return (beforeCost - afterCost)/4;
  }

  float EdgeFlipMean::compute_cost_about_halfedge(zsw::mesh::TriMesh& mesh, zsw::mesh::TriMesh::HalfedgeHandle hHandle)
  {
    zsw::mesh::TriMesh::HalfedgeHandle th = mesh.opposite_halfedge_handle(hHandle);
    if(mesh.face_handle(th).is_valid() == false) return 0;
    zsw::mesh::TriMesh::HalfedgeHandle hHandle1 = mesh.next_halfedge_handle(th);
    return mesh.property(HPropVec, hHandle).norm() 
      * acos(mesh.property(HPropNormal, hHandle).dot(mesh.property(HPropNormal, hHandle1)));
  }

}//namespace

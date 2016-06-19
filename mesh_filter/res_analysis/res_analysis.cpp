#include "res_analysis.h"

#include <fstream>
#include <iostream>
#include "zsw_flann.h"
#include "vtk.h"

#define _EXPORTING

#define DEBUG 0

namespace zsw
{
  namespace mesh{
    void  diffMeshVertex(jtf::mesh2::meshes &mesh_a, jtf::mesh2::meshes &mesh_b,
                         const std::string &output_vtk)
    {
      std::cout << "-----------" << std::endl;
      assert(mesh_a.mesh_.cols() == mesh_b.mesh_.cols());
      assert(mesh_a.node_.cols() == mesh_b.node_.cols());
      zsw::Flann<zsw::Scalar> flann(mesh_b.node_.data(), mesh_b.node_.cols());
      std::vector<size_t> indices;
      std::vector<zsw::Scalar> dist;
      flann.queryNearest(mesh_a.node_, indices, dist);
      for(size_t i=0; i<dist.size(); ++i) { dist[i]=sqrt(dist[i]); }
#if DEBUG
      std::cout << "29435 --> " << indices[29435] << " with dis:" << dist[29435] << std::endl;
#endif
      std::ofstream ofs(output_vtk);
      if(!ofs) {
        throw std::logic_error("can't open file:"+output_vtk+" for write!");
      }
      tri2vtk(ofs, mesh_a.node_.data(), mesh_a.node_.cols(),
              mesh_a.mesh_.data(), mesh_a.mesh_.cols());
      point_data(ofs, dist.begin(), dist.size(), "diff_dist");
      zsw::Scalar max_dis = 0;
      for(size_t i=0; i<dist.size(); ++i) {
        max_dis = std::max(dist[i], max_dis);
      }
      std::cout << "max_dis:" << max_dis << std::endl;;
    }

    void diffMeshNormal(const jtf::mesh2::meshes &mesh_a, const jtf::mesh2::meshes &mesh_b,
                        const std::string &output_vtk)
    {
      std::cerr << __FILE__ << " " << __LINE__ << __FUNCTION__ << " haven't implemented!" << std::endl;
    }

    struct Np{
      Eigen::Vector3d vertex_;
      size_t o_id_;
    };

    static bool ispre(Np a, Np b) {      
      for(size_t i=0; i<3; ++i) {
        if(a.vertex_[i] < b.vertex_[i]) return true;
        else if(a.vertex_[i] > b.vertex_[i]) return false;
      }
      return true;
    }

    void sortMeshVertex(Eigen::Matrix<size_t, 3, Eigen::Dynamic> &mesh, Eigen::Matrix<zsw::Scalar, 3, Eigen::Dynamic> &node)
    {
      std::vector<Np> nps(node.cols());
      for(size_t i=0; i<node.cols(); ++i) {
        nps[i].vertex_ = node.block<3,1>(0,i);
        nps[i].o_id_ = i;
      }
      std::sort(nps.begin(), nps.end(), ispre);
      std::vector<size_t> o2n(nps.size());
      for(size_t i=0; i<nps.size(); ++i) {
        o2n[nps[i].o_id_] = i;
      }

      // modify
      for(size_t i=0; i<node.cols(); ++i) {
        node.block<3,1>(0,i) = nps[i].vertex_;
      }
      for(size_t i=0; i<mesh.cols(); ++i) {
        mesh(0,i) = o2n[mesh(0,i)];
        mesh(1,i) = o2n[mesh(1,i)];
        mesh(2,i) = o2n[mesh(2,i)];
      }
    }
  }
}

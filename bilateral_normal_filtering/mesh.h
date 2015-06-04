#ifndef MESH_H
#define MESH_H

#include <Eigen/Dense>

namespace zsw {
  namespace mesh {

    typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrixd;
    typedef Eigen::Matrix<size_t, 3, Eigen::Dynamic> Matrixst;

    class edge2cell_adjacent
    {
    public:
      static edge2cell_adjacent *create(const Matrixst & mesh,
                                        const bool &show_debug_info = true);

      static edge2cell_adjacent *create_on_dirty_mesh(const Matrixst & mesh);

      size_t get_edge_idx(size_t vi, size_t vj) const;

      size_t get_edge_idx(const size_t *v) const;

      std::pair<size_t,size_t> query(size_t vi, size_t vj) const;

      std::pair<size_t,size_t> query(const size_t *v) const;

      static inline bool is_boundary_edge(const std::pair<size_t,size_t> &nb_tri_id) {
        return (nb_tri_id.first == -1) ^ (nb_tri_id.second == -1);
      }

      Matrixst get_edge(size_t id) const {
        Matrixst rtn(2);
        rtn[0] = edges_[id].first;
        rtn[1] = edges_[id].second;
        return rtn;
      }
      std::vector<std::pair<size_t,size_t> > edges_;
      std::vector<std::pair<size_t,size_t> > edge2cell_;
    private:
      edge2cell_adjacent(){}

      ///
      /// @brief init initialize the edge2cell relationship
      /// @param mesh input manifold and orientable surface mesh
      /// @param show_debug_info a switcher to show information for debugging
      /// @return return 0 if ok, or non-zeros
      ///
      int init(const Matrixst &mesh, const bool & show_debug_info);

      ///
      /// @brief init initialize the edge2cell relationship
      /// @param mesh input manifold surface mesh
      /// @return return 0 if ok, or non-zeros
      ///
      int init_on_dirty_mesh(const Matrixst &mesh);
    };

    class V2CellAdjacent
    {
    public:
      static V2CellAdjacent * createV2CellAdj(const Eigen::Matrix<size_t, 3, Eigen::Dynamic> &tri);
      const std::vector<size_t> &query(const size_t vi) const;
    };

    class TriMesh{
      Eigen::Matrix<size_t, 3, Eigen::Dynamic> tri_;
      Matrixd vertex_;
      // extra data
      std::shared_ptr<Matrixd> face_normal_;
      std::shared_ptr<Eigen::VectorXd> face_area_;
      std::shared_ptr<edge2cell_adjacent> edge2cell_adj_;
      std::shared_ptr<V2CellAdjacent> v2cell_adj_;
    };
  }
}

#endif /* MESH_H */

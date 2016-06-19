#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

struct VertexInfo
{
  double Fs_;
  double fs_;
  VertexInfo() { Fs_=1e9; fs_=1e9; }
};

struct CellInfo
{
  double min_val_;
  double max_val_;
  CellInfo() {
    min_val_=1e9;
    max_val_=-1e9;
  }
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, K>    Vb;
typedef CGAL::Triangulation_cell_base_with_info_3<CellInfo, K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                    Tds;
//Use the Fast_location tag. Default or Compact_location works too.
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Delaunay;
typedef Delaunay::Point                                             Point;

int main(int argc, char *argv[])
{

  return 0;
}

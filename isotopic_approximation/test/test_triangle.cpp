#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <vector>
typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>    Vb;
typedef CGAL::Triangulation_data_structure_3<Vb>                    Tds;
//Use the Fast_location tag. Default or Compact_location works too.
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Delaunay;
typedef Delaunay::Point                                             Point;

int main()
{
  Delaunay T;
  Delaunay::Vertex_handle vh = T.insert( Point(0,0,0));  vh->info()=0;
  vh = T.insert(Point(1,0,0)); vh->info()=1;
  vh = T.insert(Point(0,1,0)); vh->info()=2;
  vh = T.insert(Point(0,0,1)); vh->info()=3;
  vh = T.insert(Point(2,2,2)); vh->info()=4;
  vh = T.insert(Point(-1,0,1)); vh->info()=5;
  assert( T.is_valid() );
  for(Delaunay::Finite_cells_iterator cit=T.finite_cells_begin();
      cit!=T.finite_cells_end(); ++cit) {
    for(int i=0; i<4; ++i) {
      Delaunay::Vertex_handle vh = cit->vertex(i);
      std::cerr << vh->info() << " ";
    }
    std::cerr << std::endl;
  }
  return 0;
}

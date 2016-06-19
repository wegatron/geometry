#include "basic_data_structure.h"

using namespace zsw;

int main()
{
    std::vector<P_Info> p;
    p.reserve(8);
    p.push_back(std::make_pair(Point(0, 0, 0), VertexInfo()));
    p.push_back(std::make_pair(Point(0, 0, 1), VertexInfo()));
    p.push_back(std::make_pair(Point(0, 1, 0), VertexInfo()));
    p.push_back(std::make_pair(Point(1, 0, 0), VertexInfo()));
    p.push_back(std::make_pair(Point(1, 0, 1), VertexInfo()));
    p.push_back(std::make_pair(Point(1, 1, 0), VertexInfo()));
    p.push_back(std::make_pair(Point(0, 1, 1), VertexInfo()));
    p.push_back(std::make_pair(Point(1, 1, 1), VertexInfo()));
    TriangulationWapper t(p);
    TTds &td = t.delaunay_triangulation_.tds();
    for(auto i=td.edges_begin(); i!=td.edges_end(); ++i)
        t.isSatisfyLinkCondition(*i);

    return 0;
}

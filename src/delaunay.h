#ifndef DELAUNAY_H
#define DELAUNAY_H 

#include "graph.h"
#include "dcel.h"
#include "utility.h"

class delaunay {
    // incremental delaunay triangulation algorithm
public:
    delaunay(const std::vector<util::point>& points);

    const dcel& triangulation() const;

    // returns edges of the triangulation without 
    // imaginary points point_minus_1 and point_minus_2
    std::vector<util::line_segment> get_edges() const;

    // returns [min_x, max_x, min_y, max_y]
    // so we can determine the range of
    // coordinate system when drawing triangulation
    std::vector<double> range() const;

private:
    // init dcel structure with the biggest triangle
    void init_dcel();
    // init graph structure with the biggest triangle
    void init_graph();

    enum class position { strictly_interior, boundary, outside };
    // returns position of point in a triangle that is stored
    // in a node with the given node index
    position get_position(util::point point, int node_index) const;
    // returns node index that has a triangle which contains the point
    int find_node(util::point point) const;
    // add new point to the current triangulation
    // when all points are added, we will have delaunay triangulation
    void add_point(int point_index);

    // split the triangle when point is strictly inside the triangle
    // node_index - index of node in a graph that contains the triangle
    // that is to be splitted
    void split_triangle_interior(int point_index, int node_index);
    // split the triangle when point is on the triangle boundary
    // node_index - index of node in a graph that contains the triangle
    // that is to be splitted
    void split_triangle_boundary(int point_index, int node_index);

    void try_flip(dcel::edgeref<false> edge);
    void flip_edge(dcel::edgeref<false> edge);

    graph m_graph; 
    dcel m_dcel;
};

#endif /* DELAUNAY_H */

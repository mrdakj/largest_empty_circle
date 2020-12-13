#ifndef VORONOI_H
#define VORONOI_H 

#include "dcel.h"
#include "utility.h"

class voronoi {
public:
    voronoi(const dcel& triangulation);

    const dcel& graph() const;

    // returns voronoi edges
    std::vector<util::line_segment> get_edges() const;

    // returns [min_x, max_x, min_y, max_y]
    // so we can determine the range of
    // coordinate system when drawing voronoi diagram
    // skip points at infinity
    std::vector<double> range() const;

private:
    // get voronoi vertices from delaunay triangles
    void circumcenters();
    // get voronoi vertex for imaginary delaunay face
    util::point get_external_center(dcel::edgeref<true> edge, util::point centre) const;

    // add new point to the current voronoi graph
    // when all points are added, we will have voronoi graph
    void add_point(int point_index);

    // delaunay triangulation
    const dcel& m_triangulation;
    // voronoi graph
    dcel m_dcel;
};

#endif /* VORONOI_H */

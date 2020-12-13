#ifndef CONVEX_HULL_H
#define CONVEX_HULL_H 

#include "delaunay.h"
#include "utility.h"

class convex_hull {
public:
    convex_hull(const dcel& triangulation);

    const std::vector<util::line_segment> edges() const;

    // returns true if point p is inside a convex hull (including edges)
    bool inside(const util::point& p) const;

    // returns intersection points of interval (origin, destination) and convex hull edges
    // there can be up to two intersection points
    std::vector<util::point> get_inersection(const util::point& origin, const util::point& destination);

private:
    std::vector<util::point> get_vertices(const dcel& triangulation) const;
    std::vector<util::line_segment> get_edges(const std::vector<util::point>& convex_hull_vertices) const;

    // convex hull edges in positive direction
    std::vector<util::line_segment> m_edges;
};

#endif /* CONVEX_HULL_H */

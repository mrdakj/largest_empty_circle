#ifndef LARGEST_EMPTY_CIRCLE_H
#define LARGEST_EMPTY_CIRCLE_H 

#include "utility.h"
#include "dcel.h"
#include "convex_hull.h"

class largest_empty_circle {
public:
    largest_empty_circle(const dcel& delaunay, const dcel& voronoi);

    const std::vector<util::circle>& candidates() const;
    util::circle get_largest_circle() const;

private:
    // candidate empty circles, 
    // the largest empty circle is in candidates
    std::vector<util::circle> m_candidates;
};

#endif /* LARGEST_EMPTY_CIRCLE_H */

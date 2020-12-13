#include "convex_hull.h"
#include <algorithm>

convex_hull::convex_hull(const dcel& triangulation)
    : m_edges(get_edges(get_vertices(triangulation)))
{
}

std::vector<util::point> convex_hull::get_vertices(const dcel& triangulation) const
{
    std::vector<util::point> convex_hull_vertices;

    // get an edge departing from 0 point
    auto current_edge = triangulation.edge(triangulation.vertex(0).incident_edge()-1);

    // get edge 0 -> dceltype::point_minus_2
    while (!(current_edge.origin() == 1 && (current_edge | edgerelation::twin).origin() == dceltype::point_minus_2)) {
        // Get next edge departing from 0 point.
        current_edge = current_edge | edgerelation::previous | edgerelation::twin;
    }

    // get previous edge as it is the first edge in the convex hull
    current_edge = current_edge | edgerelation::previous;
    int first_edge_id = current_edge.id();

    // get convex hull vertices
    do {
        // insert next point
        convex_hull_vertices.push_back(current_edge.point());
        // get next edge
        current_edge = current_edge | edgerelation::previous | edgerelation::twin | edgerelation::previous;
        // if point is imaginary then skip the edge
        if (current_edge.origin() < 0) {
            current_edge = current_edge | edgerelation::twin | edgerelation::previous;
        }
    } while (current_edge.id() != first_edge_id);

    return convex_hull_vertices;
}

const std::vector<util::line_segment> convex_hull::edges() const
{
    return m_edges;
}

std::vector<util::line_segment> convex_hull::get_edges(const std::vector<util::point>& convex_hull_vertices) const
{
    std::vector<util::line_segment> result;
    for (int i = 0; i < (int)convex_hull_vertices.size()-1; ++i) {
        result.emplace_back(convex_hull_vertices[i], convex_hull_vertices[i+1]);
    }

    result.emplace_back(convex_hull_vertices.back(), convex_hull_vertices[0]);
    return result;
}

bool convex_hull::inside(const util::point& p) const
{
    // convex hull is given in positive direction
    return std::none_of(m_edges.cbegin(), m_edges.cend(), 
               [&](const auto& convex_hull_edge) {
                   // point is (not striclty) inside the convex hull if there is no negative turn
                   return convex_hull_edge.origin().get_direction(convex_hull_edge.destination(), p) == util::direction::negative;
               });
}

std::vector<util::point> convex_hull::get_inersection(const util::point& origin, const util::point& destination)
{
    std::vector<util::point> intersections;
    auto input_segment = util::line_segment(origin,destination);

    for (const auto& e : m_edges) {
        auto intersection = input_segment.intersection_point(util::line_segment(e.origin(), e.destination()));
        if (intersection) {
            if (*intersection == origin || *intersection == destination) {
                continue;
            }
            intersections.emplace_back(*intersection);
        }
    }

    return intersections;
}


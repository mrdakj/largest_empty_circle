#include "largest_empty_circle.h"
#include <unordered_map>
#include <algorithm>
#include <cassert>

largest_empty_circle::largest_empty_circle(const dcel& delaunay, const dcel& voronoi)
{
    convex_hull ch(delaunay);

    std::unordered_map<int, bool> inside_map;

    auto check_point = [&](int point_id, const util::point& point) {
        // returns true if point is (not strictly) inside convex hull
        // and updates candidates and map
        if (inside_map.find(point_id) != inside_map.end()) {
            return inside_map[point_id];
        }
        else {
            if ((inside_map[point_id] = ch.inside(point))) {
                // voronoi vertex id corresponds to delaunay face id
                auto point_in_face = delaunay.point(point_id);
                m_candidates.emplace_back(point, point.distance(point_in_face));
            }

            return inside_map[point_id];
        }
    };

    for (int i = 0; i < voronoi.edge_count(); i += 2) {
        auto voronoi_edge = voronoi.edge(i);
        auto origin = voronoi_edge.point();
        int origin_id = voronoi_edge.origin();
        auto destination = (voronoi_edge | edgerelation::twin).point();
        int destination_id = (voronoi_edge | edgerelation::twin).origin();

        if (!check_point(origin_id, origin) || !check_point(destination_id, destination)) {
            // check intersection only if at least one point is outside the convex hull
            // it is important to node that intersection can exist even if both points are
            // outside of the convex hull
            auto intersections = ch.get_inersection(origin, destination);
            if (!intersections.empty()) {
                auto vertices_ids_origin_face = delaunay.points_ids(origin_id);
                auto vertices_ids_destination_face = delaunay.points_ids(destination_id);
                assert(vertices_ids_origin_face.size() == 3);
                assert(vertices_ids_destination_face.size() == 3);

                // find a vertex id that belongs to both faces
                auto it = std::find_if(vertices_ids_origin_face.begin(), vertices_ids_origin_face.end(), [&](int id) { 
                              return std::any_of(vertices_ids_destination_face.begin(), vertices_ids_destination_face.end(),
                                      [&](int id2) { return id == id2; });
                          });

                assert(it != vertices_ids_origin_face.end());
                auto p = delaunay.vertex(*it-1).point();
                std::transform(intersections.cbegin(), intersections.cend(), std::back_inserter(m_candidates), 
                    [&](const auto& intersection_point) {
                        return util::circle(intersection_point, intersection_point.distance(p));
                    });
            }
        }
    }
}

const std::vector<util::circle>& largest_empty_circle::candidates() const
{
    return m_candidates;
}

util::circle largest_empty_circle::get_largest_circle() const
{
    // get a circle with maximum radius
    return *std::max_element(m_candidates.begin(), m_candidates.end(), 
            [](const auto& lhs, const auto& rhs) { return lhs.r() < rhs.r(); });
}

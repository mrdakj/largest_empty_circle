#include "voronoi.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>

#define INF (100000)

voronoi::voronoi(const dcel& triangulation)
    : m_triangulation(triangulation)
{
    m_dcel.add(dceltype::face{-1});

    // create voronoi vertices
    circumcenters();

    // Compute Voronoi area for every point in Delaunay triangulation.
    for (int i=0; i < m_triangulation.vertex_count(); ++i) {
        add_point(i);
    }
}

void voronoi::add_point(int point_index)
{
    int new_edge_id = m_dcel.edge_count()+1;
    int newVoronoiFaceId = m_dcel.face_count();

    // Get edge departing from point.
    auto current_edge = m_triangulation.edge(m_triangulation.vertex(point_index).incident_edge()-1);
    while (m_triangulation.imaginary(current_edge.face()) || current_edge.external()) {
        current_edge = current_edge | edgerelation::previous | edgerelation::twin;
    }
    auto first_processed_edge_id = current_edge.id();

    int previous_edge_id = -1;
    int first_porcessed_edge_id = -1;

    do {
        // Skip imaginary edges.
        if (!current_edge.has_negative_vertex()) {
            // Get origin and destination points of Voronoi edge.
            int voronoiOriginId = (current_edge | edgerelation::twin).face();
            int voronoiDestId = current_edge.face();

            // Get current triangulation edge data.
            int edge_destination = (current_edge | edgerelation::twin).origin();

            // Check if Voronoi edge already created.
            if (edge_destination < current_edge.origin()) {
                // Get one of the edges of existing Voronoi face.
                auto existing_edge = m_dcel.edge(m_dcel.face(edge_destination).edge() - 1);

                // Loop until existing edge is found.
                while (voronoiDestId != existing_edge.origin()) {
                    existing_edge = existing_edge | edgerelation::next;
                }

                // Get twin of existing Voronoi edge (now existing edge is in current Voronoi face!!!!).
                existing_edge = existing_edge | edgerelation::twin;

                // Update existing edge.
                existing_edge.set_previous(previous_edge_id);
                existing_edge.set_next(new_edge_id);
                existing_edge.set_face(newVoronoiFaceId);

                if (previous_edge_id != -1) {
                    m_dcel.edge(previous_edge_id-1).set_next(existing_edge.id());
                }

                // Update previous edge.
                previous_edge_id = existing_edge.id();

                // Save edge if this is first edge in face.
                if (first_porcessed_edge_id == -1) {
                    first_porcessed_edge_id = existing_edge.id();
                }
            }
            // Current and twin edge did not exist -> create new edge and its twin.
            else
            {
                // If first edge inserted -> save its id.
                if (first_porcessed_edge_id == -1) {
                    first_porcessed_edge_id = new_edge_id;
                }

                // Add edge and its twin.
                m_dcel.add(dceltype::edge{voronoiOriginId, new_edge_id + 1, previous_edge_id, new_edge_id + 2, newVoronoiFaceId});
                m_dcel.add(dceltype::edge{voronoiDestId, new_edge_id, -1, -1, -1});

                // Update points.
                m_dcel.vertex(voronoiOriginId - 1).set_incident_edge(new_edge_id);
                m_dcel.vertex(voronoiDestId - 1).set_incident_edge(new_edge_id+1);

                // Update prvious and next edges id.
                previous_edge_id = new_edge_id;
                new_edge_id = new_edge_id + 2;
            }
        }
        // Update current edge and its index.
        current_edge = current_edge | edgerelation::previous | edgerelation::twin;
    } while (current_edge.id() != first_processed_edge_id);

    // Update "previous" edge for first inserted edge.
    m_dcel.edge(first_porcessed_edge_id-1).set_previous(previous_edge_id);
    // Update "next" edge for last inserted edge.
    m_dcel.edge(previous_edge_id-1).set_next(first_porcessed_edge_id);
    m_dcel.add(dceltype::face(first_porcessed_edge_id));
}

void voronoi::circumcenters()
{
    // voronoi vertices are created from circumcentres of delaunay's trinagles
    // we will create vertices such that its ids correspond to delaunay faces ids
    // voronoi vertex with id=x (index x-1) is created from a delaunay traingle
    // that forms a face with id=x
    //
    // - if delaunay triangle is a real face (it is not an external face (0) and has only real points)
    // then voronoi vertex is a circumcenter of that trinagle
    // - if delaunay triangle is imaginary and it doesn't contain both p_minus_2 and p_minus_1
    // then voronoi vertex is created from incident real face such that it lies on a line 
    // that contains the circumcenter of incident real face and that is perpendicular to
    // common edge of the imaginary face and the incident face; point is at infinity on that line
    // - if delaunay triangle is the bottom one (constains both point_minus_2 and point_minus_1)
    // that it is not incident to any real face, and voronoi vertex for it will be invalid and it
    // will not be used; we create voronoi vertex for that traingle just to keep id relation between
    // voronoi vertices ids and delaunay faces ids
    int	last_imaginary_face = 0;

    // external face (0) has no circumcenter
    for (int face_id = 1; face_id < m_triangulation.face_count(); ++face_id) {
        if (!m_triangulation.imaginary(face_id)) {
            // face is not imaginary - it is not an external face (0)
            // and has only real points
            auto triangle_points = m_triangulation.points(face_id);
            assert(triangle_points.size() == 3);
            util::circle circle{triangle_points[0], triangle_points[1], triangle_points[2]};
            m_dcel.add(dceltype::vertex{circle.center()});
        }
        else
        {
            // face is imaginary - it is external face (0)
            // or has point_minus_1 and/or point_minus_2
            last_imaginary_face = face_id;
            // add an inifinity point now, it will be updated below
            m_dcel.add(dceltype::vertex{util::point(INF, INF)});
        }
    }

    // create imaginary faces circumcenters
    for (int face_id=1; face_id <= last_imaginary_face; ++face_id) {
        if (m_triangulation.imaginary(face_id) && !m_triangulation.bottom(face_id)) {
            // face is imaginary and doesn't contain both p_minus_2 and p_minus_1
            auto edge = m_triangulation.face_edge(face_id);

            // find incident real face
            int neighbour_face_id = 0;
            while (true) {
                // check if twin edge belongs to real delaunay face
                neighbour_face_id = (edge | edgerelation::twin).face();
                if (!m_triangulation.imaginary(neighbour_face_id)) {
                    break;
                }
                else {
                    edge = edge | edgerelation::next;
                }
            }

            // circumcenter of the incident real face
            auto center = m_dcel.vertex(neighbour_face_id-1).point();
            m_dcel.vertex(face_id-1) = {get_external_center(edge, center)};
        }
    }
}

util::point voronoi::get_external_center(dcel::edgeref<true> edge, util::point center) const
{
    auto origin = edge.point();
    auto destination = (edge | edgerelation::twin).point();
    util::point middle_point {(origin.x()+destination.x())/2.0, (origin.y()+destination.y())/2.0};
    auto direction = origin.get_direction(destination, center);

    if (direction == util::direction::collinear) {
        // center is on the edge
        assert(middle_point == center);
        util::point rotated_point = destination.rotate_90(middle_point);
        util::point vec = {rotated_point.x() - middle_point.x(), rotated_point.y() - middle_point.y()};
        return {center.x() + INF*vec.x(), center.y() + INF*vec.y()};
    }
    else {
        // center is not on the edge
        // if center is inside the triangle (negative direction) then p1 is center, p2 is middle_point
        // if center is outside the triangle (positive direction) then p1 is middle_point, p2 is center
        auto p1 = (direction == util::direction::negative) ? center : middle_point;
        auto p2 = (direction == util::direction::negative) ? middle_point : center;

        util::point vec{p2.x()-p1.x(), p2.y()-p1.y()};

        // p2 + INF(p2-p1)
        return {p2.x() + INF*vec.x(), p2.y() + INF*vec.y()};
    }
}

std::vector<util::line_segment> voronoi::get_edges() const
{
    std::vector<util::line_segment> result;

    // skip twin edges
    for (int i = 0; i < m_dcel.edge_count(); i += 2) {
        auto edge = m_dcel.edge(i);
        auto origin = edge.point();
        auto destination = (edge | edgerelation::twin).point();
        if (origin != destination) {
            result.emplace_back(origin, destination);
        }
    }

    return result;
}

std::vector<double> voronoi::range() const
{
    double xmin = INF;
    double ymin = INF;
    double xmax = -INF;
    double ymax = -INF;

    for (int i = 0; i < m_dcel.vertex_count(); ++i) {
        // vertex id = delaunay face id = i+1
        if (m_triangulation.imaginary(i+1)) {
            // skip points at infinity
            continue;
        }
        xmin = std::min(xmin, m_dcel.vertex(i).point().x());
        ymin = std::min(ymin, m_dcel.vertex(i).point().y());
        xmax = std::max(xmax, m_dcel.vertex(i).point().x());
        ymax = std::max(ymax, m_dcel.vertex(i).point().y());
    }

    return {xmin, xmax, ymin, ymax};
}

const dcel& voronoi::graph() const
{
    return m_dcel;
}

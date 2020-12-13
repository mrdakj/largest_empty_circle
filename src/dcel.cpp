#include "dcel.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

#define EXTERNAL_FACE (0)

int dceltype::point_minus_2 = -2;
int dceltype::point_minus_1 = -1;

// vertex
dceltype::vertex::vertex(util::point point, int incident_edge)
    : m_point(std::move(point))
      , m_incident_edge(incident_edge)
{
}

bool dceltype::vertex::operator<(const vertex& other) const
{
    return m_point < other.m_point;
}

bool dceltype::vertex::operator>(const vertex& other) const
{
    return m_point > other.m_point;
}

bool dceltype::vertex::operator==(const vertex& other) const
{
    return m_point == other.m_point;
}

bool dceltype::vertex::operator!=(const vertex& other) const
{
    return m_point != other.m_point;
}

bool dceltype::vertex::operator<=(const vertex& other) const
{
    return m_point <= other.m_point;
}

bool dceltype::vertex::operator>=(const vertex& other) const
{
    return m_point >= other.m_point;
}

double dceltype::vertex::x() const
{
    return m_point.x();
}

double dceltype::vertex::y() const
{
    return m_point.y();
}

void dceltype::vertex::set_incident_edge(int incident_edge)
{
    m_incident_edge = incident_edge;
}

int dceltype::vertex::incident_edge() const
{
    return m_incident_edge;
}

util::point dceltype::vertex::point() const
{
    return m_point;
}

std::ostream& operator<<(std::ostream& out, const dceltype::vertex& v)
{
    return out << "(" << v.point() << "," << v.incident_edge() << ")";
}

// edge
dceltype::edge::edge(int origin, int twin, int previous, int next, int face)
    : m_origin(origin)
    , m_twin(twin)
    , m_previous(previous)
    , m_next(next)
    , m_face(face)
{
}

int dceltype::edge::origin() const
{
    return m_origin;
}

int dceltype::edge::twin() const
{
    return m_twin;
}

int dceltype::edge::previous() const
{
    return m_previous;
}

int dceltype::edge::next() const
{
    return m_next;
}

int dceltype::edge::face() const
{
    return m_face;
}

void dceltype::edge::set_origin(int origin)
{
    m_origin = origin;
}

void dceltype::edge::set_twin(int twin)
{
    m_twin = twin;
}

void dceltype::edge::set_previous(int previous)
{
    m_previous = previous;
}

void dceltype::edge::set_next(int next)
{
    m_next = next;
}

void dceltype::edge::set_face(int face)
{
    m_face = face;
}

std::ostream& operator<<(std::ostream& out, const dceltype::edge& e)
{
    return out << "(" << e.origin() << "," << e.twin() << "," << e.previous() << "," << e.next() << "," << e.face() << ")";
}

// face
dceltype::face::face(int edge_id)
    : m_edge(edge_id)
{
}

void dceltype::face::set_edge(int edge_id)
{
    m_edge = edge_id;
}

int dceltype::face::edge() const
{
    return m_edge;
}

std::ostream& operator<<(std::ostream& out, const dceltype::face& f)
{
    return out << "(" << f.edge() << ")";
}

// dcel
dcel::dcel(const std::vector<util::point>& points)
{
    if (points.size() < 3) {
        throw std::runtime_error("Minimum number of points is 3");
    }
    std::transform(points.cbegin(), points.cend(), std::back_inserter(m_vertices), [](util::point p) { return dceltype::vertex(p); });
}

void dcel::set_highest_first()
{
    assert(m_edges.empty() && m_faces.empty() && !m_vertices.empty());
    int index_of_highest = get_highest_vertex_index();
    std::swap(m_vertices[0], m_vertices[index_of_highest]);
}

int dcel::vertex_count() const
{
    return m_vertices.size();
}

int dcel::edge_count() const
{
    return m_edges.size();
}

int dcel::face_count() const
{
    return m_faces.size();
}

const dceltype::vertex& dcel::vertex(int i) const
{
    return m_vertices[i];
}

const dceltype::face& dcel::face(int i) const
{
    return m_faces[i];
}

dceltype::vertex& dcel::vertex(int i)
{
    return m_vertices[i];
}

dceltype::face& dcel::face(int i)
{
    return m_faces[i];
}

const std::vector<dceltype::vertex>& dcel::vertices() const
{
    return m_vertices;
}

void dcel::add(dceltype::vertex v)
{
    m_vertices.emplace_back(std::move(v));
}

void dcel::add(dceltype::edge e)
{
    m_edges.emplace_back(std::move(e));
}

void dcel::add(dceltype::face f)
{
    m_faces.emplace_back(std::move(f));
}

bool dcel::external_edge(int edge_index) const
{
    auto e = edge(edge_index);
    return e.face() == EXTERNAL_FACE || (e | edgerelation::twin).face() == EXTERNAL_FACE;
}

bool dcel::has_negative_vertex(int edge_index) const
{
    auto e = edge(edge_index);
    return e.origin() < 0 || (e | edgerelation::twin).origin() < 0;
}

int dcel::get_highest_vertex_index() const
{
    auto it = std::max_element(m_vertices.cbegin(), m_vertices.cend());
    return std::distance(m_vertices.cbegin(), it);
}

util::direction dcel::get_direction(util::point p, int source_point_id, int destination_point_id) const
{
    if (source_point_id > 0 && destination_point_id > 0) {
        // normal source and normal destination point
        // direction negative - point is not in the triangle
        // direction positive - point can be in the triangle (need to check other edges)
        // direction collinear - point is on the line segment of triangle edge
        return m_vertices[source_point_id-1].point().get_direction(m_vertices[destination_point_id-1].point(), p);
    }

    if (source_point_id > 0 && destination_point_id == dceltype::point_minus_2) {
        // if p is above source point, or at the same level and right of the source point,
        // then direction is negative, otherwise it is positive
        // note that we can choose p_minus_2 such that this holds for each point in input set, 
        // though we don't need to explicitly say where it is located
        //
        // p_minus_2 x              p_minus_2 x
        //            \     x p                \
        //             \   /                    \
        //              \ /                      \
        //        source x                 source x ---- x p
        return (p > m_vertices[source_point_id-1].point()) ? util::direction::negative : util::direction::positive;
    }

    if (source_point_id > 0 && destination_point_id == dceltype::point_minus_1) {
        // if p is above source point, or at the same level and right of the source point,
        // then direction is positive, otherwise it is negative
        // note that we can choose p_minus_1 such that this holds for each point in input set, 
        // though we don't need to explicitly say where it is located
        //
        // source x                 source x
        //         \     x p                \
        //          \   /                    \
        //           \ /                      \
        //  p_minus_1 x              p_minus_1 x ---- x p
        return (p > m_vertices[source_point_id-1].point()) ? util::direction::positive : util::direction::negative;
    }

    // source point is not a normal point - it is either p_minus_2 or p_minus_1

    if (source_point_id == dceltype::point_minus_1) {
        // destination cannot be p_minus_2, so destination is normal point
        return (p > m_vertices[destination_point_id-1].point()) ? util::direction::negative : util::direction::positive;
    }

    // source point is p_minus_2
    return (destination_point_id == dceltype::point_minus_1 || p > m_vertices[destination_point_id-1].point()) ? 
            util::direction::positive : util::direction::negative;
}

int dcel::collinear_edge_id(int point_index, int face_id) const
{
    auto face_edge = edge(face(face_id).edge()-1);
    int id1 = (face_edge | edgerelation::previous).origin();
    int id2 = face_edge.origin();
    int id3 = (face_edge | edgerelation::next).origin();

    auto p = m_vertices[point_index].point();

    return // point is collinear with the previous edge
           (get_direction(p, id1, id2) == util::direction::collinear) ?
            face_edge.previous() :
           // point is collinear with the input edge
           (get_direction(p, id2, id3) == util::direction::collinear) ?
            face_edge.id() :
           // point is collinear with the next edge
           (get_direction(p, id3, id1) == util::direction::collinear) ?
            face_edge.next() :
           // point is not collinear with any edge
            -1;
}

bool dcel::imaginary(int face_id) const
{
    auto vertices_ids = points_ids(face_id);
    return face_id == EXTERNAL_FACE || 
           std::any_of(vertices_ids.begin(), vertices_ids.end(), [](int id) { return id < 0; });
}

std::vector<int> dcel::points_ids(int face_id) const
{
    std::vector<int> result;
    auto edge_in_face = edge(face(face_id).edge()-1);
    int first_edge_id = edge_in_face.id();

    do {
        result.emplace_back(edge_in_face.origin());
        edge_in_face = edge_in_face | edgerelation::next;
    } while (edge_in_face.id() != first_edge_id);

    return result;
}

std::vector<util::point> dcel::points(int face_id) const
{
    std::vector<util::point> result;
    auto edge_in_face = edge(face(face_id).edge()-1);
    int first_edge_id = edge_in_face.id();

    do {
        if (edge_in_face.origin() > 0) {
            result.emplace_back(edge_in_face.point());
        }
        edge_in_face = edge_in_face | edgerelation::next;
    } while (edge_in_face.id() != first_edge_id);

    return result;
}

util::point dcel::point(int face_id) const
{
    auto edge_in_face = edge(face(face_id).edge()-1);
    int first_edge_id = edge_in_face.id();

    do {
        if (edge_in_face.origin() > 0) {
            return edge_in_face.point();
        }
        edge_in_face = edge_in_face | edgerelation::next;
    } while (edge_in_face.id() != first_edge_id);

    // we should not reach this
    assert(false);
    return util::point();
}

bool dcel::bottom(int face_id) const
{
    auto ids = points_ids(face_id);
    int imaginary_points_count = std::count_if(ids.begin(), ids.end(), [](int id) { return id < 0; });
    return imaginary_points_count == 2;
}

// edgeref
template<>
void dcel::edgeref<false>::set_origin(int origin)
{
    d->m_edges[edge_index].set_origin(origin);
}

template<>
void dcel::edgeref<false>::set_twin(int twin)
{
    d->m_edges[edge_index].set_twin(twin);
}

template<>
void dcel::edgeref<false>::set_previous(int previous)
{
    d->m_edges[edge_index].set_previous(previous);
}

template<>
void dcel::edgeref<false>::set_next(int next)
{
    d->m_edges[edge_index].set_next(next);
}

template<>
void dcel::edgeref<false>::set_face(int face)
{
    d->m_edges[edge_index].set_face(face);
}


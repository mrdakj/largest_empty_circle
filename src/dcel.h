#ifndef DCEL_H
#define DCEL_H 

#include <iostream>
#include <vector>
#include <type_traits>
#include <cassert>

#include "utility.h"

namespace dceltype {
    // Highest intput point (point with id 1), point_minus_2 and point_minus_1 
    // will form the first triangle. The triangle should contain all input points.
    // point_minus_2 has the property that for each input point A and input point B holds that:
    // if B > A, then point_minus_2 - A - B has positive orientation, negative otherwise.
    // point_minus_1 has the property that for each input point A and input point B holds that:
    // if B > A, then point_minus_2 - A - B has negative orientation, positive otherwise.
    // Note that we don't explicitly say where point_minus_1 and point_minus_2 are.
    extern int point_minus_2;
    extern int point_minus_1;

    class vertex {
    public:
        vertex(util::point data = {0,0}, int incident_edge = -1);

        bool operator<(const vertex& other) const;
        bool operator>(const vertex& other) const;
        bool operator==(const vertex& other) const;
        bool operator!=(const vertex& other) const;
        bool operator<=(const vertex& other) const;
        bool operator>=(const vertex& other) const;

        double x() const;
        double y() const;

        void set_incident_edge(int incident_edge);
        int incident_edge() const;
        util::point point() const;

    private:
        util::point m_point;
        // id of incident edge
        int m_incident_edge;
    };

    class edge {
    public:
        edge(int origin, int twin, int previous, int next, int face);

        int origin() const;
        int twin() const;
        int previous() const;
        int next() const;
        int face() const;

        void set_origin(int origin);
        void set_twin(int twin);
        void set_previous(int previous);
        void set_next(int next);
        void set_face(int face);

    private:
        int m_origin;
        int m_twin;
        int m_previous;
        int m_next;
        int m_face;
    };

    class face {
    public:
        face(int edge_id = -1);

        void set_edge(int edge_id);
        int edge() const;

    private:
        int m_edge;
    };

    std::ostream& operator<<(std::ostream& out, const dceltype::vertex& v);
    std::ostream& operator<<(std::ostream& out, const dceltype::edge& e);
    std::ostream& operator<<(std::ostream& out, const dceltype::face& f);
}

enum class edgerelation { twin, previous, next };

class dcel {
public:
    dcel() = default;
    dcel(const std::vector<util::point>& points);

    // sets the highest point to be the first in the vector
    // if there are more points with the same y coordiante
    // take the one with the maximum x coordinate
    void set_highest_first();

    // gets index of the highest point
    // if there are more points with the same y coordiante
    // take the one with the maximum x coordinate
    int get_highest_vertex_index() const;
    // returns direction of p - source point - destination point
    // direction can be positive, negative or collinear
    util::direction get_direction(util::point p, int source_point_id, int destination_point_id) const;

    int vertex_count() const;
    int edge_count() const;
    int face_count() const;

    const dceltype::vertex& vertex(int i) const;
    const dceltype::face& face(int i) const;

    dceltype::vertex& vertex(int i);
    dceltype::face& face(int i);

    const std::vector<dceltype::vertex>& vertices() const;

    void add(dceltype::vertex v);
    void add(dceltype::edge e);
    void add(dceltype::face f);

    // returns true if face with the given id is external (0)
    // or has point_minus_1 and/or point_minus_2
    bool imaginary(int face_id) const;
    // returns true if face with the given id constains both
    // point_minus_1 and point_minus_2
    bool bottom(int face_id) const;

    // returns edge id of edge that point lies on
    // edge is in face with given id
    // returns -1 if no such edge exists
    int collinear_edge_id(int point_index, int face_id) const;

    // returns all points in a face with the given id
    std::vector<util::point> points(int face_id) const;
    // returns vertices ids in a face with the given id
    std::vector<int> points_ids(int face_id) const;
    // returns some vertex point in a face with the given id
    util::point point(int face_id) const;

    // edgeref is used for syntax sugar
    template <bool t>
    struct edgeref {
    public:
        edgeref(int edge_index1, std::conditional_t<t==true, const dcel*, dcel*> d1)
            : edge_index(edge_index1)
            , d(d1)
        {}

        int origin() const
        {
            return d->m_edges[edge_index].origin();
        }

        int twin() const
        {
            return d->m_edges[edge_index].twin();
        }

        int previous() const
        {
            return d->m_edges[edge_index].previous();
        }

        int next() const
        {
            return d->m_edges[edge_index].next();
        }

        int face() const
        {
            return d->m_edges[edge_index].face();
        }

        int id() const
        {
            return edge_index+1;
        }

        void set_origin(int origin);
        void set_twin(int twin);
        void set_previous(int previous);
        void set_next(int next);
        void set_face(int face);

        edgeref operator|(edgerelation relation)
        {
            return get_related_edge(relation);
        }

        const edgeref operator|(edgerelation relation) const
        {
            return get_related_edge(relation);
        }

        bool external() const
        {
            return d->external_edge(edge_index);
        }

        bool has_negative_vertex() const
        {
            return d->has_negative_vertex(edge_index);
        }

        util::point point() const
        {
            return d->vertex(origin()-1).point();
        }

    private:
        edgeref get_related_edge(edgerelation relation) const
        {
            switch (relation) {
                case edgerelation::twin:
                    return {d->m_edges[edge_index].twin()-1, d};
                    break;
                case edgerelation::previous:
                    return {d->m_edges[edge_index].previous()-1, d};
                    break;
                case edgerelation::next:
                    return {d->m_edges[edge_index].next()-1, d};
                    break;
                default:
                    break;
            }

            // we should not reach this
            assert(false);
            return {-1, d};
        }

        int edge_index;
        // use pointer instead of reference so we can have operator=
        std::conditional_t<t==true, const dcel*, dcel*> d;
    };

    edgeref<false> edge(int i)
    {
        return {i, this};
    }

    edgeref<true> edge(int i) const
    {
        return {i, this};
    }

    edgeref<false> face_edge(int face_id)
    {
        return edge(face(face_id).edge()-1);
    }

    edgeref<true> face_edge(int face_id) const
    {
        return edge(face(face_id).edge()-1);
    }

private:
    // returns true if edge is incident to the external face
    bool external_edge(int edge_index) const;
    // returns true if edge constains point_minus_1 or point_minus_2
    bool has_negative_vertex(int edge_index) const;

    std::vector<dceltype::vertex> m_vertices;
    std::vector<dceltype::edge> m_edges;
    std::vector<dceltype::face> m_faces;
};

#endif /* DCEL_H */

#include "utility.h"
#include <cmath>
#include <iostream>
#include <cassert>

#define EPS (0.0001)
#define COLLINEAR_THRESHOLD (0.000001)

// point
util::point::point(double x, double y)
    : m_x(x)
    , m_y(y)
{
}

double util::point::distance(point p) const
{
    return std::sqrt(std::pow(p.x() - m_x, 2) + std::pow(p.y() - m_y, 2));
}

bool util::point::operator<(point other) const
{
    return (m_y < other.m_y) || (m_y == other.m_y && m_x < other.m_x);
}

bool util::point::operator>(point other) const
{
    return other < *this;
}

bool util::point::operator==(point other) const
{
    return std::fabs(m_x-other.m_x) < EPS && std::fabs(m_y-other.m_y) < EPS;
}

bool util::point::operator!=(point other) const
{
    return !(*this == other);
}

bool util::point::operator<=(point other) const
{
    return *this < other || *this == other;
}

bool util::point::operator>=(point other) const
{
    return *this > other || *this == other;
}

double util::point::x() const
{
    return m_x;
}

double util::point::y() const
{
    return m_y;
}

bool util::point::in_circle(point a, point b, point c) const
{
    // Let d be a determinant
    // | ax-x   ay-y  (ax-x)²+(ay-y)² |
    // | bx-x   by-y  (bx-x)²+(by-y)² |
    // | cx-x   cy-y  (cx-x)²+(cy-y)² |
    // if d = 0, then (x,y) is on the circle,
    // if d > 0, then (x,y) is in the circle,
    // if d < 0, then (x,y) is outside the circle

    double d11 = a.x()-m_x;
    double d12 = a.y()-m_y;
    double d13 = std::pow(a.x()-m_x,2) + std::pow(a.y()- m_y,2);

    double d21 = b.x()-m_x;
    double d22 = b.y()-m_y;
    double d23 = std::pow(b.x()-m_x,2) + std::pow(b.y()- m_y,2);

    double d31 = c.x()-m_x;
    double d32 = c.y()-m_y;
    double d33 = pow(c.x()-m_x,2) + pow(c.y()-m_y,2);

    return d11*d22*d33 + d12*d23*d31 + d13*d21*d32 - d13*d22*d31 - d23*d32*d11 - d33*d12*d21 > 0.0;
}

double util::point::signed_area(point p, point q) const
{
    // calculate signed area of triangle (this, p, q) using cross product
    return  (-p.m_x*m_y + q.m_x*m_y + m_x*p.m_y - q.m_x*p.m_y - m_x*q.m_y + p.m_x*q.m_y)/2;
}

util::direction util::point::get_direction(point p, point q) const
{
    double area = signed_area(p, q);
    return (area > COLLINEAR_THRESHOLD) ? 
            direction::positive :
           (area < -COLLINEAR_THRESHOLD) ? 
            direction::negative : 
            direction::collinear;
}

util::point util::point::rotate_90(point a) const
{
    double x = m_x;
    double y = m_y;

    // translate point for vector AO (O is origin)
    x -= a.x();
    y -= a.y();

    // rotate (x,y) for 90 degrees around origin
    double x_rotated = -y;
    double y_rotated = x;

    // translate point for vector OA (O is origin)
    x_rotated += a.x();
    y_rotated += a.y();

    return {x_rotated, y_rotated};
}

std::ostream& operator<<(std::ostream& out, const util::point& p)
{
    return out << "(" << p.x() << "," << p.y() << ")";
}

// circle
util::circle::circle(util::point center, double r)
    : m_center(std::move(center))
    , m_r(r)
{
}

util::circle::circle(util::point a, util::point b, util::point c)
    : m_center(get_center(a,b,c))
    , m_r(m_center.distance(a))
{
}

double util::circle::r() const
{
    return m_r;
}

util::point util::circle::center() const
{
    return m_center;
}

bool util::circle::operator==(const circle& other) const
{
    return m_center == other.m_center && std::fabs(m_r - other.m_r) < EPS;
}

util::point util::circle::get_center(util::point a, util::point b, util::point c) const
{
    // point a, b and c shouldn't be collinear
    assert(a.get_direction(b,c) != util::direction::collinear);

    util::point ab_middle{(a.x() + b.x())/2, (a.y() + b.y())/2};
    util::point bc_middle{(b.x() + c.x())/2, (b.y() + c.y())/2};
    util::point rotated_1 = b.rotate_90(ab_middle);
    util::point rotated_2 = c.rotate_90(bc_middle);
    auto center = util::line(ab_middle, rotated_1).intersection_point(util::line(bc_middle, rotated_2));
    //
    // there should be an intersection point
    assert(center);
    return *center;
}

// line_segment
util::line_segment::line_segment(util::point origin, util::point destination)
    : m_origin(std::move(origin))
    , m_destination(std::move(destination))
{
}

util::point util::line_segment::origin() const
{
    return m_origin;
}

util::point util::line_segment::destination() const
{
    return m_destination;
}

bool util::line_segment::operator==(const line_segment& other) const
{
    return (m_origin == other.m_origin && m_destination == other.m_destination) || 
           (m_origin == other.m_destination && m_destination == other.m_origin);
}

std::optional<util::point> util::line_segment::intersection_point(util::line_segment other) const
{
    double o1x = m_origin.x();
    double o1y = m_origin.y();
    double d1x = m_destination.x();
    double d1y = m_destination.y();

    double o2x = other.m_origin.x();
    double o2y = other.m_origin.y();
    double d2x = other.m_destination.x();
    double d2y = other.m_destination.y();

    if (std::fabs((d2y-o2y)*(d1x-o1x) - (d2x-o2x)*(d1y-o1y)) < EPS) {
        // cross product of parallel vectors is 0 vector
        auto between = [](util::point a, util::point b, util::point c) {
            // returns true if c is between a and b 
            return // if corss product of ab and ac is 0 vector then a, b and c are collinear
                   std::fabs((b.x()-a.x())*(c.y()-a.y()) - (c.x()-a.x())*(b.y()-a.y())) < EPS && 
                   // if x-projection of c is in x-projection of ab
                   c.x() <= std::max(a.x(),b.x()) && c.x() >= std::min(a.x(), b.x()) && 
                   // if y-projection of c is in y-projection of ab
                   c.y() <= std::max(a.y(),b.y()) && c.y() >= std::min(a.y(), b.y());
        };

        if (between(other.m_origin, other.m_destination, m_origin)) {
            // return just one point as inersection, 
            // although there can be infinitely many points that belong to both line segments
            return std::optional<util::point>{m_origin};
        }

        if (between(other.m_origin, other.m_destination, m_destination)) {
            // return just one point as inersection, 
            // although there can be infinitely many points that belong to both line segments
            return std::optional<util::point>{m_destination};
        }

        if (between(m_origin, m_destination, other.m_origin)) {
            // return just one point as inersection, 
            // although there can be infinitely many points that belong to both line segments
            return std::optional<util::point>{other.m_origin};
        }

        if (between(m_origin, m_destination, other.m_destination)) {
            // return just one point as inersection, 
            // although there can be infinitely many points that belong to both line segments
            return std::optional<util::point>{other.m_destination};
        }

        // there is no intersection
        return std::nullopt;
    }

    // lines are not parallel

    // solve the system:
    // intersection_point = o1 + t1*(d1-o1)
    // intersection_point = o2 + t2*(d2-o2)
    // system should have the solution since lines are not parallel
    // if t1 in [0,1] and t2 in [0,1] then intersection point
    // belongs to both line segments
    double t1,t2;

    if (std::fabs(d1x - o1x) < EPS) {
        t2 = (o1x-o2x) / (d2x-o2x);
        t1 = (o2y-o1y + t2*(d2y-o2y)) / (d1y-o1y);
    }
    else {
        t2 = ((o1y-o2y)*(d1x-o1x) - (o1x-o2x)*(d1y-o1y)) / ((d2y-o2y)*(d1x-o1x) - (d2x-o2x)*(d1y-o1y));
        t1 = (o2x-o1x + t2*(d2x-o2x)) / (d1x-o1x);
    }

    if (t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1) {
        assert(util::point(o1x + t1*(d1x-o1x), o1y + t1*(d1y-o1y)) == util::point(o2x + t2*(d2x-o2x), o2y + t2*(d2y-o2y)));
        // o1 + t1*(d1-o1)
        return std::optional<util::point>{util::point(o1x + t1*(d1x-o1x), o1y + t1*(d1y-o1y))};
    }

    // there is no intersection
    return std::nullopt;
}

// line
util::line::line(util::point origin, util::point destination)
    : m_origin(std::move(origin))
    , m_destination(std::move(destination))
{
}

util::point util::line::origin() const
{
    return m_origin;
}

util::point util::line::destination() const
{
    return m_destination;
}

std::optional<util::point> util::line::intersection_point(util::line other) const
{
    double o1x = m_origin.x();
    double o1y = m_origin.y();
    double d1x = m_destination.x();
    double d1y = m_destination.y();

    double o2x = other.m_origin.x();
    double o2y = other.m_origin.y();
    double d2x = other.m_destination.x();
    double d2y = other.m_destination.y();

    if (std::fabs((d2y-o2y)*(d1x-o1x) - (d2x-o2x)*(d1y-o1y)) < EPS) {
        // cross product of parallel vectors is 0 vector
        // there is no intersection
        return std::nullopt;
    }

    // lines are not parallel

    // solve the system:
    // intersection_point = o1 + t1*(d1-o1)
    // intersection_point = o2 + t2*(d2-o2)
    double t1,t2;

    if (std::fabs(d1x - o1x) < EPS) {
        t2 = (o1x-o2x) / (d2x-o2x);
        t1 = (o2y-o1y + t2*(d2y-o2y)) / (d1y-o1y);
    }
    else {
        t2 = ((o1y-o2y)*(d1x-o1x) - (o1x-o2x)*(d1y-o1y)) / ((d2y-o2y)*(d1x-o1x) - (d2x-o2x)*(d1y-o1y));
        t1 = (o2x-o1x + t2*(d2x-o2x)) / (d1x-o1x);
    }

    assert(util::point(o1x + t1*(d1x-o1x), o1y + t1*(d1y-o1y)) == util::point(o2x + t2*(d2x-o2x), o2y + t2*(d2y-o2y)));
    // o1 + t1*(d1-o1)
    return std::optional<util::point>{util::point(o1x + t1*(d1x-o1x), o1y + t1*(d1y-o1y))};
}


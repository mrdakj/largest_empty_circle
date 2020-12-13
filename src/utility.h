#ifndef UTILITY_H
#define UTILITY_H 

#include <vector>
#include <optional>
#include <iostream>

namespace util {
    enum class direction { positive, negative, collinear };

    class point {
    public:
        point(double x = 0, double y = 0);

        double x() const;
        double y() const;

        bool operator<(point other) const;
        bool operator>(point other) const;
        bool operator==(point other) const;
        bool operator!=(point other) const;
        bool operator<=(point other) const;
        bool operator>=(point other) const;

        double distance(point p) const;
        double signed_area(point p, point q) const;
        // returns direction of this - p - q
        direction get_direction(point p, point q) const;

        // returns true if this point is strictly
        // inside the circumcicle of (a,b,c)
        bool in_circle(point a, point b, point c) const;

        // rotate point around point a for 90 degree
        // in positive direction
        point rotate_90(point a) const;

    private:
        double m_x;
        double m_y;
    };

    class circle {
    public:
        circle(point center = {}, double r = 0);
        circle(util::point a, util::point b, util::point c);

        point center() const;
        double r() const;

        bool operator==(const circle& other) const;

    private:
        point get_center(util::point a, util::point b, util::point c) const;

        point m_center;
        double m_r;
    };

    class line_segment {
    public:
        line_segment(util::point origin, util::point destination);

        std::optional<util::point> intersection_point(line_segment other) const;

        bool operator==(const line_segment& other) const;

        util::point origin() const;
        util::point destination() const;
    private:
        util::point m_origin;
        util::point m_destination;
    };

    class line {
    public:
        line(util::point origin, util::point destination);

        std::optional<util::point> intersection_point(line other) const;

        util::point origin() const;
        util::point destination() const;
    private:
        util::point m_origin;
        util::point m_destination;
    };
}

std::ostream& operator<<(std::ostream& out, const util::point& p);

#endif /* UTILITY_H */

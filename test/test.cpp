#include <iostream>
#include <cmath>
#include <gtest/gtest.h>
#include "../src/dcel.h"
#include "../src/graph.h"
#include "../src/delaunay.h"
#include "../src/voronoi.h"
#include "../src/convex_hull.h"
#include "../src/largest_empty_circle.h"

#define EPS (0.0001)

// test: point 
TEST(point, comparison) 
{
    auto check_equal = [&](util::point a, util::point b) {
        // a == b
        ASSERT_EQ(a,b); // a == b
        ASSERT_LE(a,b); // a <= b
        ASSERT_LE(b,a); // b <= a
        ASSERT_GE(a,b); // a >= b
        ASSERT_GE(b,a); // b >= a
        ASSERT_FALSE(a < b);
        ASSERT_FALSE(b < a);
        ASSERT_FALSE(a != b);
        ASSERT_FALSE(b != a);
    };

    check_equal({1,2}, {1,2});
    check_equal({1.3,2}, {1.3,2});
    check_equal({1.3,2.5}, {1.3,2.5});

    auto check_less_than = [&](util::point lhs, util::point rhs) {
        // lhs < rhs
        ASSERT_LT(lhs,rhs); // lhs < rhs
        ASSERT_GT(rhs,lhs); // rhs > lhs
        ASSERT_LE(lhs,rhs); // lhs <= rhs
        ASSERT_GE(rhs,lhs); // rhs >= lhs
        ASSERT_NE(lhs,rhs); // lhs != rhs
        ASSERT_NE(rhs,lhs); // rhs != lhs
    };

    check_less_than({1,1}, {1,2});
    check_less_than({1,1}, {3,2});
    check_less_than({3,1}, {1,2});
    check_less_than({1,2}, {3,2});
}

TEST(point, signed_area) 
{
    ASSERT_EQ(util::point(0,0).signed_area({2,0}, {0,2}), 2);
    ASSERT_EQ(util::point(0,2).signed_area({0,0}, {2,0}), 2);
    ASSERT_EQ(util::point(2,0).signed_area({0,2}, {0,0}), 2);

    ASSERT_EQ(util::point(0,0).signed_area({0,2}, {2,0}), -2);
    ASSERT_EQ(util::point(2,0).signed_area({0,0}, {0,2}), -2);
    ASSERT_EQ(util::point(0,2).signed_area({2,0}, {0,0}), -2);

    // signed area of 3 same points is 0
    ASSERT_EQ(util::point(0,0).signed_area({0,0}, {0,0}), 0);
    // signed area of 2 same points is 0
    ASSERT_EQ(util::point(0,0).signed_area({1,0}, {1,0}), 0);
    // signed area of collinear points is 0
    ASSERT_EQ(util::point(0,0).signed_area({2,0}, {3,0}), 0);

    // checked signed area with double values
    ASSERT_TRUE(std::fabs(util::point(4.5,8).signed_area({2.8,1}, {10,3.4}) - 23.16) < EPS);
    ASSERT_TRUE(std::fabs(util::point(4.5,8).signed_area({10,3.4}, {2.8,1}) - (-23.16)) < EPS);
}

TEST(point, direction) 
{
    ASSERT_EQ(util::point(0,0).get_direction({2,0}, {0,2}), util::direction::positive);
    ASSERT_EQ(util::point(0,0).get_direction({0,2}, {2,0}), util::direction::negative);
    ASSERT_EQ(util::point(2,0).get_direction({0,2}, {0,0}), util::direction::positive);
    ASSERT_EQ(util::point(2,0).get_direction({0,0}, {0,2}), util::direction::negative);

    // direction of collinear points is collinear
    ASSERT_EQ(util::point(0,0).get_direction({2,0}, {3,0}), util::direction::collinear);
    ASSERT_EQ(util::point(0,0).get_direction({0,2}, {0,3}), util::direction::collinear);
    ASSERT_EQ(util::point(0,0).get_direction({2,2}, {3,3}), util::direction::collinear);

    // direction of 3 same points is collinear
    ASSERT_EQ(util::point(0,0).get_direction({0,0}, {0,0}), util::direction::collinear);
    // direction of 2 same points is collinear
    ASSERT_EQ(util::point(0,0).get_direction({1,1.5}, {1,1.5}), util::direction::collinear);

    // direction with double points
    ASSERT_EQ(util::point(4.5,8).get_direction({2.8,1}, {10,3.4}), util::direction::positive);
    ASSERT_EQ(util::point(4.5,8).get_direction({10,3.4}, {2.8,1}), util::direction::negative);
}

// test:line_segment
TEST(line_segment, intersection_point) 
{
    // not parallel lines
    // intersection point at line segment endings
    // o
    // |
    // |
    // x---------o
    ASSERT_EQ(*util::line_segment({0,0}, {1,0}).intersection_point(util::line_segment({0,0}, {0,1})), util::point(0,0));
    ASSERT_EQ(*util::line_segment({0,0}, {1,0}).intersection_point(util::line_segment({0,1}, {0,0})),  util::point(0,0));
    //           o
    //           |
    //           |
    // o---------x
    ASSERT_EQ(*util::line_segment({0,0}, {1,0}).intersection_point(util::line_segment({1,0}, {1,2})),  util::point(1,0));
    ASSERT_EQ(*util::line_segment({0,0}, {1,0}).intersection_point(util::line_segment({1,2}, {1,0})),  util::point(1,0));
    //          
    //             1 o
    //   1,2 x     
    //      /
    //     /
    //    /
    // 2 o
    ASSERT_EQ(*util::line_segment({3,3}, {6,4}).intersection_point(util::line_segment({0,0}, {3,3})),  util::point(3,3));

    // intersection point at one line segment ending
    //      o
    //      |
    //      |
    // o----x----o
    ASSERT_EQ(*util::line_segment({0,0}, {1,0}).intersection_point(util::line_segment({0.5,0}, {1,1})),  util::point(0.5,0));
    ASSERT_EQ(*util::line_segment({0,0}, {1,0}).intersection_point(util::line_segment({1,1}, {0.5,0})),  util::point(0.5,0));

    // intersection point inside both line segments
    //   o     
    //    \     o
    //     \   /
    //      \ /
    //       x
    //      / \
    //     /   o
    //    o
    //
    ASSERT_EQ(*util::line_segment({1,1}, {3,8}).intersection_point(util::line_segment({2,0.5}, {-2,5})),  util::point(42.0/37,109.0/74));

    // parallel lines
    // intersection at line segment endings
    // 1       1,2          2
    // o--------x-----------o
    ASSERT_EQ(*util::line_segment({0,0}, {1,0}).intersection_point(util::line_segment({1,0}, {5,0})),  util::point(1,0));
    ASSERT_EQ(*util::line_segment({0,0}, {1,0}).intersection_point(util::line_segment({5,0}, {1,0})),  util::point(1,0));
    // 2       2,1          1
    // o--------x-----------o
    ASSERT_EQ(*util::line_segment({1,0}, {5,0}).intersection_point(util::line_segment({0,0}, {1,0})),  util::point(1,0));
    // 1      2    1      2
    // o------o----x------o
    ASSERT_EQ(*util::line_segment({0,0}, {1,0}).intersection_point(util::line_segment({0.5,0}, {2,0})),  util::point(1,0));
    // 2    1      2   1
    // o----x------o---o
    ASSERT_EQ(*util::line_segment({0.5,0}, {2,0}).intersection_point(util::line_segment({0,0}, {1,0})),  util::point(0.5,0));
    // 2    1      1   2
    // o----x------o---o
    ASSERT_EQ(*util::line_segment({0.5,0}, {1,0}).intersection_point(util::line_segment({0,0}, {2,0})),  util::point(0.5,0));
    // 1    2      2   1
    // o----x------o---o
    ASSERT_EQ(*util::line_segment({0,0}, {2,0}).intersection_point(util::line_segment({0.5,0}, {1,0})),  util::point(0.5,0));
    //       1 o
    //        /
    //   1,2 x
    //      /
    //     /
    //    /
    // 2 o
    ASSERT_EQ(*util::line_segment({3,3}, {4,4}).intersection_point(util::line_segment({0,0}, {3,3})),  util::point(3,3));
    //       1 o
    //        /
    //     2 o
    //      /
    //   1 x
    //    /
    // 2 o
    ASSERT_EQ(*util::line_segment({2,2}, {4,4}).intersection_point(util::line_segment({0,0}, {3,3})),  util::point(2,2));
}

TEST(line_segment, no_intersection_point) 
{
    // not parallel lines
    // o
    // |
    // |
    // o o---------o
    ASSERT_FALSE(util::line_segment({0.2,0}, {1,0}).intersection_point(util::line_segment({0,0}, {0,1})));
    ASSERT_FALSE(util::line_segment({0.2,0}, {1,0}).intersection_point(util::line_segment({0,1}, {0,0})));
    //             o
    //             |
    //             |
    // o---------o o
    ASSERT_FALSE(util::line_segment({0,0}, {1,0}).intersection_point(util::line_segment({1.2,0}, {2,0})));
    ASSERT_FALSE(util::line_segment({0,0}, {1,0}).intersection_point(util::line_segment({2,0}, {1.2,0})));
    //              1 o
    //         1 o     
    //     2 o     
    //      /
    //     /
    //    /
    // 2 o
    ASSERT_FALSE(util::line_segment({3.2,3}, {6,4}).intersection_point(util::line_segment({0,0}, {3,3})));

    //    o
    //    |
    //    |
    //    o
    // o---------o
    ASSERT_FALSE(util::line_segment({0,0}, {1,0}).intersection_point(util::line_segment({0.5,0.2}, {1,1})));
    ASSERT_FALSE(util::line_segment({0,0}, {1,0}).intersection_point(util::line_segment({1,1}, {0.5,0.2})));

    //          1 o
    //           /
    //          /
    //         /
    //        /  
    //       /    
    //    1 o  
    //       
    // 2 o
    //    \
    //     \
    //      \
    //     2 o
    ASSERT_FALSE(util::line_segment({1,1}, {3,8}).intersection_point(util::line_segment({2, -5}, {-2,-4})));

    // parallel lines
    // intersection at line segment endings
    // 1       1 2          2
    // o-------o o----------o
    ASSERT_FALSE(util::line_segment({0,0}, {1,0}).intersection_point(util::line_segment({1.2,0}, {5,0})));
    ASSERT_FALSE(util::line_segment({0,0}, {1,0}).intersection_point(util::line_segment({5,0}, {1.2,0})));
    // 2       2 1          1
    // o-------o o----------o
    ASSERT_FALSE(util::line_segment({1.2,0}, {5,0}).intersection_point(util::line_segment({0,0}, {1,0})));
    // 1       1
    // o-------o
    //            2          2
    //            o----------o
    ASSERT_FALSE(util::line_segment({0,0}, {1,0}).intersection_point(util::line_segment({1.5,-1}, {2,-1})));
    // 1       1
    // o-------o
    //       2          2
    //       o----------o
    ASSERT_FALSE(util::line_segment({0,0}, {1,0}).intersection_point(util::line_segment({0.5,-1}, {2,-1})));
    //          1 o
    //           /
    //          /
    //       1 o
    //       
    //     2 o
    //      /
    //     /
    //    /
    // 2 o
    ASSERT_FALSE(util::line_segment({4,4}, {5,5}).intersection_point(util::line_segment({0,0}, {3,3})));
    //     2 o    o 1
    //      /    /
    //     /    o 1
    //    /
    // 2 o
    ASSERT_FALSE(util::line_segment({2,1}, {4,3}).intersection_point(util::line_segment({0,0}, {3,3})));
}

// line
TEST(line, intersection_point) 
{
    // not parallel lines
    //   |
    //   o
    //   |
    //   |
    // --x---------o---
    //   |
    ASSERT_EQ(*util::line({0,0}, {1,0}).intersection_point(util::line({0,0}, {0,1})), util::point(0,0));
    ASSERT_EQ(*util::line({0,0}, {1,0}).intersection_point(util::line({0,1}, {0,0})),  util::point(0,0));
    //             |
    //             o
    //             |
    //             |
    // --o---------x--
    //             |
    ASSERT_EQ(*util::line({0,0}, {1,0}).intersection_point(util::line({1,0}, {1,2})),  util::point(1,0));
    ASSERT_EQ(*util::line({0,0}, {1,0}).intersection_point(util::line({1,2}, {1,0})),  util::point(1,0));
    //   o
    //   |
    //   |
    // --x-o---------o--
    //   |
    ASSERT_EQ(*util::line({0.2,0}, {1,0}).intersection_point(util::line({0,0}, {0,1})), util::point(0,0));
    //         /
    //        /    1 o
    //   1,2 x     
    //      /
    //     /
    //    /
    // 2 o
    //  /
    ASSERT_EQ(*util::line({3,3}, {6,4}).intersection_point(util::line({0,0}, {3,3})),  util::point(3,3));
    //  \
    //   o       /
    //    \     o
    //     \   /
    //      \ /
    //       x
    //      / \
    //     /   o
    //    o     \
    //   /
    ASSERT_EQ(*util::line({1,1}, {3,8}).intersection_point(util::line({2,0.5}, {-2,5})),  util::point(42.0/37,109.0/74));
    //          1 o
    //           /
    //          /
    //         /
    //        /  
    //       /    
    //    1 o  
    //       
    // 2 o
    //    \
    //     \
    //      \
    //     2 o
    ASSERT_EQ(*util::line({1,1}, {3,8}).intersection_point(util::line({2, -5}, {-2,-4})), util::point(-8.0/15,-131.0/30));
}

TEST(line, no_intersection_point) 
{
    // parallel lines
    // 1       1,2          2
    // o--------x-----------o
    ASSERT_FALSE(util::line({0,0}, {1,0}).intersection_point(util::line({1,0}, {5,0})));
    ASSERT_FALSE(util::line({0,0}, {1,0}).intersection_point(util::line({5,0}, {1,0})));
    // 2       2,1          1
    // o--------x-----------o
    ASSERT_FALSE(util::line({1,0}, {5,0}).intersection_point(util::line({0,0}, {1,0})));
    // 1      2    1      2
    // o------o----x------o
    ASSERT_FALSE(util::line({0,0}, {1,0}).intersection_point(util::line({0.5,0}, {2,0})));
    // 2    1      2   1
    // o----x------o---o
    ASSERT_FALSE(util::line({0.5,0}, {2,0}).intersection_point(util::line({0,0}, {1,0})));
    // 2    1      1   2
    // o----x------o---o
    ASSERT_FALSE(util::line({0.5,0}, {1,0}).intersection_point(util::line({0,0}, {2,0})));
    // 1    2      2   1
    // o----x------o---o
    ASSERT_FALSE(util::line({0,0}, {2,0}).intersection_point(util::line({0.5,0}, {1,0})));
    //       1 o
    //        /
    //   1,2 x
    //      /
    //     /
    //    /
    // 2 o
    ASSERT_FALSE(util::line({3,3}, {4,4}).intersection_point(util::line({0,0}, {3,3})));
    //       1 o
    //        /
    //     2 o
    //      /
    //   1 x
    //    /
    // 2 o
    ASSERT_FALSE(util::line({2,2}, {4,4}).intersection_point(util::line({0,0}, {3,3})));
    //          1 o
    //           /
    //          /
    //       1 o
    //       
    //     2 o
    //      /
    //     /
    //    /
    // 2 o
    ASSERT_FALSE(util::line({4,4}, {5,5}).intersection_point(util::line({0,0}, {3,3})));
    //     2 o    o 1
    //      /    /
    //     /    o 1
    //    /
    // 2 o
    ASSERT_FALSE(util::line({2,1}, {4,3}).intersection_point(util::line({0,0}, {3,3})));
}

// test: circle
TEST(circle, creation) 
{
    //    3
    //    |\
    //    | \
    //    |  x O
    //    |   \
    //    |    \
    //  1 ------ 2
    util::circle c({0,0}, {6,0}, {0,8});
    ASSERT_EQ(c, util::circle({3,4}, 5));

    c = util::circle ({1,10}, {2,5}, {3,6});
    ASSERT_EQ(c, util::circle({2.0/3,22.0/3}, 2.6874));
}

// test: vertex
TEST(vertex, comparison) 
{
    auto check_equal = [&](dceltype::vertex a, dceltype::vertex b) {
        // a == b
        ASSERT_EQ(a,b); // a == b
        ASSERT_LE(a,b); // a <= b
        ASSERT_LE(b,a); // b <= a
        ASSERT_GE(a,b); // a >= b
        ASSERT_GE(b,a); // b >= a
        ASSERT_FALSE(a < b);
        ASSERT_FALSE(b < a);
        ASSERT_FALSE(a != b);
        ASSERT_FALSE(b != a);
    };

    check_equal({{1,2}}, {{1,2}});
    check_equal({{1.3,2}}, {{1.3,2}});
    check_equal({{1.3,2.5}}, {{1.3,2.5}});

    auto check_less_than = [&](dceltype::vertex lhs, dceltype::vertex rhs) {
        // lhs < rhs
        ASSERT_LT(lhs,rhs); // lhs < rhs
        ASSERT_GT(rhs,lhs); // rhs > lhs
        ASSERT_LE(lhs,rhs); // lhs <= rhs
        ASSERT_GE(rhs,lhs); // rhs >= lhs
        ASSERT_NE(lhs,rhs); // lhs != rhs
        ASSERT_NE(rhs,lhs); // rhs != lhs
    };

    check_less_than({{1,1}}, {{1,2}});
    check_less_than({{1,1}}, {{3,2}});
    check_less_than({{3,1}}, {{1,2}});
    check_less_than({{1,2}}, {{3,2}});
}

// test: dcel
TEST(dcel, external_edge) 
{
    dcel d({{1,2}, {4,5}, {5,6}});
    d.add(dceltype::edge{1,4,7,8,1});
    d.add(dceltype::edge{dceltype::point_minus_2,6,9,10,2});
    d.add(dceltype::edge{dceltype::point_minus_1,5,11,12,3});
    d.add(dceltype::edge{dceltype::point_minus_2,1,6,5,0});
    d.add(dceltype::edge{1,3,4,6,0});
    d.add(dceltype::edge{dceltype::point_minus_1,2,5,4,0});

    d.add(dceltype::edge{2,12,8,1,1});
    d.add(dceltype::edge{dceltype::point_minus_2,9,1,7,1});
    d.add(dceltype::edge{2,8,10,2,2});
    d.add(dceltype::edge{dceltype::point_minus_1,11,2,9,2});
    d.add(dceltype::edge{2,10,12,3,3});
    d.add(dceltype::edge{1,7,3,11,3});

    for (int i = 0; i < 6; ++i) {
        ASSERT_TRUE(d.edge(i).external());
    }

    for (int i = 6; i < d.edge_count(); ++i) {
        ASSERT_FALSE(d.edge(i).external());
    }
}

TEST(dcel, has_negative_vertex) 
{
    dcel d({{1,2}, {4,5}, {5,6}});
    d.add(dceltype::edge{1,4,7,8,1});
    d.add(dceltype::edge{dceltype::point_minus_2,6,9,10,2});
    d.add(dceltype::edge{dceltype::point_minus_1,5,11,12,3});
    d.add(dceltype::edge{dceltype::point_minus_2,1,6,5,0});
    d.add(dceltype::edge{1,3,4,6,0});
    d.add(dceltype::edge{dceltype::point_minus_1,2,5,4,0});

    d.add(dceltype::edge{2,12,8,1,1});
    d.add(dceltype::edge{dceltype::point_minus_2,9,1,7,1});
    d.add(dceltype::edge{2,8,10,2,2});
    d.add(dceltype::edge{dceltype::point_minus_1,11,2,9,2});
    d.add(dceltype::edge{2,10,12,3,3});
    d.add(dceltype::edge{1,7,3,11,3});

    ASSERT_FALSE(d.edge(6).has_negative_vertex());
    ASSERT_FALSE(d.edge(11).has_negative_vertex());

    for (int i = 0; i < d.edge_count(); ++i) {
        if (i == 6 || i == 11) {
            continue;
        }
        ASSERT_TRUE(d.edge(i).has_negative_vertex());
    }
}

TEST(dcel, get_highest_vertex_index) 
{
    ASSERT_EQ(dcel({{1,2}, {4,5}, {5,6}}).get_highest_vertex_index(), 2);
    ASSERT_EQ(dcel({{1,2}, {5,6}, {4,5}}).get_highest_vertex_index(), 1);
    ASSERT_EQ(dcel({{5,6}, {1,2}, {4,5}}).get_highest_vertex_index(), 0);

    ASSERT_EQ(dcel({{5,6}, {6,6}, {4,5}}).get_highest_vertex_index(), 1);
    ASSERT_EQ(dcel({{5,6}, {6,6}, {10,6}}).get_highest_vertex_index(), 2);

    ASSERT_EQ(dcel({{-5,6}, {1,2}, {4,5}}).get_highest_vertex_index(), 0);

    dcel d({{5,6}, {60,5}, {70,3}});
    ASSERT_EQ(d.get_highest_vertex_index(), 0);
    d.add(dceltype::vertex{{7,4}});
    ASSERT_EQ(d.get_highest_vertex_index(), 0);
    d.add(dceltype::vertex{{7,9}});
    ASSERT_EQ(d.get_highest_vertex_index(), d.vertices().size()-1);
}

TEST(dcel, get_direction) 
{
    // id:    1      2      3      4      5      6      7      8        9         10      11
    dcel d({{0,0}, {2,0}, {0,2}, {2,2}, {3,0}, {0,3}, {3,3}, {1,1.5}, {4.5,8}, {2.8,1}, {10,3.4}});

    // normal source and destination point
    ASSERT_EQ(d.get_direction({0,0}, 2, 3), util::direction::positive);
    ASSERT_EQ(d.get_direction({0,0}, 3, 2), util::direction::negative);
    ASSERT_EQ(d.get_direction({2,0}, 3, 1), util::direction::positive);
    ASSERT_EQ(d.get_direction({2,0}, 1, 3), util::direction::negative);

    // direction of collinear points is collinear
    ASSERT_EQ(d.get_direction({0,0}, 2, 5), util::direction::collinear);
    ASSERT_EQ(d.get_direction({0,0}, 3, 6), util::direction::collinear);
    ASSERT_EQ(d.get_direction({0,0}, 4, 7), util::direction::collinear);

    // direction of 3 same points is collinear
    ASSERT_EQ(d.get_direction({0,0}, 1, 1), util::direction::collinear);
    // direction of 2 same points is collinear
    ASSERT_EQ(d.get_direction({0,0}, 8, 8), util::direction::collinear);

    // direction with double points
    ASSERT_EQ(d.get_direction({4.5,8}, 10, 11), util::direction::positive);
    ASSERT_EQ(d.get_direction({4.5,8}, 11, 10), util::direction::negative);

    // id:      1       2      3
    d = dcel({{0,10}, {4,8}, {5,7}});

    // normal source, destination point dceltype::point_minus_2
    // point (5,9) is above point (4,8) so direction is negative
    ASSERT_EQ(d.get_direction({5,9}, 2, dceltype::point_minus_2), util::direction::negative);
    // point (5,7) is below point (4,8) so direction is positive
    ASSERT_EQ(d.get_direction({5,7}, 2, dceltype::point_minus_2), util::direction::positive);
    // point (5,8) is greater than point (4,8) so direction is negative
    ASSERT_EQ(d.get_direction({5,8}, 2, dceltype::point_minus_2), util::direction::negative);
    // point (3,8) is less than point (4,8) so direction is positive
    ASSERT_EQ(d.get_direction({3,8}, 2, dceltype::point_minus_2), util::direction::positive);

    // normal source, destination point dceltype::point_minus_1
    // point (5,9) is above point (4,8) so direction is positive
    ASSERT_EQ(d.get_direction({5,9}, 2, dceltype::point_minus_1), util::direction::positive);
    // point (5,7) is below point (4,8) so direction is negative
    ASSERT_EQ(d.get_direction({5,7}, 2, dceltype::point_minus_1), util::direction::negative);
    // point (5,8) is greater than point (4,8) so direction is positive
    ASSERT_EQ(d.get_direction({5,8}, 2, dceltype::point_minus_1), util::direction::positive);
    // point (3,8) is less than point (4,8) so direction is negative
    ASSERT_EQ(d.get_direction({3,8}, 2, dceltype::point_minus_1), util::direction::negative);

    // source point dceltype::point_minus_2, destination point dceltype::point_minus_1 -> direction is always positive
    ASSERT_EQ(d.get_direction({3,8}, dceltype::point_minus_2, dceltype::point_minus_1), util::direction::positive);

    // source point dceltype::point_minus_2, normal destination point
    // point (5,9) is above point (4,8) so direction is positive
    ASSERT_EQ(d.get_direction({5,9}, dceltype::point_minus_2, 2), util::direction::positive);
    // point (5,7) is below point (4,8) so direction is negative
    ASSERT_EQ(d.get_direction({5,7}, dceltype::point_minus_2, 2), util::direction::negative);
    // point (5,8) is greater than point (4,8) so direction is positive
    ASSERT_EQ(d.get_direction({5,8}, dceltype::point_minus_2, 2), util::direction::positive);
    // point (3,8) is less than point (4,8) so direction is negative
    ASSERT_EQ(d.get_direction({3,8}, dceltype::point_minus_2, 2), util::direction::negative);

    // source point dceltype::point_minus_1, normal destination point
    // point (5,9) is above point (4,8) so direction is negative
    ASSERT_EQ(d.get_direction({5,9}, dceltype::point_minus_1, 2), util::direction::negative);
    // point (5,7) is below point (4,8) so direction is positive
    ASSERT_EQ(d.get_direction({5,7}, dceltype::point_minus_1, 2), util::direction::positive);
    // point (5,8) is greater than point (4,8) so direction is negative
    ASSERT_EQ(d.get_direction({5,8}, dceltype::point_minus_1, 2), util::direction::negative);
    // point (3,8) is less than point (4,8) so direction is positive
    ASSERT_EQ(d.get_direction({3,8}, dceltype::point_minus_1, 2), util::direction::positive);
}

// test: graph
TEST(graph, creation_and_face_map_check) 
{
    graph g;
    g.add({{1, dceltype::point_minus_2, dceltype::point_minus_1}, 1});
    g.add({{dceltype::point_minus_2, 2, 1}, 1});
    g.add({{dceltype::point_minus_1, 2, dceltype::point_minus_2}, 2});
    g.add({{1, 2, dceltype::point_minus_1}, 3});

    ASSERT_EQ(g.size(), 4);
    for (int i = 0; i < g.size(); ++i) {
        ASSERT_EQ(g[i].children_count(), 0);
    }

    g[0].set_children({1,2,3});

    ASSERT_EQ(g[0].children_count(), 3);
    for (int i = 1; i < g.size(); ++i) {
        ASSERT_EQ(g[i].children_count(), 0);
    }

    ASSERT_EQ(g.get_node(1), 1);
    ASSERT_EQ(g.get_node(2), 2);
    ASSERT_EQ(g.get_node(3), 3);
}

TEST(graph, creation_and_face_map_check_multiconnections) 
{
    graph g;
    g.add({{1, dceltype::point_minus_2, dceltype::point_minus_1}, 1});
    g.add({{dceltype::point_minus_2, 2, 1}, 1});
    g.add({{dceltype::point_minus_1, 2, dceltype::point_minus_2}, 2});
    g.add({{1, 2, dceltype::point_minus_1}, 3});
    g.add({{dceltype::point_minus_2, 3, 2}, 2});
    g.add({{dceltype::point_minus_1, 3, dceltype::point_minus_2}, 4});
    g.add({{2, 3, dceltype::point_minus_1}, 5});
    g.add({{2, 1, 3}, 2});
    g.add({{dceltype::point_minus_2, 3, 1}, 1});

    ASSERT_EQ(g.size(), 9);
    for (int i = 0; i < g.size(); ++i) {
        ASSERT_EQ(g[i].children_count(), 0);
    }

    g[0].set_children({1,2,3});
    g[1].set_children({7,8});
    g[2].set_children({4,5,6});
    g[4].set_children({7,8});

    ASSERT_EQ(g[0].children_count(), 3);
    ASSERT_EQ(g[1].children_count(), 2);
    ASSERT_EQ(g[2].children_count(), 3);
    ASSERT_EQ(g[3].children_count(), 0);
    ASSERT_EQ(g[4].children_count(), 2);
    ASSERT_EQ(g[5].children_count(), 0);
    ASSERT_EQ(g[6].children_count(), 0);
    ASSERT_EQ(g[7].children_count(), 0);
    ASSERT_EQ(g[8].children_count(), 0);

    ASSERT_EQ(g.get_node(1), 8);
    ASSERT_EQ(g.get_node(2), 7);
    ASSERT_EQ(g.get_node(3), 3);
    ASSERT_EQ(g.get_node(4), 5);
    ASSERT_EQ(g.get_node(5), 6);
}

// test: delaunay
TEST(delaunay, triangulation1) 
{
    delaunay del{{{0, 0}, {0, 1}, {1, 0}, {1, 1}, {0.5,0.5}}};
    auto edges = del.get_edges();
    ASSERT_TRUE(edges.size() == 12);
    ASSERT_TRUE(edges[0] == util::line_segment(util::point(0.5,0.5),util::point(1,1)));
    ASSERT_TRUE(edges[1] == util::line_segment(util::point(1,1),util::point(0,1)));
    ASSERT_TRUE(edges[2] == util::line_segment(util::point(0,1),util::point(0.5,0.5)));
    ASSERT_TRUE(edges[3] == util::line_segment(util::point(0,0),util::point(0.5,0.5)));
    ASSERT_TRUE(edges[4] == util::line_segment(util::point(0.5,0.5),util::point(0,1)));
    ASSERT_TRUE(edges[5] == util::line_segment(util::point(0,1),util::point(0,0)));
    ASSERT_TRUE(edges[6] == util::line_segment(util::point(0.5,0.5),util::point(0,0)));
    ASSERT_TRUE(edges[7] == util::line_segment(util::point(0,0),util::point(1,0)));
    ASSERT_TRUE(edges[8] == util::line_segment(util::point(1,0),util::point(0.5,0.5)));
    ASSERT_TRUE(edges[9] == util::line_segment(util::point(1,1),util::point(0.5,0.5)));
    ASSERT_TRUE(edges[10] == util::line_segment(util::point(0.5,0.5),util::point(1,0)));
    ASSERT_TRUE(edges[11] == util::line_segment(util::point(1,0),util::point(1,1)));
}

TEST(delaunay, triangulation2) 
{
    delaunay del = {{{0, 0}, {0, 1}, {1, 0}, {1, 1}}};
    auto edges = del.get_edges();
    ASSERT_TRUE(edges.size() == 6);
    ASSERT_TRUE(edges[0] == util::line_segment(util::point(1,0),util::point(1,1)));
    ASSERT_TRUE(edges[1] == util::line_segment(util::point(1,1),util::point(0,1)));
    ASSERT_TRUE(edges[2] == util::line_segment(util::point(0,1),util::point(1,0)));
    ASSERT_TRUE(edges[3] == util::line_segment(util::point(0,1),util::point(0,0)));
    ASSERT_TRUE(edges[4] == util::line_segment(util::point(0,0),util::point(1,0)));
    ASSERT_TRUE(edges[5] == util::line_segment(util::point(1,0),util::point(0,1)));
}

TEST(delaunay, triangulation3) 
{
    delaunay del = {{{0, 0}, {0, 1.1}, {1, 0}, {1, 1}, {2,3.4}, {0.2,1.6}}};
    auto edges = del.get_edges();
    ASSERT_TRUE(edges.size() == 15);
    ASSERT_TRUE(edges[0] == util::line_segment(util::point(1,0),util::point(2,3.4)));
    ASSERT_TRUE(edges[1] == util::line_segment(util::point(2,3.4),util::point(1,1)));
    ASSERT_TRUE(edges[2] == util::line_segment(util::point(1,1),util::point(1,0)));
    ASSERT_TRUE(edges[3] == util::line_segment(util::point(1,1),util::point(0,0)));
    ASSERT_TRUE(edges[4] == util::line_segment(util::point(0,0),util::point(1,0)));
    ASSERT_TRUE(edges[5] == util::line_segment(util::point(1,0),util::point(1,1)));
    ASSERT_TRUE(edges[6] == util::line_segment(util::point(0.2,1.6),util::point(1,1)));
    ASSERT_TRUE(edges[7] == util::line_segment(util::point(1,1),util::point(2,3.4)));
    ASSERT_TRUE(edges[8] == util::line_segment(util::point(2,3.4),util::point(0.2,1.6)));
    ASSERT_TRUE(edges[9] == util::line_segment(util::point(0,0),util::point(1,1)));
    ASSERT_TRUE(edges[10] == util::line_segment(util::point(1,1),util::point(0,1.1)));
    ASSERT_TRUE(edges[11] == util::line_segment(util::point(0,1.1),util::point(0,0)));
    ASSERT_TRUE(edges[12] == util::line_segment(util::point(1,1),util::point(0.2,1.6)));
    ASSERT_TRUE(edges[13] == util::line_segment(util::point(0.2,1.6),util::point(0,1.1)));
    ASSERT_TRUE(edges[14] == util::line_segment(util::point(0,1.1),util::point(1,1)));
}

TEST(delaunay, triangulation4) 
{
    delaunay del = {{{0, 0}, {0, 1.1}, {1, 0}, {1, 1}}};
    auto edges = del.get_edges();
    ASSERT_TRUE(edges.size() == 6);
    ASSERT_TRUE(edges[0] == util::line_segment(util::point(1,1),util::point(0,0)));
    ASSERT_TRUE(edges[1] == util::line_segment(util::point(0,0),util::point(1,0)));
    ASSERT_TRUE(edges[2] == util::line_segment(util::point(1,0),util::point(1,1)));
    ASSERT_TRUE(edges[3] == util::line_segment(util::point(0,0),util::point(1,1)));
    ASSERT_TRUE(edges[4] == util::line_segment(util::point(1,1),util::point(0,1.1)));
    ASSERT_TRUE(edges[5] == util::line_segment(util::point(0,1.1),util::point(0,0)));
}

TEST(delaunay, triangulation5) 
{
    delaunay del = {{{9, 1}, {2, 1.9}, {2, 0}, {0, 1.54}, {2,3.4}, {0.2,1.6}, {1.2,2.6}, {24, 12}}};
    auto edges = del.get_edges();
    ASSERT_TRUE(edges.size() == 27);
    ASSERT_TRUE(edges[0] == util::line_segment(util::point(2,1.9),util::point(1.2,2.6)));
    ASSERT_TRUE(edges[1] == util::line_segment(util::point(1.2,2.6),util::point(0.2,1.6)));
    ASSERT_TRUE(edges[2] == util::line_segment(util::point(0.2,1.6),util::point(2,1.9)));
    ASSERT_TRUE(edges[3] == util::line_segment(util::point(9,1),util::point(2,3.4)));
    ASSERT_TRUE(edges[4] == util::line_segment(util::point(2,3.4),util::point(2,1.9)));
    ASSERT_TRUE(edges[5] == util::line_segment(util::point(2,1.9),util::point(9,1)));
    ASSERT_TRUE(edges[6] == util::line_segment(util::point(9,1),util::point(2,1.9)));
    ASSERT_TRUE(edges[7] == util::line_segment(util::point(2,1.9),util::point(2,0)));
    ASSERT_TRUE(edges[8] == util::line_segment(util::point(2,0),util::point(9,1)));
    ASSERT_TRUE(edges[9] == util::line_segment(util::point(2,3.4),util::point(9,1)));
    ASSERT_TRUE(edges[10] == util::line_segment(util::point(9,1),util::point(24,12)));
    ASSERT_TRUE(edges[11] == util::line_segment(util::point(24,12),util::point(2,3.4)));
    ASSERT_TRUE(edges[12] == util::line_segment(util::point(0.2,1.6),util::point(2,0)));
    ASSERT_TRUE(edges[13] == util::line_segment(util::point(2,0),util::point(2,1.9)));
    ASSERT_TRUE(edges[14] == util::line_segment(util::point(2,1.9),util::point(0.2,1.6)));
    ASSERT_TRUE(edges[15] == util::line_segment(util::point(1.2,2.6),util::point(0,1.54)));
    ASSERT_TRUE(edges[16] == util::line_segment(util::point(0,1.54),util::point(0.2,1.6)));
    ASSERT_TRUE(edges[17] == util::line_segment(util::point(0.2,1.6),util::point(1.2,2.6)));
    ASSERT_TRUE(edges[18] == util::line_segment(util::point(2,0),util::point(0.2,1.6)));
    ASSERT_TRUE(edges[19] == util::line_segment(util::point(0.2,1.6),util::point(0,1.54)));
    ASSERT_TRUE(edges[20] == util::line_segment(util::point(0,1.54),util::point(2,0)));
    ASSERT_TRUE(edges[21] == util::line_segment(util::point(1.2,2.6),util::point(2,1.9)));
    ASSERT_TRUE(edges[22] == util::line_segment(util::point(2,1.9),util::point(2,3.4)));
    ASSERT_TRUE(edges[23] == util::line_segment(util::point(2,3.4),util::point(1.2,2.6)));
    ASSERT_TRUE(edges[24] == util::line_segment(util::point(0,1.54),util::point(1.2,2.6)));
    ASSERT_TRUE(edges[25] == util::line_segment(util::point(1.2,2.6),util::point(2,3.4)));
    ASSERT_TRUE(edges[26] == util::line_segment(util::point(2,3.4),util::point(0,1.54)));
}

TEST(delaunay, triangulation6) 
{
    delaunay del = {{{9, 1}, {2, 1.9}, {2, 0}, {5,5}, {1.2,4.9}, {2,3.4}, {0.2,1.6}, {1.2,2.6}, {24, 12}}};
    auto edges = del.get_edges();
    ASSERT_TRUE(edges.size() == 33);
    ASSERT_TRUE(edges[0] == util::line_segment(util::point(2,1.9),util::point(1.2,2.6)));
    ASSERT_TRUE(edges[1] == util::line_segment(util::point(1.2,2.6),util::point(0.2,1.6)));
    ASSERT_TRUE(edges[2] == util::line_segment(util::point(0.2,1.6),util::point(2,1.9)));
    ASSERT_TRUE(edges[3] == util::line_segment(util::point(2,1.9),util::point(9,1)));
    ASSERT_TRUE(edges[4] == util::line_segment(util::point(9,1),util::point(5,5)));
    ASSERT_TRUE(edges[5] == util::line_segment(util::point(5,5),util::point(2,1.9)));
    ASSERT_TRUE(edges[6] == util::line_segment(util::point(5,5),util::point(9,1)));
    ASSERT_TRUE(edges[7] == util::line_segment(util::point(9,1),util::point(24,12)));
    ASSERT_TRUE(edges[8] == util::line_segment(util::point(24,12),util::point(5,5)));
    ASSERT_TRUE(edges[9] == util::line_segment(util::point(24,12),util::point(1.2,4.9)));
    ASSERT_TRUE(edges[10] == util::line_segment(util::point(1.2,4.9),util::point(5,5)));
    ASSERT_TRUE(edges[11] == util::line_segment(util::point(5,5),util::point(24,12)));
    ASSERT_TRUE(edges[12] == util::line_segment(util::point(9,1),util::point(2,1.9)));
    ASSERT_TRUE(edges[13] == util::line_segment(util::point(2,1.9),util::point(2,0)));
    ASSERT_TRUE(edges[14] == util::line_segment(util::point(2,0),util::point(9,1)));
    ASSERT_TRUE(edges[15] == util::line_segment(util::point(2,1.9),util::point(0.2,1.6)));
    ASSERT_TRUE(edges[16] == util::line_segment(util::point(0.2,1.6),util::point(2,0)));
    ASSERT_TRUE(edges[17] == util::line_segment(util::point(2,0),util::point(2,1.9)));
    ASSERT_TRUE(edges[18] == util::line_segment(util::point(1.2,2.6),util::point(1.2,4.9)));
    ASSERT_TRUE(edges[19] == util::line_segment(util::point(1.2,4.9),util::point(0.2,1.6)));
    ASSERT_TRUE(edges[20] == util::line_segment(util::point(0.2,1.6),util::point(1.2,2.6)));
    ASSERT_TRUE(edges[21] == util::line_segment(util::point(2,3.4),util::point(2,1.9)));
    ASSERT_TRUE(edges[22] == util::line_segment(util::point(2,1.9),util::point(5,5)));
    ASSERT_TRUE(edges[23] == util::line_segment(util::point(5,5),util::point(2,3.4)));
    ASSERT_TRUE(edges[24] == util::line_segment(util::point(2,3.4),util::point(5,5)));
    ASSERT_TRUE(edges[25] == util::line_segment(util::point(5,5),util::point(1.2,4.9)));
    ASSERT_TRUE(edges[26] == util::line_segment(util::point(1.2,4.9),util::point(2,3.4)));
    ASSERT_TRUE(edges[27] == util::line_segment(util::point(1.2,2.6),util::point(2,1.9)));
    ASSERT_TRUE(edges[28] == util::line_segment(util::point(2,1.9),util::point(2,3.4)));
    ASSERT_TRUE(edges[29] == util::line_segment(util::point(2,3.4),util::point(1.2,2.6)));
    ASSERT_TRUE(edges[30] == util::line_segment(util::point(1.2,4.9),util::point(1.2,2.6)));
    ASSERT_TRUE(edges[31] == util::line_segment(util::point(1.2,2.6),util::point(2,3.4)));
    ASSERT_TRUE(edges[32] == util::line_segment(util::point(2,3.4),util::point(1.2,4.9)));
}

// convex hull
TEST(convex_hull, edges1) 
{
    delaunay del = {{{0, 0}, {0, 1.1}, {1, 0}, {1, 1}, {2,3.4}, {0.2,1.6}}};
    convex_hull ch{del.triangulation()};
    auto edges = ch.edges();
    ASSERT_EQ(edges.size(), 5);
    ASSERT_TRUE(edges[0] == util::line_segment(util::point(0.2,1.6),util::point(0,1.1)));
    ASSERT_TRUE(edges[1] == util::line_segment(util::point(0,1.1),util::point(0,0)));
    ASSERT_TRUE(edges[2] == util::line_segment(util::point(0,0),util::point(1,0)));
    ASSERT_TRUE(edges[3] == util::line_segment(util::point(1,0),util::point(2,3.4)));
    ASSERT_TRUE(edges[4] == util::line_segment(util::point(2,3.4),util::point(0.2,1.6)));
}

TEST(convex_hull, edges2) 
{
    delaunay del = {{{0, 0}, {0, 1}, {1, 0}, {1, 1}, {0.5,0.5}}};
    convex_hull ch{del.triangulation()};
    auto edges = ch.edges();
    ASSERT_EQ(edges.size(), 4);
    ASSERT_TRUE(edges[0] == util::line_segment(util::point(0,1),util::point(0,0)));
    ASSERT_TRUE(edges[1] == util::line_segment(util::point(0,0),util::point(1,0)));
    ASSERT_TRUE(edges[2] == util::line_segment(util::point(1,0),util::point(1,1)));
    ASSERT_TRUE(edges[3] == util::line_segment(util::point(1,1),util::point(0,1)));
}

TEST(convex_hull, edges3) 
{
    delaunay del = {{{9, 1}, {2, 1.9}, {2, 0}, {0, 1.54}, {2,3.4}, {0.2,1.6}, {1.2,2.6}, {24, 12}}};
    convex_hull ch{del.triangulation()};
    auto edges = ch.edges();
    ASSERT_EQ(edges.size(), 5);
    ASSERT_TRUE(edges[0] == util::line_segment(util::point(2,3.4),util::point(0,1.54)));
    ASSERT_TRUE(edges[1] == util::line_segment(util::point(0,1.54),util::point(2,0)));
    ASSERT_TRUE(edges[2] == util::line_segment(util::point(2,0),util::point(9,1)));
    ASSERT_TRUE(edges[3] == util::line_segment(util::point(9,1),util::point(24,12)));
    ASSERT_TRUE(edges[4] == util::line_segment(util::point(24,12),util::point(2,3.4)));
}

TEST(convex_hull, edges4) 
{
    delaunay del = {{{4,4},{3,12},{15,4},{12,17},{2,19},{6,10},{2,13},{4,12},{14,13},{16,7},{12,4},{3,8},{3,3},{13,19},{3,16},{15,2},{16,17},{13,14},{6,4},{3,11}}};
    convex_hull ch{del.triangulation()};
    auto edges = ch.edges();
    ASSERT_EQ(edges.size(), 7);
    ASSERT_TRUE(edges[0] == util::line_segment(util::point(2,19),util::point(2,13)));
    ASSERT_TRUE(edges[1] == util::line_segment(util::point(2,13),util::point(3,3)));
    ASSERT_TRUE(edges[2] == util::line_segment(util::point(3,3),util::point(15,2)));
    ASSERT_TRUE(edges[3] == util::line_segment(util::point(15,2),util::point(16,7)));
    ASSERT_TRUE(edges[4] == util::line_segment(util::point(16,7),util::point(16,17)));
    ASSERT_TRUE(edges[5] == util::line_segment(util::point(16,17),util::point(13,19)));
    ASSERT_TRUE(edges[6] == util::line_segment(util::point(13,19),util::point(2,19)));
}

TEST(convex_hull, point_inside) 
{
    delaunay del = {{{0, 0}, {0, 1.1}, {1, 0}, {1, 1}, {2,3.4}, {0.2,1.6}}};
    convex_hull ch{del.triangulation()};
    // all vertices of convex hull are inside convex hull
    ASSERT_TRUE(ch.inside({0,0}));
    ASSERT_TRUE(ch.inside({0,1.1}));
    ASSERT_TRUE(ch.inside({1,0}));
    ASSERT_TRUE(ch.inside({1,1}));
    ASSERT_TRUE(ch.inside({2,3.4}));
    ASSERT_TRUE(ch.inside({0.2,1.6}));

    // point on edge is inside convex hull
    ASSERT_TRUE(ch.inside({0,0.5}));

    // point strictly in is inside convex hull
    ASSERT_TRUE(ch.inside({0.5,0.5}));
}

TEST(convex_hull, interval_intersection) 
{
    delaunay del = {{{0, 0}, {0, 1.1}, {1, 0}, {1, 1}, {2,3.4}, {0.2,1.6}}};
    convex_hull ch{del.triangulation()};

    // one intersection point
    ASSERT_EQ(ch.get_inersection({0.5,0.5}, {0.5,-1}), std::vector<util::point>{util::point(0.5,0)});

    // two intersection points
    std::vector<util::point> intersection_points{{0,0.5}, {0.5,0}};
    ASSERT_EQ(ch.get_inersection({-1,1.5}, {1.5,-1}), intersection_points);

    // no intersection with one touching point - interval is considered, not segment
    ASSERT_TRUE(ch.get_inersection({0,0}, {-1,0}).empty());

    // no intersection, interval is inside
    ASSERT_TRUE(ch.get_inersection({0.1,0.1}, {0.2, 0.2}).empty());

    // no intersection, interval is outside
    ASSERT_TRUE(ch.get_inersection({-1,-1}, {-4, -1}).empty());
}

// voronoi
TEST(voronoi, edge_count1) 
{
    delaunay del{{{0, 0}, {0, 1}, {1, 0}, {1, 1}, {0.5,0.5}}};
    voronoi vor{del.triangulation()};
    auto voronoi_edges = vor.get_edges();
    ASSERT_EQ(voronoi_edges.size(), 8);
}

TEST(voronoi, edge_count2) 
{
    delaunay del{{{0, 0}, {0, 1}, {1, 0}, {1, 1}}};
    voronoi vor{del.triangulation()};
    auto voronoi_edges = vor.get_edges();
    ASSERT_EQ(voronoi_edges.size(), 4);
}

TEST(voronoi, edge_count3) 
{
    delaunay del{{{0, 0}, {0, 1.1}, {1, 0}, {1, 1}, {2,3.4}, {0.2,1.6}}};
    voronoi vor{del.triangulation()};
    auto voronoi_edges = vor.get_edges();
    ASSERT_EQ(voronoi_edges.size(), 10);
}

TEST(voronoi, edge_count4) 
{
    delaunay del{{{0, 0}, {0, 1.1}, {1, 0}, {1, 1}}};
    voronoi vor{del.triangulation()};
    auto voronoi_edges = vor.get_edges();
    ASSERT_EQ(voronoi_edges.size(), 5);
}

TEST(voronoi, edge_count5) 
{
    delaunay del{{{9, 1}, {2, 1.9}, {2, 0}, {0, 1.54}, {2,3.4}, {0.2,1.6}, {1.2,2.6}, {24, 12}}};
    voronoi vor{del.triangulation()};
    auto voronoi_edges = vor.get_edges();
    ASSERT_EQ(voronoi_edges.size(), 16);
}

TEST(voronoi, edge_count6) 
{
    delaunay del{{{9, 1}, {2, 1.9}, {2, 0}, {5,5}, {1.2,4.9}, {2,3.4}, {0.2,1.6}, {1.2,2.6}, {24, 12}}};
    voronoi vor{del.triangulation()};
    auto voronoi_edges = vor.get_edges();
    ASSERT_EQ(voronoi_edges.size(), 19);
}

TEST(voronoi, edge_count7) 
{
    delaunay del{{{15,6},{0,3},{17,6},{4,0},{18,5},{9,17},{4,7},{4,12},{10,4},{16,13}}};
    voronoi vor{del.triangulation()};
    auto voronoi_edges = vor.get_edges();
    ASSERT_EQ(voronoi_edges.size(), 21);
}

TEST(voronoi, edge_count8) 
{
    delaunay del{{{4,4},{3,12},{15,4},{12,17},{2,19},{6,10},{2,13},{4,12},{14,13},{16,7},{12,4},{3,8},{3,3},{13,19},{3,16},{15,2},{16,17},{13,14},{6,4},{3,11}}};
    voronoi vor{del.triangulation()};
    auto voronoi_edges = vor.get_edges();
    ASSERT_EQ(voronoi_edges.size(), 50);
}

TEST(largest_empty_circle, circle_and_candidates) 
{
    auto check_circle_and_candidates = [](std::vector<util::point> points, int expected_candidate_size, const util::circle& expected_largest_circle) {
        delaunay del{std::move(points)};
        voronoi vor{del.triangulation()};

        largest_empty_circle lec(del.triangulation(), vor.graph());
        auto candidates = lec.candidates();
        auto largest_circle = lec.get_largest_circle();
        ASSERT_EQ(candidates.size(), expected_candidate_size);
        ASSERT_EQ(largest_circle, expected_largest_circle);
    };

    check_circle_and_candidates(
        {{0, 0}, {0, 1}, {1, 0}, {1, 1}}, // points
        6,                                // expected candidate size
        {{0.5, 0.5}, 0.707107});          // expected largest circle
    check_circle_and_candidates(
        {{0, 0}, {0, 1.1}, {1, 0}, {1, 1}}, 
        6,
        {{0.45, 0.55}, 0.710634});
    check_circle_and_candidates(
        {{0, 0}, {0, 1.1}, {1, 0}, {1, 1}, {2,3.4}, {0.2,1.6}}, 
        10,
        {{1.32857, 2.27143}, 1.3132});
    check_circle_and_candidates(
        {{9, 1}, {2, 1.9}, {2, 0}, {0, 1.54}, {2,3.4}, {0.2,1.6}, {1.2,2.6}, {24, 12}}, 
        14,
        {{15.0363, 8.496}, 9.62427});
    check_circle_and_candidates(
        {{9, 1}, {2, 1.9}, {2, 0}, {5,5}, {1.2,4.9}, {2,3.4}, {0.2,1.6}, {1.2,2.6}, {24, 12}}, 
        16,
        {{14.6108, 9.07617}, 9.83391});
    check_circle_and_candidates(
        {{15,6},{0,3},{17,6},{4,0},{18,5},{9,17},{4,7},{4,12},{10,4},{16,13}}, 
        18,
        {{10.1765, 10.3824}, 6.38479});
    check_circle_and_candidates(
        {{4,4},{3,12},{15,4},{12,17},{2,19},{6,10},{2,13},{4,12},{14,13},{16,7},{12,4},{3,8},{3,3},{13,19},{3,16},{15,2},{16,17},{13,14},{6,4},{3,11}}, 
        38, 
        {{7.22222, 19}, 5.17949});
    check_circle_and_candidates(
        {{0, 0}, {0, 1}, {1, 0}, {1, 1}, {0.5,0.5}}, 
        4,
        {{0.5, 1}, 0.5});
}

TEST(largest_empty_circle, triangle) 
{
    auto check_circle = [](util::point a, util::point b, util::point c) {
        delaunay del{{a,b,c}};
        voronoi vor{del.triangulation()};

        largest_empty_circle lec(del.triangulation(), vor.graph());
        auto largest_circle = lec.get_largest_circle();

        ASSERT_EQ(largest_circle, util::circle(a,b,c));
    };

    // center of circumcircle should be in the triangle
    // so the largest and the circumcircle are the same
    check_circle({0,0}, {0,1}, {1,0});
    check_circle({0,0}, {5,5}, {10,0});
}

int main(int argc, char** argv) 
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}

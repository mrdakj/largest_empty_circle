find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable (test test.cpp)
target_link_libraries (test delaunay voronoi convex_hull largest_empty_circle ${GTEST_BOTH_LIBRARIES})

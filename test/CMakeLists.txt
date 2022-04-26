find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable (tests test.cpp)
target_link_libraries (tests delaunay voronoi convex_hull largest_empty_circle ${GTEST_BOTH_LIBRARIES})

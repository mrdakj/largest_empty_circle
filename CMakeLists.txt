cmake_minimum_required (VERSION 3.14)
# set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(THE_EMPTIEST_CIRCLE)

include(CTest)
add_subdirectory(src) 
add_subdirectory(test)
add_test(UnitTests ./test/test)

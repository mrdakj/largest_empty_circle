# The Largest Empty Circle 

The largest empty circle problem is the problem of finding the largest circle that contains no points in a given set and its center is inside a convex hull of the set.

# Algorithm

First Delaunay triangulation is found using incremental approach. Then, Voronoi diagram is constructed from the Delaunay triangulation. Center of the circle is either Voronoi vertex or the intersection point of a convex hull and a Voronoi edge.

![Alt text](out/4_circle_convex_hull.png?raw=true "Largest empty circle")

![Alt text](out/7_delaunay.png?raw=true "Delaunay triangulation")

![Alt text](out/7_voronoi.png?raw=true "Voronoi diagram")

![Alt text](out/7_circle.png?raw=true "Largest empty circle")

***
## :package: Installation
:exclamation: Requirements: C++17, OpenGL, cmake, GoogleTest

1. To install cmake using pacman package manager:

    ```sh
    pacman -S cmake

    ```

2. To install GoogleTest using pacman package manager:

    ```sh
    pacman -S gtest

    ```

### Manual

1. Clone this repository somewhere on your machine.

    ```sh
    git clone https://github.com/mrdakj/largest_empty_circle.git

    ```
2. Compile

    ```sh
    cd largest_empty_circle
    mkdir build
    cd build
    cmake ..
    (cmake -DCMAKE_BUILD_TYPE=Debug  ..) 
    cmake --build .

    ```

3. Run tests

    ```sh
    ctest

    ```

3. Run program

    ```sh
    ./src/main input_file [delaunay] [voronoi] [convex_hull] [circle] [all_circles]
    ./src/main ../input/1.txt delaunay circle 

    ```


<table>
  <tr>
    <th colspan="2">Options</th>
  </tr>
  <tr>
    <td>delaunay</td><td>constructs Delaunay triangulation</td>
  </tr>
  <tr>
    <td>voronoi</td><td>constructs Voronoi diagram</td>
  </tr>
  <tr>
    <td>convex_hull</td><td>constructs convex hull</td>
  </tr>
  <tr>
    <td>circle</td><td>finds the largest empty circle</td>
  </tr>
  <tr>
    <td>all_circles</td><td>displays all candidate circles</td>
  </tr>
</table>

#include <GL/glut.h>
#include <iostream>
#include <cmath>
#include <unordered_set>
#include <fstream>

#include "utility.h"
#include "dcel.h"
#include "delaunay.h"
#include "voronoi.h"
#include "convex_hull.h"
#include "largest_empty_circle.h"

enum class option { delaunay, voronoi, circle, all_circles, convex_hull, unknown };

option get_option(const std::string& option_string)
{
    return (option_string == "delaunay") ? option::delaunay : 
           (option_string == "voronoi") ? option::voronoi : 
           (option_string == "circle") ? option::circle : 
           (option_string == "all_circles") ? option::all_circles : 
           (option_string == "convex_hull") ? option::convex_hull : 
           option::unknown;
}

std::unordered_set<option> enabled_options;

std::vector<util::point> points;
std::vector<util::line_segment> delaunay_edges;
std::vector<util::line_segment> voronoi_edges;
std::vector<util::line_segment> convex_hull_edges;
std::vector<util::circle> candidates;
util::circle largest_circle;

void draw_circle(const util::circle& circle)
{
    glBegin(GL_LINE_LOOP);
    for(int i = 0; i < 300; ++i) {
        double angle =  2 * 3.14159 * i / 300;
        double x = std::cos(angle) * circle.r();
        double y = std::sin(angle) * circle.r();
        glVertex2d(circle.center().x()+x,circle.center().y()+y);
    }
    glEnd();
}

void draw_points()
{
    // green
    glColor3f(0.0f, 1.0f, 0.0f);
    glBegin(GL_POINTS);
    for (auto p : points) {
        glVertex2f(p.x(), p.y());
    }
    glEnd();
}

void draw_delaunay()
{
    // white
    glColor3f(1.0f, 1.0f, 1.0f);
    glBegin(GL_LINES);
    for (auto e : delaunay_edges) {
        glVertex2f(e.origin().x(), e.origin().y());
        glVertex2f(e.destination().x(), e.destination().y());
    }
    glEnd();

    // green
    glColor3f(0.0f, 1.0f, 0.0f);
    glBegin(GL_POINTS);
    for (auto e : delaunay_edges) {
        glVertex2f(e.origin().x(), e.origin().y());
    }
    glEnd();
}

void draw_voronoi()
{
    // white
    glColor3f(1.0f, 1.0f, 1.0f);
    glBegin(GL_LINES);
    for (auto e : voronoi_edges) {
        glVertex2f(e.origin().x(), e.origin().y());
        glVertex2f(e.destination().x(), e.destination().y());
    }
    glEnd();

    // red
    glColor3f(1.0f, 0.0f, 0.0f);
    glBegin(GL_POINTS);
    for (auto e : voronoi_edges) {
        glVertex2f(e.origin().x(), e.origin().y());
    }
    glEnd();
}

void draw_convex_hull()
{
    // white
    glColor3f(1.0f, 1.0f, 1.0f);
    glBegin(GL_LINES);
    for (auto e : convex_hull_edges) {
        glVertex2f(e.origin().x(), e.origin().y());
        glVertex2f(e.destination().x(), e.destination().y());
    }
    glEnd();
}

void draw_all_circles()
{
    // blue
    glColor3f(0, 0, 0.5);
    for (const auto& c : candidates) {
        draw_circle(c);
    }
}

void draw_largest_circle()
{
    // blue
    glColor3f(0, 0, 1);
    glBegin(GL_POINTS);
        glVertex2f(largest_circle.center().x(), largest_circle.center().y());
    glEnd();
    draw_circle(largest_circle);
}

void display() 
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glEnable(GL_POINT_SMOOTH);
    glPointSize(5.0);

    draw_points();

    // delaunay
    if (enabled_options.find(option::delaunay) != enabled_options.end()) {
        draw_delaunay();
    }

    // voronoi
    if (enabled_options.find(option::voronoi) != enabled_options.end()) {
        draw_voronoi();
    }

    // convex hull
    if (enabled_options.find(option::convex_hull) != enabled_options.end()) {
        draw_convex_hull();
    }
    
    // all candidate circles
    if (enabled_options.find(option::all_circles) != enabled_options.end()) {
        draw_all_circles();
    }

    // largest circle
    if (enabled_options.find(option::circle) != enabled_options.end()) {
        draw_largest_circle();
    }

    glFlush();
}

std::vector<util::point> read_points(std::ifstream& input_file)
{
    std::vector<util::point> points;

    if (input_file) {
        double x,y;
        char comma;
        while (input_file >> x >> comma >> y) {
            points.emplace_back(x,y);
        }
    }

    return points;
}

void init_window()
{
    glutInitWindowSize(600, 600);
    glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH)-600)/2, (glutGet(GLUT_SCREEN_HEIGHT)-600)/2);
    glutCreateWindow("The largest empty circle");
}

int main(int argc, char** argv) 
{
    if (argc < 2) {
        std::cout << "usage: ./main input_file [delaunay] [voronoi] [convex_hull] [circle] [all_circles]" << std::endl;
        return -1;
    }

    std::ifstream input_file(argv[1]);
    if (!input_file) {
        std::cout << "file not found" << std::endl;
        return -1;
    }

    points = read_points(input_file);

    for (int i = 2; i < argc; ++i) {
        enabled_options.emplace(get_option(argv[i]));
    }

    delaunay del{points};

    delaunay_edges = del.get_edges();
    auto range = del.range();

    voronoi vor{del.triangulation()};
    voronoi_edges = vor.get_edges();
    auto voronoi_range = vor.range();

    convex_hull ch{del.triangulation()};
    convex_hull_edges = ch.edges();

    largest_empty_circle lec(del.triangulation(), vor.graph());
    candidates = lec.candidates();
    largest_circle = lec.get_largest_circle();

    glutInit(&argc, argv);
    init_window();

    // use the same for x and y
    double min_number = std::min(range[0], range[2]) - 1;
    double max_maxnumber = std::max(range[1], range[3]) + 1;
    // set coordinate system
    gluOrtho2D(min_number, max_maxnumber, min_number, max_maxnumber);

    glutDisplayFunc(display);

    glutMainLoop();
    return 0;
}


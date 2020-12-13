#include "delaunay.h"
#include <algorithm>
#include <cassert>
#include <iostream>

delaunay::delaunay(const std::vector<util::point>& points)
    : m_dcel(points)
{
    init_dcel();
    init_graph();

    // first point is already added as it is a part of the biggest triangle
    for (int i = 1; i < m_dcel.vertex_count(); ++i) {
        add_point(i);
    }
}

void delaunay::init_dcel()
{
    // first point (id = 1, index = 0) is the highest one at it will be
    // usded to form the biggest triangle
    m_dcel.set_highest_first();
    m_dcel.vertex(0).set_incident_edge(1);

    // add 6 edges of the biggest triangle (1, point_minus_2, point_minus_1)
    // that contains all the points; there are 6 edges because
    // each edge has its twin edge in dcel structure
    // edge = {origin id, twin id, previous edge id, next edge id, face id}
    m_dcel.add({1, 4, 3, 2, 1});
    m_dcel.add({dceltype::point_minus_2, 6, 1, 3, 1});
    m_dcel.add({dceltype::point_minus_1, 5, 2, 1, 1});
    m_dcel.add({dceltype::point_minus_2, 1, 6, 5, 0});
    m_dcel.add({1, 3, 4, 6, 0});
    m_dcel.add({dceltype::point_minus_1, 2, 5, 4, 0});

    // add external face 
    // face = {incident edge id}
    m_dcel.add(dceltype::face{4});
    // add first internal face
    m_dcel.add(dceltype::face{1});
}

void delaunay::init_graph()
{
    // add first node in the graph - 
    // the biggest triangle (1, point_minus_2, point_minus_1) that constains
    // all the poinst, and it represents the first internal face (1)
    // node = {(point1 id, point2 id, point3 id), face id}
    m_graph.add(node({1, dceltype::point_minus_2, dceltype::point_minus_1}, 1));
}

delaunay::position delaunay::get_position(util::point point, int node_index) const
{
    // triangle points ids
    auto ids = m_graph[node_index].vertices();

    auto turn1 = m_dcel.get_direction(point, ids[0], ids[1]);
    auto turn2 = m_dcel.get_direction(point, ids[1], ids[2]);
    auto turn3 = m_dcel.get_direction(point, ids[2], ids[0]);

    return // all turns are positive
           (turn1 == util::direction::positive && turn2 == util::direction::positive && turn3 == util::direction::positive) ?
            position::strictly_interior :
           // there is no negative turn
           (turn1 != util::direction::negative && turn2 != util::direction::negative && turn3 != util::direction::negative) ?
            position::boundary :
           // there is a negative turn
            position::outside;
}

int delaunay::find_node(util::point point) const
{
    int current_index = 0;
    while (!m_graph[current_index].leaf()) {
        const auto& children = m_graph[current_index].children();
        auto it = std::find_if(children.cbegin(), children.cend(), [&](int child_index) { 
                      return get_position(point, child_index) != position::outside;
                  });
        // there should be a triangle that constains the point at each level
        assert(it != children.end());
        current_index = *it;
    }

    return current_index;
}

void delaunay::add_point(int point_index)
{
    auto point = m_dcel.vertex(point_index).point();
    // get node index that contains the point
    int node_index = find_node(point);
    auto position = get_position(point, node_index);
    assert(position != position::outside);

    if (position == position::strictly_interior) {
        // point is strictly in the triangle
        split_triangle_interior(point_index, node_index);
    }
    else {
        // point is on the triangle edge
        split_triangle_boundary(point_index, node_index);
    }
}

void delaunay::split_triangle_interior(int point_index, int node_index)
{
    // D is a new point
    //
    //         C            
    //         /\
    //        /  \
    //       /    \
    //      /  Dx  \
    //     /        \
    //  A ------------ B
    //
    //  get 3 new triangles -> BDA, CDB, ADC
    //  BDA is old face
    //  CDB is the first new face
    //  ADC is the second new face

    int face_id = m_graph[node_index].face();

    auto face_edge = m_dcel.face_edge(face_id);               // AB
    auto previous_edge = face_edge | edgerelation::previous;  // CA
    auto next_edge = face_edge | edgerelation::next;          // BC

    int new_edge_id = m_dcel.edge_count()+1;
    int new_face_id = m_dcel.face_count();

    // add a new edge: new_edge_id (DA)
    m_dcel.add(dceltype::edge{
        point_index+1,          // origin - D
        new_edge_id+5,          // twin - AD
        new_edge_id+1,          // previous - BD
        face_edge.id(),         // next - AB
        face_id});              // old face

    // add a new edge: new_edge_id+1 (BD)
    m_dcel.add(dceltype::edge{
        next_edge.origin(),     // origin - B
        new_edge_id+2,          // twin - DB
        face_edge.id(),         // previous - AB
        new_edge_id,            // next - DA
        face_id});              // old face

    // add a new edge: new_edge_id+2 (DB)
    m_dcel.add(dceltype::edge{
        point_index+1,          // origin - D
        new_edge_id+1,          // twin - BD
        new_edge_id+3,          // previous - CD
        next_edge.id(),         // next - BC
        new_face_id});          // the first new face

    // add a new edge: new_edge_id+3 (CD)
    m_dcel.add(dceltype::edge{
        previous_edge.origin(), // origin - C
        new_edge_id+4,          // twin - DC
        next_edge.id(),         // previous - BC
        new_edge_id+2,          // next - DB
        new_face_id});          // the first new face

    // add a new edge: new_edge_id+4 (DC)
    m_dcel.add(dceltype::edge{
        point_index+1,          // origin - D
        new_edge_id+3,          // twin - CD
        new_edge_id+5,          // previous - AD
        previous_edge.id(),     // next - CA
        new_face_id+1});        // the second new face

    // add a new edge: new_edge_id+5 (AD)
    m_dcel.add(dceltype::edge{
        face_edge.origin(),     // origin - A
        new_edge_id,            // twin - DA
        previous_edge.id(),     // previous - CA
        new_edge_id+4,          // next - DC
        new_face_id+1});        // the second new face

    // incident edge of D is DA
    m_dcel.vertex(point_index).set_incident_edge(new_edge_id);

    // update existing edges
    // AB previous is DA
    face_edge.set_previous(new_edge_id);
    // AB next is BD
    face_edge.set_next(new_edge_id+1);

    // BC previous is DB
    next_edge.set_previous(new_edge_id+2);
    // BC next is CD
    next_edge.set_next(new_edge_id+3);
    // BC face is the first new face
    next_edge.set_face(new_face_id);

    // CA previous is DC
    previous_edge.set_previous(new_edge_id+4);
    // CA next is AD
    previous_edge.set_next(new_edge_id+5);
    // CA face is the second new face
    previous_edge.set_face(new_face_id+1);

    // incident edge of the first new face is DB
    m_dcel.add(dceltype::face{new_edge_id+2});
    // incident edge of the second new face is DC
    m_dcel.add(dceltype::face{new_edge_id+4});

    // update the graph
    int new_node_id = m_graph.size();
    m_graph[node_index].set_children({new_node_id, new_node_id+1, new_node_id+2});

    // insert three new nodes
    auto node0 = node({
        (face_edge /*AB*/ | edgerelation::next /*BD*/).origin(),          // B
        (face_edge /*AB*/ | edgerelation::previous /*DA*/).origin(),      // D
        face_edge.origin()},                                              // A
        face_id);                                                         // old face

    auto node1 = node({
        (next_edge /*BC*/ | edgerelation::next /*CD*/).origin(),          // C
        (next_edge /*BC*/ | edgerelation::previous /*DB*/).origin(),      // D
        next_edge.origin()},                                              // B
        new_face_id);                                                     // the first new face

    auto node2 = node({
        (previous_edge /*CA*/ | edgerelation::next /*AD*/).origin(),      // A
        (previous_edge /*CA*/ | edgerelation::previous /*DC*/).origin(),  // D
        previous_edge.origin()},                                          // C
        new_face_id+1);                                                   // the second new face

    m_graph.add(std::move(node0));
    m_graph.add(std::move(node1));
    m_graph.add(std::move(node2));

    // flip edges if needed
    // there is a theorem that newly added edges cannot be flipped now, so we don't need to check them now
    // when one edge gets filpped, we will recursively check for other edges
    // AB
    try_flip(face_edge);
    // CA
    try_flip(previous_edge);
    // BC
    try_flip(next_edge);
}

void delaunay::split_triangle_boundary(int point_index, int node_index)
{
    // E is a new point
    //
    //  D ------- C     D ------- C
    //    |\    |         |\   /|
    //    | \   |         | \ / |
    //    |  x  |   -->   |  x  |
    //    | E \ |         | / \ |
    //    |    \|         |/   \|
    //  A ------- B     A ------- B
    //
    //  get 4 new triangles - ABE, EDA, CDE, EBC
    //  ABE - the first new face
    //  EDA - the first old face
    //  CDE - the second old face
    //  EBC - the second new face

    int new_edge_id = m_dcel.edge_count()+1;
    int new_face_id = m_dcel.face_count();

    // get edge id where new point is collinear
    int collinear_edge_id = m_dcel.collinear_edge_id(point_index, m_graph[node_index].face());
    assert(collinear_edge_id != -1);

    auto collinear_edge = m_dcel.edge(collinear_edge_id-1);     // BD
    auto collinear_edge2 = collinear_edge | edgerelation::twin; // DB

    int first_old_face = collinear_edge.face();
    int second_old_face = collinear_edge2.face();

    // nodes ids that will be updated
    int old_node1 = m_graph.get_node(first_old_face);  // ABD
    int old_node2 = m_graph.get_node(second_old_face); // BCD

    // ------------------------------- update ABD ----------------------------------------
    //  D ------- C     D ------- C
    //    |\    |         |\    |
    //    | \   |         | \   |
    //    |  x  |   -->   |  x  |
    //    | E \ |         | / \ |
    //    |    \|         |/   \|
    //  A ------- B     A ------- B
    auto prev_edge = collinear_edge | edgerelation::previous; // AB
    auto next_edge = collinear_edge | edgerelation::next;     // DA

    // add a new edge: new_edge_id (ED)
    m_dcel.add(dceltype::edge(
        point_index+1,          // origin - E
        collinear_edge2.id(),   // twin - DB (it will be updated to DE later)
        new_edge_id+1,          // previous - AE
        next_edge.id(),         // next - DA
        first_old_face));       // the first old face

    // add a new edge: new_edge_id+1 (AE)
    m_dcel.add(dceltype::edge(
        prev_edge.origin(),    // origin - A
        new_edge_id+2,         // twin - EA
        next_edge.id(),        // previous - DA
        new_edge_id,           // next - ED
        first_old_face));      // the first old face

    // add a new edge: new_edge_id+2 (EA)
    m_dcel.add(dceltype::edge(
        point_index+1,         // origin - E
        new_edge_id+1,         // twin - AE
        collinear_edge.id(),   // previous - BD (it will be updated to BE later)
        prev_edge.id(),        // next - AB
        new_face_id));         // the new face

    // incident edge of E is EA
    m_dcel.vertex(point_index).set_incident_edge(new_edge_id+2);

    // update BD to BE by changing next edge of BD to be EA
    collinear_edge.set_next(new_edge_id+2);
    // BE twin is EB
    collinear_edge.set_twin(new_edge_id+3);
    // BE belongs to new face now
    collinear_edge.set_face(new_face_id);

    // DA next is AE
    next_edge.set_next(new_edge_id+1);
    // DA previous is ED
    next_edge.set_previous(new_edge_id);

    // AB previous is EA
    prev_edge.set_previous(new_edge_id+2);
    // AB belongs to new face now
    prev_edge.set_face(new_face_id);

    // update the face
    // incident edge of the old face is AE
    m_dcel.face(first_old_face).set_edge(new_edge_id+1);

    // add new face
    // incident edge of the new face is EA
    m_dcel.add(dceltype::face(new_edge_id+2));

    // update the graph
    m_graph[old_node1].set_children({m_graph.size(), m_graph.size()+1});

    m_graph.add(node(
        {(collinear_edge /*BE*/ | edgerelation::previous /*AB*/).origin(),  // A
         collinear_edge.origin(),                                           // B
         (collinear_edge /*BE*/ | edgerelation::next /*EA*/).origin()},     // E
        new_face_id));                                                      // the first new face

    m_graph.add(node(
        {(next_edge /*DA*/ | edgerelation::previous /*ED*/).origin(),       // E
          next_edge.origin(),                                               // D
         (next_edge /*DA*/ | edgerelation::next /*AE*/).origin()},          // A
        first_old_face));                                                   // the first old face


    // ------------------------------- update BCD ----------------------------------------
    //  D ------- C     D ------- C
    //    |\    |         |\   /|
    //    | \   |         | \ / |
    //    |  x  |   -->   |  x  |
    //    | E \ |         |   \ |
    //    |    \|         |    \|
    //  A ------- B     A ------- B
    auto prev_edge2 = collinear_edge2 | edgerelation::previous; // CD
    auto next_edge2 = collinear_edge2 | edgerelation::next;     // BC

    // add a new edge: new_edge_id+3 (EB)
    m_dcel.add(dceltype::edge(
        point_index+1,        // origin - E
        collinear_edge.id(),  // twin - BE
        new_edge_id+4,        // previous - CE
        next_edge2.id(),      // next - BC
        new_face_id+1));      // the second new face

    // add a new edge: new_edge_id+4 (CE)
    m_dcel.add(dceltype::edge(
        prev_edge2.origin(),  // origin - C
        new_edge_id+5,        // twin - EC
        next_edge2.id(),      // previous - BC
        new_edge_id+3,        // next - EB
        new_face_id+1));      // the second new face

    // add a new edge: new_edge_id+5 (EC)
    m_dcel.add(dceltype::edge(
        point_index+1,        // origin - E
        new_edge_id+4,        // twin - CE
        collinear_edge2.id(), // previous - DB (it will be updated to DE later)
        prev_edge2.id(),      // next - CD
        second_old_face));    // the second old face

    // update DB to DE by changing next edge of DB to be EC
    collinear_edge2.set_next(new_edge_id+5);
    // DE twin is ED
    collinear_edge2.set_twin(new_edge_id);

    //  BC previous is EB
    next_edge2.set_previous(new_edge_id+3);
    //  BC next is CE
    next_edge2.set_next(new_edge_id+4);
    //  BC face is the second new face
    next_edge2.set_face(new_face_id+1);

    // CD prevous is EC
    prev_edge2.set_previous(new_edge_id+5);

    // update the face
    // edge in the second old face is EC 
    m_dcel.face(second_old_face).set_edge(new_edge_id+5);

    // add new face
    // edge in the second new face is CE
    m_dcel.add(dceltype::face{new_edge_id+4});

    // update the graph
    m_graph[old_node2].set_children({m_graph.size(), m_graph.size()+1});

    m_graph.add(node(
        {(collinear_edge2 /*DE*/ | edgerelation::previous /*CD*/).origin(), // C
         collinear_edge2.origin(),                                          // D
         (collinear_edge2 /*DE*/ | edgerelation::next /*EC*/).origin()},    // E
        second_old_face));                                                  // the second old face

    m_graph.add(node(
        {(next_edge2 /*BC*/ | edgerelation::previous /*EB*/).origin(),      // E
         next_edge2.origin(),                                               // B
         (next_edge2 /*BC*/ | edgerelation::next /*CE*/).origin()},         // C
        new_face_id+1));                                                    // the second new face

    // flip edges if needed
    // there is a theorem that newly added edges cannot be flipped now, so we don't need to check them now
    // when one edge gets filpped, we will recursively check for other edges
    // AB
    try_flip(prev_edge);
    // DA
    try_flip(next_edge);
    // CD
    try_flip(prev_edge2);
    // BC
    try_flip(next_edge2);
}

void delaunay::try_flip(dcel::edgeref<false> edge)
{
    if (edge.external()) {
        // external edge cannot be flipped
        return;
    }

    //        C
    //       /|\
    //      / | \
    //     /  |  \
    //  D /   |   \ B
    //    \   |   /
    //     \  |  /
    //      \ | /
    //       \|/
    //        A
    //
    // edge - AC

    bool flip = false;

    int a = edge.origin();
    int c = (edge | edgerelation::twin).origin();
    int d = (edge | edgerelation::previous).origin();
    int b = (edge | edgerelation::twin | edgerelation::previous).origin();

    if (edge.has_negative_vertex()) {
        // A or C are negative (point_minus_2 or point_minus_1)
        // cannot both of them be negative since point_minus_2 - point_minus_1 is an external edge
        if (d > 0 &&  b > 0) {
            auto d_point = (edge | edgerelation::previous).point();
            auto b_point = (edge | edgerelation::twin | edgerelation::previous).point();

            // set d_point to have highest y-coordinate
            if (d_point.y() < b_point.y()) {
                std::swap(d_point, b_point);
            }

            if (a < 0) {
                // A is negative -> C is positive
                auto c_point = (edge | edgerelation::twin).point();

                flip = (a == dceltype::point_minus_2) ? 
                    //        d
                    //         |\   c
                    //         | \ /
                    //         |  \
                    //         | / \
                    // p_minus_2----b
                    // don't let p_minus_2-c to go inside convex hull, 
                    // flip p_minus_2-c to get b-d
                    (d_point.get_direction(b_point, c_point) == util::direction::positive) : 
                    // don't let p_minus_1-c to go inside convex hull, 
                    // flip p_minus_1-c to get b-d
                    (d_point.get_direction(b_point, c_point) == util::direction::negative);
            }
            else {
                // C is negative -> A is positive
                auto a_point = edge.point();

                flip = (c == dceltype::point_minus_2) ?
                    // don't let p_minus_2-a to go inside convex hull, 
                    // flip p_minus_2-a to get b-d
                    (d_point.get_direction(b_point, a_point) == util::direction::positive) :
                    // don't let p_minus_1-a to go inside convex hull, 
                    // flip p_minus_1-a to get b-d
                    (d_point.get_direction(b_point, a_point) == util::direction::negative);
            }
        }
    }
    else {
        // A and C are positive
        if (d > 0 && b > 0) {
            // all points are positive - in circle check
            auto a_point = edge.point();
            auto c_point = (edge /*AC*/ | edgerelation::twin /*CA*/).point();
            auto d_point = (edge | edgerelation::previous).point();
            auto b_point = (edge | edgerelation::twin | edgerelation::previous).point();

            // check if b is strictly in the circle defined by a-c-d
            flip = b_point.in_circle(a_point,c_point,d_point);
        }
    }

    if (flip) {
        flip_edge(edge);
    }
}

void delaunay::flip_edge(dcel::edgeref<false> edge)
{
    //        C                        C
    //       /|\                      / \
    //      / | \                    /   \
    //     /  |  \                  /     \
    //  D /   |   \ B    -->     D /-------\ B
    //    \   |   /                \       /
    //     \  |  /                  \     /
    //      \ | /                    \   /
    //       \|/                      \ /
    //        A                        A
    //
    // edge - AC
    
    auto twin = edge | edgerelation::twin; // CA

    int old_node1 = m_graph.get_node(edge.face());
    int old_node2 = m_graph.get_node(twin.face());

    if (edge.origin() > 0) {
        // if A is positive update its incident edge to be AB
        m_dcel.vertex(edge.origin()-1).set_incident_edge(twin.next());
    }
    if ((edge | edgerelation::twin).origin() > 0) {
        // if C is positive update its incident edge to be CD
        m_dcel.vertex((edge | edgerelation::twin).origin()-1).set_incident_edge(edge.next());
    }

    // AC will become BD
    // CA will become DB
    // for now, set only origins
    auto temp = (edge | edgerelation::previous).origin(); // D
    // AC origin is B
    edge.set_origin((twin | edgerelation::previous).origin());
    // CA origin is D
    twin.set_origin(temp);

    // update next edges
    // CD next is DA (it will be DB)
    (edge | edgerelation::next).set_next(edge.twin());
    // AB next is BC (it will be BD)
    (twin | edgerelation::next).set_next(edge.id());
    // DA next is AB
    (edge | edgerelation::previous).set_next(twin.next());
    // BC next is CD
    (twin | edgerelation::previous).set_next(edge.next());
    // BC next is DA, so it is BD now
    edge.set_next(edge.previous());
    // DA next is BC, so it is DB now
    twin.set_next(twin.previous());

    // update previous
    // BD previous is AB
    edge.set_previous((edge | edgerelation::next).next());
    // DB previous is CD
    twin.set_previous((twin | edgerelation::next).next());
    // DA previous is BD
    (edge | edgerelation::next).set_previous(edge.id());
    // BC previous is DB
    (twin | edgerelation::next).set_previous(twin.id());
    // AB previous is DA
    (edge | edgerelation::previous).set_previous(edge.next());
    // CD previous is BC
    (twin | edgerelation::previous).set_previous(twin.next());

    // update faces
    //
    //        C                        C
    //       /|\                      / \
    //      / | \                    /   \
    //     /  |  \                  /  2  \
    //  D / 1 | 2 \ B    -->     D /-------\ B
    //    \   |   /                \       /
    //     \  |  /                  \  1  /
    //      \ | /                    \   /
    //       \|/                      \ /
    //        A                        A
    //
    // two faces - 1 and 2

    // AB face is face 1
    (edge | edgerelation::previous).set_face(edge.face());
    // CD face is face 2
    (twin | edgerelation::previous).set_face(twin.face());

    // incident edge of face 1 is BD
    m_dcel.face(edge.face()).set_edge(edge.id());
    // incident edge of face 2 is DB
    m_dcel.face(twin.face()).set_edge(edge.twin());

    // update the graph
    m_graph[old_node1].set_children({m_graph.size(), m_graph.size()+1});
    m_graph[old_node2].set_children({m_graph.size(), m_graph.size()+1});

    // insert two new nodes
    m_graph.add(node({
        (edge | edgerelation::previous).origin(),  // A
        edge.origin(),                             // B
        (edge | edgerelation::next).origin()},     // D
        edge.face()));                             // face 1

    m_graph.add(node(
        {(twin | edgerelation::previous).origin(), // C
        (edge | edgerelation::twin).origin(),      // D
        (twin | edgerelation::next).origin()},     // B
        twin.face()));                             // face 2

    // recursively check edges that can be illegal
    // AB
    try_flip(edge | edgerelation::previous);
    // BC
    try_flip(twin | edgerelation::next);
}

std::vector<double> delaunay::range() const
{
    auto xmin_it = std::min_element(m_dcel.vertices().begin(), m_dcel.vertices().end(), 
            [](const auto& a, const auto& b) { return a.point().x() < b.point().x(); });
    auto xmax_it = std::max_element(m_dcel.vertices().begin(), m_dcel.vertices().end(), 
            [](const auto& a, const auto& b) { return a.point().x() < b.point().x(); });
    auto ymin_it = std::min_element(m_dcel.vertices().begin(), m_dcel.vertices().end(), 
            [](const auto& a, const auto& b) { return a.point().y() < b.point().y(); });
    auto ymax_it = std::max_element(m_dcel.vertices().begin(), m_dcel.vertices().end(), 
            [](const auto& a, const auto& b) { return a.point().y() < b.point().y(); });

    return {xmin_it->point().x(), xmax_it->point().x(), ymin_it->point().y(), ymax_it->point().y()};
}

std::vector<util::line_segment> delaunay::get_edges() const
{
    std::vector<util::line_segment> result;
    // face 0 is external face
    for (int face_index = 1; face_index < m_dcel.face_count(); ++face_index) {
        auto current_edge = m_dcel.edge(m_dcel.face(face_index).edge()-1);
        int start_id = current_edge.id();

        bool valid_face = true;
        std::vector<util::line_segment> edges;

        do {
            if (current_edge.origin() > 0 && (current_edge | edgerelation::twin).origin() > 0) {
                edges.emplace_back(current_edge.point(), (current_edge | edgerelation::twin).point());
            }
            else {
                valid_face = false;
            }
            current_edge = current_edge | edgerelation::next;
        } while (start_id != current_edge.id() && valid_face);

        if (valid_face) {
            result.insert(result.end(), std::make_move_iterator(edges.begin()), std::make_move_iterator(edges.end()));
        }
    }

    return result;
}


const dcel& delaunay::triangulation() const
{
    return m_dcel;
}

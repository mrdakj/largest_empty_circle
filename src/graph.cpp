#include "graph.h"
#include <algorithm>
#include <cassert>
#include <iostream>

// node
node::node(std::vector<int> vertices, int face_id)
    : m_vertices(std::move(vertices))
    , m_face(face_id)
{}

void node::set_children(std::vector<int> children)
{
    m_children = std::move(children);
}

int node::children_count() const
{
    return m_children.size();
}

int node::face() const
{
    return m_face;
}

bool node::leaf() const
{
    return children_count() == 0;
}

const std::vector<int>& node::children() const
{
    return m_children;
}

const std::vector<int>& node::vertices() const
{
    return m_vertices;
}

// graph
int graph::size() const
{
    return m_nodes.size();
}

const node& graph::operator[](int i) const
{
    return m_nodes[i];
}

node& graph::operator[](int i)
{
    return m_nodes[i];
}

const std::vector<node>& graph::nodes() const
{
    return m_nodes;
}

void graph::add(node n)
{
    if (n.face() > 0)
    {
        m_face_to_node[n.face()] = m_nodes.size();
    }
    m_nodes.emplace_back(std::move(n));
}

int graph::get_node(int face_id) const
{
    auto it = m_face_to_node.find(face_id); 
    assert(it != m_face_to_node.end());
    return it->second;
}

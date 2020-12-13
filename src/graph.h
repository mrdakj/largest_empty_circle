#ifndef GRAPH_H
#define GRAPH_H 

#include <vector>
#include <unordered_map>

class node {
public:
    node(std::vector<int> vertices, int face_id);

    const std::vector<int>& children() const;
    const std::vector<int>& vertices() const;
    int face() const;

    void set_children(std::vector<int> children);
    int children_count() const;
    bool leaf() const;

private:
    // vertices ids
    std::vector<int> m_vertices;
    // children nodes ids
    std::vector<int> m_children;
    // face id
    int m_face;
};

class graph
{
public:
    const std::vector<node>& nodes() const;
    const node& operator[](int i) const;
    node& operator[](int i);

    // get node id for the given face
    int get_node(int face_id) const;

    void add(node n);
    int size() const;

private:
    // map: face id -> node id
    // last node id with face id
    std::unordered_map<int,int> m_face_to_node;
    std::vector<node> m_nodes;
};

#endif /* GRAPH_H */

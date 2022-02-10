
#include <string>

#include "Graph.hpp"


namespace raven {

Graph::Node::Node(const biosoup::NucleicAcid& sequence) : id(num_objects++),
                                                          sequence(sequence),
                                                          count(1),
                                                          is_unitig(),
                                                          is_circular(),
                                                          is_polished(),
                                                          transitive(),
                                                          inedges(),
                                                          outedges(),
                                                          pair() {}


Graph::Node::Node(Graph::Node* begin, Graph::Node* end) : id(num_objects++),
                                                          sequence(),
                                                          count(),
                                                          is_unitig(),
                                                          is_circular(begin == end),
                                                          is_polished(),
                                                          transitive(),
                                                          inedges(),
                                                          outedges(),
                                                          pair() {
    std::string data{};

    auto it = begin;
    while (true) {
        data += it->outedges.front()->Label();
        count += it->count;
        if ((it = it->outedges.front()->head) == end) {
            break;
        }
    }
    if (begin != end) {
        data += end->sequence.InflateData();
        count += end->count;
    }

    is_unitig = count > 5 && data.size() > 9999;

    sequence = biosoup::NucleicAcid(
            (is_unitig ? "Utg" : "Ctg") + std::to_string(id & (~1UL)),
            data);
}

Graph::Edge::Edge(Graph::Node* tail, Graph::Node* head, std::uint32_t length) : id(num_objects++), length(length), weight(0), tail(tail),
                                                                                head(head), pair() {
    tail->outedges.emplace_back(this);
    head->inedges.emplace_back(this);
}

Graph::Graph() : stage(-5), piles(), nodes(), edges() {

    Graph::Node::num_objects = 0;
    Graph::Edge::num_objects = 0;
}

std::atomic<std::uint32_t> Graph::Node::num_objects{0};
std::atomic<std::uint32_t> Graph::Edge::num_objects{0};

}

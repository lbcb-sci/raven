#include <map>
#include "raven/graph/serialization/graph_repr.h"

namespace raven {

void PrintGfa(const Graph& graph, const std::string& path) {
  if (path.empty()) {
    return;
  }

  std::ofstream os(path);
  for (const auto& it : graph.nodes) {
    if ((it == nullptr) || it->is_rc() ||
        (it->count == 1 && it->outdegree() == 0 && it->indegree() == 0)) {
      continue;
    }
    os << "S\t" << it->sequence.name << "\t" << it->sequence.InflateData()
       << "\tLN:i:" << it->sequence.inflated_len << "\tRC:i:" << it->count
       << std::endl;
    if (it->is_circular) {
      os << "L\t" << it->sequence.name << "\t" << '+' << "\t"
         << it->sequence.name << "\t" << '+' << "\t0M" << std::endl;
    }
  }
  for (const auto& it : graph.edges) {
    if (it == nullptr || it->is_rc()) {
      continue;
    }
    os << "L\t" << it->tail->sequence.name << "\t"
       << (it->tail->is_rc() ? '-' : '+')  // NOLINT
       << "\t" << it->head->sequence.name << "\t"
       << (it->head->is_rc() ? '-' : '+')  // NOLINT
       << "\t" << it->tail->sequence.inflated_len - it->length << 'M'
       << std::endl;
  }
  os.close();
}

void PrintCsv(const Graph& graph, const std::string& path) {
  if (path.empty()) {
    return;
  }

  std::ofstream os(path);
  for (const auto& it : graph.nodes) {
    if ((it == nullptr) || it->is_rc() ||
        (it->count == 1 && it->outdegree() == 0 && it->indegree() == 0)) {
      continue;
    }
    os << it->id << " [" << it->id / 2 << "]"
       << " LN:i:" << it->sequence.inflated_len << " RC:i:" << it->count << ","
       << it->pair->id << " [" << it->pair->id / 2 << "]"
       << " LN:i:" << it->pair->sequence.inflated_len
       << " RC:i:" << it->pair->count << ",0,-" << std::endl;
  }
  for (const auto& it : graph.edges) {
    if (it == nullptr) {
      continue;
    }
    os << it->tail->id << " [" << it->tail->id / 2 << "]"
       << " LN:i:" << it->tail->sequence.inflated_len
       << " RC:i:" << it->tail->count << "," << it->head->id << " ["
       << it->head->id / 2 << "]"
       << " LN:i:" << it->head->sequence.inflated_len
       << " RC:i:" << it->head->count << ",1," << it->id << " " << it->length
       << " " << it->weight << std::endl;
  }
  for (const auto& it : graph.nodes) {  // circular edges TODO(rvaser): check
    if (it == nullptr || !it->is_circular) {
      continue;
    }
    os << it->id << " [" << it->id / 2 << "]"
       << " LN:i:" << it->sequence.inflated_len << " RC:i:" << it->count << ","
       << it->id << " [" << it->id / 2 << "]"
       << " LN:i:" << it->sequence.inflated_len << " RC:i:" << it->count
       << ",1,-" << std::endl;
  }
  os.close();
}

void PrintJson(const raven::Graph& graph, const std::string& path) {
  if (path.empty()) {
    return;
  }

  std::ofstream os(path);
  cereal::JSONOutputArchive archive(os);

  for (const auto& it : graph.piles) {
    if (it->is_invalid()) {
      continue;
    }

    archive(cereal::make_nvp(std::to_string(it->id()), *(it.get())));
  }
}

void splitString(const std::string &inputString, char delimiter, std::vector<std::string> &elems) {
    std::stringstream stringStream;
    stringStream.str(inputString);
    std::string item;
    while (std::getline(stringStream, item, delimiter)) {
        elems.push_back(item);
    }
}

Node* findNodeForSequenceName(std::string   tailSequenceName, std::vector<std::unique_ptr<Node>> &nodes) {
  for (auto &node : nodes) {
    if (node == nullptr) continue;
    if (node->sequence.name == tailSequenceName) return node.get();
  }
  return nullptr;
}

Graph LoadGfa(const std::string& path) {
  Graph graph = Graph();
  
  if (path.empty()) {
    return graph;
  }

  std::uint32_t currentNodeEvenId = 0;
  std::uint32_t currentNodeOddId  = 1;

  std::uint32_t currentEdgeId     = 0;

  std::ifstream is(path);
  std::string inputLine;
  
  // TODO (adolmac) - maybe it would be a good idea to read all nodes first, then read all edges, since I need all nodes in order to create edges
  while(getline(is, inputLine)) {

    std::cout << "Input line: " << inputLine << std::endl;

    std::vector<std::string> rowValues;
    splitString(inputLine, '\t', rowValues);

    if (rowValues[0] == "S") { // this is a node
      
      std::cout << "Creating new node" << std::endl;

      std::string   sequenceName                = rowValues[1];
      std::string   sequenceInflatedData        = rowValues[2];
      std::uint32_t sequenceInflatedDataLength  = stol(rowValues[3].substr(5));
      std::uint32_t count                       = stol(rowValues[4].substr(5));

      biosoup::NucleicAcid sequence = biosoup::NucleicAcid(sequenceName, sequenceInflatedData);
      std::unique_ptr<Node> newNode(new Node());
      
      newNode->count       = count;
      newNode->is_unitig   = false;
      newNode->is_circular = false;
      newNode->is_polished = false;
      newNode->sequence    = sequence;
      
      graph.nodes.push_back(std::move(newNode));

    } else if (rowValues[0] == "L") { // this is an egde

      std::cout << "Creating new edge" << std::endl;

      std::string   tailSequenceName                  = rowValues[1];
      std::string   isTailReverseComplement           = rowValues[2];
      std::string   headSequenceName                  = rowValues[3];
      std::string   isHeadReverseComplement           = rowValues[4];
      std::uint32_t tailInflatedLengthMinusEdgeLength = stol(rowValues[5].substr(0, rowValues[5].size() - 1));

      Node* tail;
      std::uint32_t edgeLength = 0;
      tail = findNodeForSequenceName(tailSequenceName, graph.nodes);
      if (tail != nullptr) {
        if (isTailReverseComplement == "+") {
          tail->id = currentNodeEvenId++;
        } else {
          tail->id = currentNodeOddId++;
        }
        edgeLength = tail->sequence.inflated_len - tailInflatedLengthMinusEdgeLength;
      }

      Node* head;
      head = findNodeForSequenceName(headSequenceName, graph.nodes);
      if (tail != nullptr) {
        if (isHeadReverseComplement == "+") {
          head->id = currentNodeEvenId++;
        } else {
          head->id = currentNodeOddId++;
        }
      }

      std::unique_ptr<Edge> newEdge(new Edge(tail, head, edgeLength));
      newEdge->id = currentEdgeId; // (adolmac) have no idea if this is Ok or not since I do not have info about this in GFA

      graph.edges.push_back(std::move(newEdge));

    } else {
      std::cout << "Unknown element: " << inputLine << std::endl;
    }
  
  }

  graph.stage = -3;

  is.close();

  return graph;
}

}  // namespace raven

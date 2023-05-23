#include "raven/graph/serialization/graph_repr.h"
#include <string>
#include <vector>
#include "edlib.h"
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

    os << "S"; 
    os << "\t"; 
    os << it->sequence.name;
    os << "\t";
    os << it->sequence.InflateData();
    os << "\t"; 
    os << "LN:i:" << it->sequence.inflated_len;
    os << "\t";
    os << "RC:i:" << it->count;
    os << "\t";
    os << "dp:f:" << it->coverage;
    os << std::endl;
    
    if (it->is_circular) {
      os << "L";
      os << "\t";
      os << it->sequence.name;
      os << "\t";
      os << '+';
      os << "\t";
      os << it->sequence.name;
      os << "\t";
      os << '+';
      os << "\t";
      os << "0M"; 
      os << std::endl;
    }
  }

  for (const auto& it : graph.edges) {
    if (it == nullptr || it->is_rc()) {
      continue;
    }
    
    os << "L";
    os << "\t";
    os << it->tail->sequence.name;
    os << "\t";
    os << (it->tail->is_rc() ? '-' : '+');  // NOLINT
    os << "\t";
    os << it->head->sequence.name;
    os << "\t";
    os << (it->head->is_rc() ? '-' : '+');  // NOLINT
    os << "\t";
    os << it->tail->sequence.inflated_len - it->length << 'M';
    os << std::endl;
  }

  os.close();
}

std::vector<std::string> getGfa(const Graph& graph) {
  std::vector<std::string> resultVector;
  std::string line;

  for (const auto& it : graph.nodes) {
    if ((it == nullptr) || it->is_rc() ||
        (it->count == 1 && it->outdegree() == 0 && it->indegree() == 0)) {      
      continue;
    }

    line.clear();

    line += "S"; 
    line += "\t"; 
    line += it->sequence.name;
    line += "\t";
    line += it->sequence.InflateData();
    line += "\t"; 
    line += "LN:i:" + std::to_string(it->sequence.inflated_len);
    line += "\t";
    line += "RC:i:" + std::to_string(it->count);
    
    resultVector.push_back(line);

    if (it->is_circular) {
      line.clear();

      line += "L";
      line += "\t";
      line += it->sequence.name;
      line += "\t";
      line += '+';
      line += "\t";
      line += it->sequence.name;
      line += "\t";
      line += '+';
      line += "\t";
      line += "0M"; 
      
      resultVector.push_back(line);
    }
  }

  for (const auto& it : graph.edges) {
    if (it == nullptr || it->is_rc()) {
      continue;
    }
    line.clear();

    line += "L";
    line += "\t";
    line += it->tail->sequence.name;
    line += "\t";
    line += (it->tail->is_rc() ? '-' : '+');  // NOLINT
    line += "\t";
    line += it->head->sequence.name;
    line += "\t";
    line += (it->head->is_rc() ? '-' : '+');  // NOLINT
    line += "\t";
    line += std::to_string(it->tail->sequence.inflated_len - it->length) + 'M';

    resultVector.push_back(line);
  }

  return resultVector;
}

void PrintCsv(const Graph& graph, const std::string& path, bool printSequenceName, bool printPileBeginEnd, bool printEdgeSimilarity) {
  if (path.empty()) {
    return;
  }

  std::ofstream os(path);

  for (const auto& it : graph.nodes) {
    if ((it == nullptr) || it->is_rc() ||
        (it->count == 1 && it->outdegree() == 0 && it->indegree() == 0)) {
      continue;
    }

    bool addDashAtTheEnd = true;

    os << it->id;
    os << " [" << it->id / 2 << "]";
    os << " LN:i:" << it->sequence.inflated_len;
    os << " RC:i:" << it->count;
    os << ",";
    os << it->pair->id;
    os << " [" << it->pair->id / 2 << "]";
    os << " LN:i:" << it->pair->sequence.inflated_len;
    os << " RC:i:" << it->pair->count;
    os << ",";
    os << "0";
    os << ",";

    if (printSequenceName) {
      addDashAtTheEnd = false;
      os << it->sequence.name;
      os << " ";
    }

    if (printPileBeginEnd && it->sequence.id < graph.piles.size()) {
      addDashAtTheEnd = false;
      os << graph.piles[it->sequence.id]->begin();
      os << " ";
      os << graph.piles[it->sequence.id]->end();
    } 

    if (addDashAtTheEnd) {
      os << "-";
    }

    os << std::endl;
  }

  for (const auto& it : graph.edges) {
    if (it == nullptr) {
      continue;
    }

    os << it->tail->id;
    os << " [" << it->tail->id / 2 << "]";
    os << " LN:i:" << it->tail->sequence.inflated_len;
    os << " RC:i:" << it->tail->count;
    os << ",";
    os << it->head->id;
    os << " [" << it->head->id / 2 << "]";
    os << " LN:i:" << it->head->sequence.inflated_len;
    os << " RC:i:" << it->head->count;
    os << ",";
    os << "1";
    os << ",";
    os << it->id;
    os << " ";
    os << it->length;
    os << " ";
    os << it->weight;

    if (printEdgeSimilarity) {
      std::string lhs{it->tail->sequence.InflateData(it->length)};
      std::string rhs{it->head->sequence.InflateData(0, lhs.size())};
      EdlibAlignResult result = edlibAlign(
        lhs.c_str(), lhs.size(),
        rhs.c_str(), rhs.size(),
        edlibDefaultAlignConfig());
      double score = 1 - result.editDistance / static_cast<double>(lhs.size());

      os << " ";
      os << score;
    }

    os << std::endl;
  }

  for (const auto& it : graph.nodes) {  // circular edges TODO(rvaser): check
    if (it == nullptr || !it->is_circular) {
      continue;
    }
    os << it->id;
    os << " [" << it->id / 2 << "]";
    os << " LN:i:" << it->sequence.inflated_len;
    os << " RC:i:" << it->count;
    os << ",";
    os << it->id;
    os << " [" << it->id / 2 << "]";
    os << " LN:i:" << it->sequence.inflated_len;
    os << " RC:i:" << it->count;
    os << ",";
    os << "1";
    os << ",";
    os << "-";
    os << std::endl;
  }

  os.close();
}

std::vector<std::string> getCsv(const Graph& graph, bool printSequenceName, bool printPileBeginEnd, bool printEdgeSimilarity) {
  std::vector<std::string> resultVector;
  std::string line;  

  
  for (const auto& it : graph.nodes) {
    if ((it == nullptr) || it->is_rc() ||
        (it->count == 1 && it->outdegree() == 0 && it->indegree() == 0)) {
      continue;
    }

    line.clear();
    bool addDashAtTheEnd = true;

    line += std::to_string(it->id);
    line += " [" + std::to_string(it->id / 2) + "]";
    line += " LN:i:" + std::to_string(it->sequence.inflated_len);
    line += " RC:i:" + std::to_string(it->count);
    line += ",";
    line += std::to_string(it->pair->id);
    line += " [" + std::to_string(it->pair->id / 2) + "]";
    line += " LN:i:" + std::to_string(it->pair->sequence.inflated_len);
    line += " RC:i:" + std::to_string(it->pair->count);
    line += ",";
    line += "0";
    line += ",";

    if (printSequenceName) {
      addDashAtTheEnd = false;
      line += it->sequence.name;
      line += " ";
    }

    if (printPileBeginEnd && it->sequence.id < graph.piles.size()) {
      addDashAtTheEnd = false;
      line += std::to_string(graph.piles[it->sequence.id]->begin());
      line += " ";
      line += std::to_string(graph.piles[it->sequence.id]->end());
    } 

    if (addDashAtTheEnd) {
      line += "-";
    }

    resultVector.push_back(line);
  }

  for (const auto& it : graph.edges) {
    if (it == nullptr) {
      continue;
    }

    line.clear();

    line += std::to_string(it->tail->id);
    line += " [" + std::to_string(it->tail->id / 2) + "]";
    line += " LN:i:" + std::to_string(it->tail->sequence.inflated_len);
    line += " RC:i:" + std::to_string(it->tail->count);
    line += ",";
    line += std::to_string(it->head->id);
    line += " [" + std::to_string(it->head->id / 2) + "]";
    line += " LN:i:" + std::to_string(it->head->sequence.inflated_len);
    line += " RC:i:" + std::to_string(it->head->count);
    line += ",";
    line += "1";
    line += ",";
    line += std::to_string(it->id);
    line += " ";
    line += std::to_string(it->length);
    line += " ";
    line += std::to_string(it->weight);

    if (printEdgeSimilarity) {
      std::string lhs{it->tail->sequence.InflateData(it->length)};
      std::string rhs{it->head->sequence.InflateData(0, lhs.size())};
      EdlibAlignResult result = edlibAlign(
        lhs.c_str(), lhs.size(),
        rhs.c_str(), rhs.size(),
        edlibDefaultAlignConfig());
      double score = 1 - result.editDistance / static_cast<double>(lhs.size());

      line += " ";
      line += std::to_string(score);
    }

    resultVector.push_back(line);
  }

  for (const auto& it : graph.nodes) {  // circular edges TODO(rvaser): check
    if (it == nullptr || !it->is_circular) {
      continue;
    }

    line.clear();

    line += std::to_string(it->id);
    line += " [" + std::to_string(it->id / 2) + "]";
    line += " LN:i:" + std::to_string(it->sequence.inflated_len);
    line += " RC:i:" + std::to_string(it->count);
    line += ",";
    line += std::to_string(it->id);
    line += " [" + std::to_string(it->id / 2) + "]";
    line += " LN:i:" + std::to_string(it->sequence.inflated_len);
    line += " RC:i:" + std::to_string(it->count);
    line += ",";
    line += "1";
    line += ",";
    line += "-";

    resultVector.push_back(line);
  }

  return resultVector;
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

  std::uint32_t currentNodeId = 0;
  std::uint32_t currentEdgeId = 0;

  std::ifstream is(path);
  std::string inputLine;
  
  while(getline(is, inputLine)) {

    std::vector<std::string> rowValues;
    splitString(inputLine, '\t', rowValues);

    if (rowValues[0] == "S") { // this is a node
      
      std::string   sequenceName                = rowValues[1];
      std::string   sequenceInflatedData        = rowValues[2];
      std::uint32_t count                       = stol(rowValues[4].substr(5));

      biosoup::NucleicAcid sequence = biosoup::NucleicAcid(sequenceName, sequenceInflatedData);
      std::unique_ptr<Node> newNode(new Node());
      
      newNode->id          = currentNodeId;
      currentNodeId += 2; // since node is_rc() method is based on node id and PrintGfa only prints !is_rc() Nodes all ids are set to even numbers

      newNode->count       = count;
      newNode->is_unitig   = false;
      newNode->is_circular = false;
      newNode->is_polished = false;
      newNode->sequence    = sequence;
      
      graph.nodes.push_back(std::move(newNode));

    } else if (rowValues[0] == "L") { // this is an egde

      std::string   tailSequenceName                  = rowValues[1];
      std::string   isTailReverseComplement           = rowValues[2];
      std::string   headSequenceName                  = rowValues[3];
      std::string   isHeadReverseComplement           = rowValues[4];
      std::uint32_t tailInflatedLengthMinusEdgeLength = stol(rowValues[5].substr(0, rowValues[5].size() - 1));

      if (tailInflatedLengthMinusEdgeLength == 0 && tailSequenceName == headSequenceName) { // circular

        Node* node;
        node = findNodeForSequenceName(headSequenceName, graph.nodes);
        if (node != nullptr) {
          node->is_circular = true;
        }

      } else {

        Node* tail;
        std::uint32_t edgeLength = 0;
        tail = findNodeForSequenceName(tailSequenceName, graph.nodes);
        if (tail != nullptr) {
          edgeLength = tail->sequence.inflated_len - tailInflatedLengthMinusEdgeLength;
        }

        Node* head;
        head = findNodeForSequenceName(headSequenceName, graph.nodes);

        std::unique_ptr<Edge> newEdge(new Edge(tail, head, edgeLength));
        newEdge->id = currentEdgeId; // (adolmac) have no idea if this is Ok or not since I do not have info about egde Id in GFA
        currentEdgeId += 2;          // since edge is_rc() method is based on edge id and PrintGfa only prints !is_rc() Edges all ids are set to even numbers

        if (tail != nullptr) {
          tail->outedges.push_back(newEdge.get());
        }

        if (head != nullptr) {
          head->inedges.push_back(newEdge.get());
        }

        graph.edges.push_back(std::move(newEdge));

      }


    } else {
      std::cout << "Unknown element: " << inputLine << std::endl;
    }
  
  }

  graph.stage = -3;

  is.close();

  return graph;
}

}  // namespace raven

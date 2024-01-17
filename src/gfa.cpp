#include "centrolign/gfa.hpp"

#include <stdexcept>

namespace centrolign {

using namespace std;

BaseGraph read_gfa(std::istream& in, bool encode,
                   std::unordered_map<std::string, int>* int_tags) {
    
    // records of (first node ID, last node ID)
    vector<pair<uint64_t, uint64_t>> decompacted_nodes;
    
    BaseGraph graph;
    
    string line;
    while (in) {
        getline(in, line);
    
        if (line.empty()) {
            continue;
        }
        auto tokens = tokenize(line);
        if (tokens.front() == "S") {
            if (tokens.size() != 3) {
                throw runtime_error("GFA S line does not have 3 tab-separated columns: " + line);
            }
            // make a non-branching path of nodes
            auto gfa_id = parse_int(tokens[1]);
            while (decompacted_nodes.size() <= gfa_id) {
                decompacted_nodes.emplace_back(-1, -1);
            }
            auto& seq = tokens[2];
            if (encode) {
                seq = encode_seq(seq);
            }
            auto node_id = graph.add_node(seq.front());
            decompacted_nodes[gfa_id] = make_pair(node_id, node_id);
            for (size_t i = 1; i < seq.size(); ++i) {
                auto next_id = graph.add_node(seq[i]);
                decompacted_nodes[gfa_id].second = next_id;
                graph.add_edge(node_id, next_id);
                node_id = next_id;
            }
        }
        else if (tokens.front() == "L") {
            if (tokens.size() != 6) {
                throw runtime_error("GFA L line does not have 6 tab-separated columns: " + line);
            }
            if (tokens[5] != "*" && tokens[5] != "0M") {
                throw runtime_error("GFA cannot contain overlapped edges: " + line);
            }
            if (tokens[2] != tokens[4]) {
                throw runtime_error("GFA cannot contain reversing edges: " + line);
            }
            else if (tokens[2] == "-") {
                // forward-ize the edge
                swap(tokens[1], tokens[3]);
            }
            // add an edge between the endpoints of the nodes' paths
            graph.add_edge(decompacted_nodes[parse_int(tokens[1])].second,
                           decompacted_nodes[parse_int(tokens[3])].first);
        }
        else if (tokens.front() == "P") {
            if (tokens.size() != 4) {
                throw runtime_error("GFA P line does not have 4 tab-separated columns: " + line);
            }
            // FIXME: it could also contain 0M,0M,...
            if (tokens[3] != "*") {
                throw runtime_error("GFA cannot contain overlapped paths: " + line);
            }
            // add a path
            auto path_id = graph.add_path(tokens[1]);
            auto path_nodes = tokenize(tokens[2], ',');
            for (auto& node : path_nodes) {
                if (node.back() != '+') {
                    throw runtime_error("GFA path " + tokens[1] + " has reverse strand node: " + node);
                }
                auto gfa_id = parse_int(node.substr(0, node.size() - 1));
                auto node_id = decompacted_nodes[gfa_id].first;
                graph.extend_path(path_id, node_id);
                while (node_id != decompacted_nodes[gfa_id].second) {
                    node_id = graph.next(node_id).front();
                    graph.extend_path(path_id, node_id);
                }
            }
        }
        else if (tokens.front() == "H") {
            // parse any tags from the header
            for (size_t i = 1; i < tokens.size(); ++i) {
                auto tag_tokens = tokenize(tokens[i], ':');
                if (tag_tokens.size() != 3 || tag_tokens.front().size() != 2
                    || tag_tokens[1].size() != 1) {
                    throw runtime_error("GFA header tag " + tokens[i] + " is not formatted 'XX:X:value'");
                }
                if (tag_tokens[1] == "i" && int_tags) {
                    (*int_tags)[tag_tokens.front()] = parse_int(tag_tokens.back());
                }
            }
            
        }
    }
    
    return graph;
}

}

#ifndef centrolign_three_connected_components_hpp
#define centrolign_three_connected_components_hpp

#include <vector>
#include <map>
#include <algorithm>
#include <utility>
#include <iostream>
#include <cassert>

#include "centrolign/labeled_graph.hpp"
#include "centrolign/connected_components.hpp"

namespace centrolign {

// based heavily on Adam Novak's implementation in
// https://github.com/vgteam/vg/blob/master/src/algorithms/three_edge_connected_components.cpp

// get the 3-connected components of the graph. the graph must be
// 2-connected for this algorithm to be valid
template<class Graph>
std::vector<std::vector<uint64_t>> three_edge_connected_components(const Graph& graph) {
    
    struct Node {
        /// When in the DFS were we first visited?
        size_t dfs_counter;
        /// When in the DFS were we last visited?
        /// Needed for finding replacement neighbors to implement path range
        /// absorption in part 1.3, when we're asked for a range to a neighbor
        /// that got eaten.
        size_t dfs_exit;
        /// What is our "low point" in the search. This is the earliest
        /// dfs_counter for a node that this node or any node in its DFS
        /// subtree has a back-edge to.
        size_t low_point;
        /// What is the effective degree of this node in the graph with all the
        /// absorb-eject modifications applied?
        size_t effective_degree = 0;
        /// What node has the continuation of this node's path? If equal to -1,
        /// the path ends after here. The node's path is the path from this node,
        /// into its DFS subtree, to (one of) the nodes in the subtree that has the
        /// back-edge that caused this node's low point to be so low. Basically a
        /// low point traceback.
        size_t path_tail = -1;
        /// Is this node actually on its own path?
        /// Nodes can be removed from their paths if those nodes don't matter
        /// any more (i.e. got absorbed) but their paths still need to be tails
        /// for other paths.
        bool is_on_path;
        /// Has the node been visited yet? Must be 0. TODO: Move to its own
        /// vector to make zeroing them all free-ish with page table
        /// shenanigans.
        bool visited = false;

    };
    
    // We need a DFS stack that we manage ourselves, to avoid stack-overflowing
    // as we e.g. walk along big cycles.
    struct DFSStackFrame {
        /// Track the node that this stack frame represents
        size_t current;
        /// Track all the neighbors left to visit.
        /// When we visit a neighbor we pop it off the back.
        std::vector<size_t> neighbors;
        /// When we look at the neighbors, we need to be able to tell the tree
        /// edge to the parent from further back edges to the parent. So we
        /// have a flag for whether we have seen the parent tree edge already,
        /// and the first neighbors entry that is our parent will get called
        /// the tree edge.
        bool saw_parent_tree_edge = false;
        /// Track whether we made a recursive DFS call into the last neighbor
        /// or not. If we did, we need to do some work when we come out of it
        /// and return to this frame.
        bool recursing = false;
    };

    
    // a simple graph to keep graph of merges that we've done (with 1 byte of unncessary overhead)
    LabeledGraph<uint8_t> merge_graph;
    for (uint64_t node_id = 0; node_id < graph.node_size(); ++node_id) {
        merge_graph.add_node(0);
    }
    
    
    // We need to have all the nodes pre-allocated, so node references don't
    // invalidate when we follow edges.
    std::vector<Node> nodes(graph.node_size());
    
    // We need to say how to absorb-eject along a whole path.
    //
    // We let you specify the node to absorb into; if it isn't
    // -1, it is assumed to be the first node, and
    // actually on the path, and path_start (if itself on its path) is also
    // absorbed into it. This lets you absorb into a path with something
    // prepended, without constructing the path.
    //
    // Similarly, we let you specify a past end to stop before. If this isn't
    // -1, we stop and don't absorb the specified
    // node, if we reach it. This lets us implement absorbing a range of a
    // path, as called for in the algorithm.
    //
    // If you specify a past_end, and we never reach it, but also don't have
    // just a single-node, no-edge "null" path, then something has gone wrong
    // and we've violated known truths about the algorithm.
    auto absorb_all_along_path = [&](size_t into, size_t path_start, size_t path_past_end) {
        
        // Set this to false as soon as we cross an edge
        bool path_null = true;
        
        size_t here = path_start;
        while (here != path_past_end) {
            // Until we hit the end of the path
            
            if (here == -1) {
                // We hit the end of the path and never saw path_past_end.
               
                // Only allowed if the path was actually edge-free and no merges needed to happen.
                assert(path_null);
                
                // Stop now.
                break;
            }
            
            // Find the node we are at
            auto& here_node = nodes[here];
            
            if (here_node.is_on_path) {
                // We're actually on the path.
                
                if (into == -1) {
                    // We haven't found a first node to merge into yet; it is
                    // this one.
                    into = here;
                }
                else {
                    // We already have a first node to merge into, so merge.
                    
                    // We are doing a merge! We'd better actually find the
                    // ending range bound, or something is wrong with our
                    // implementation of the algorithm.
                    path_null = false;
                    
                    // Update the effective degrees as if we merged this node
                    // with the connected into node.
                    nodes[into].effective_degree = (nodes[into].effective_degree +
                                                    here_node.effective_degree - 2);
                    
                    // Merge us into the same 3 edge connected component
                    merge_graph.add_edge(into, here);
                }
            }
            
            // Advance to the tail of the path
            here = here_node.path_tail;
        }
    };
    
    
    std::vector<DFSStackFrame> stack;
    
    // We need a way to produce unvisited nodes when we run out of nodes in a
    // connected component. This will always point to the next unvisited node
    // in order. If it points to graph.node_size(), all nodes are visited. When we
    // fisit this node, we have to scan ahead for the next unvisited node, in
    // number order.
    size_t next_unvisited = 0;
    
    // We also keep a global DFS counter, so we don't have to track parent
    // relationships when filling it in on the nodes.
    //
    // The paper starts it at 1, so we do too.
    size_t dfs_counter = 1;
    
    while (next_unvisited != graph.node_size()) {
        // We haven't visited everything yet.
        if (!nodes[0].visited) {
            // If possible start at the suggested root
            stack.emplace_back();
            stack.back().current = 0;
        } else {
            // Stack up the next unvisited node.
            stack.emplace_back();
            stack.back().current = next_unvisited;
        }
        
        while (!stack.empty()) {
            // While there's still nodes on the DFS stack from the last component we broke into
            // Grab the stack frame.
            // Note that this reference will be invalidated if we add stuff to the stack!
            auto& frame = stack.back();
            // And the current node
            auto& node = nodes[frame.current];
            
            if (!node.visited) {
                // This is the first time we are in this stack frame. We need
                // to do the initial visit of the node and set up the frame
                // with the list of edges to do.
                node.visited = true;
                
                if (frame.current == next_unvisited) {
                    // We need to find the next unvisited node, if any, since
                    // we just visited what it used to be.
                    do {
                        next_unvisited++;
                    }
                    while (next_unvisited != graph.node_size() && nodes[next_unvisited].visited);
                }
                
                node.dfs_counter = dfs_counter;
                dfs_counter++;
                node.low_point = node.dfs_counter;
                // Make sure the node's path is just itself
                node.path_tail = -1;
                node.is_on_path = true;
                
                // Stack up all the edges to follow.
                for (bool left : {true, false}) {
                    for (auto node_id : (left ? graph.previous(frame.current) : graph.next(frame.current))) {
                        frame.neighbors.push_back(node_id);
                    }
                }
                
                // Now we're in a state where we can process edges.
                // So kick back to the work loop as if we just processed an edge.
                continue;
            }
            else {
                // We have (possibly 0) edges left to do for this node.
                if (!frame.neighbors.empty()) {
                    
                    // We have an edge to do!
                    // Look up the neighboring node.
                    size_t neighbor_number = frame.neighbors.back();
                    auto& neighbor = nodes[neighbor_number];
                    
                    if (!frame.recursing) {
                        // This is the first time we are thinking about this neighbor.
                        
                        // Increment degree of the node we're coming from
                        node.effective_degree++;
                        
                        if (!neighbor.visited) {
                            // We need to recurse on this neighbor.
                            
                            // So remember we are recursing.
                            frame.recursing = true;
                            // And set up the recursive frame.
                            stack.emplace_back();
                            stack.back().current = neighbor_number;
                            // Kick back to the work loop; we will see the
                            // unvisited node on top of the stack and do its
                            // visit and add its edges to its to do list.
                        } else {
                            // No need to recurse.This is either a back-edge or the back side of the tree edge to the parent.
                            
                            if (stack.size() > 1 && neighbor_number == stack[stack.size() - 2].current && !frame.saw_parent_tree_edge) {
                                // This is the edge we took to get here (tree edge)
                                
                                // For tree edges, since they aren't either kind of back edge, neither 1.2 nor 1.3 fires.
                                // But the next edge to the parent will be a back edge.
                                frame.saw_parent_tree_edge = true;
                            } else if (neighbor.dfs_counter < node.dfs_counter) {
                                // The edge to the neighbor is an outgoing
                                // back-edge (i.e. the neighbor was visited
                                // first). Paper step 1.2.
                                
                                if (neighbor.dfs_counter < node.low_point) {
                                    // The neighbor is below our low point.
                                    
                                    // Absorb along our whole path.
                                    absorb_all_along_path(-1, frame.current, -1);
                                    
                                    // Adopt the neighbor's DFS counter as our
                                    // new, lower low point.
                                    node.low_point = neighbor.dfs_counter;
                                    
                                    // Our path is now just us.
                                    node.is_on_path = true;
                                    node.path_tail = -1;
                                    
                                }
                            }
                            else if (node.dfs_counter < neighbor.dfs_counter) {
                                // The edge to the neighbor is an incoming
                                // back-edge (i.e. we were visited first, but
                                // we recursed into something that got us to
                                // this neighbor already). Paper step 1.3.
                                
                                // Drop our effective degree by 2 (I think
                                // we're closing a cycle or something?)
                                node.effective_degree -= 2;
                                
                                // Now, the algorithm says to absorb
                                // "P_w[w..u]", a notation that it does not
                                // rigorously define. w is here, and u is the
                                // neighbor. The neighbor is not necessarily
                                // actually *on* our path at this point, not
                                // least of which because the neighbor may have
                                // already been eaten and merged into another
                                // node, which in theory adopted the back edge
                                // we are looking at. In practice we don't have
                                // the data structure to find that node. So
                                // here's the part where we have to do
                                // something clever to "allow certain paths
                                // that we track to traverse the stolen edges".
                                
                                // What we have to do is find the node that
                                // *is* along our path that either is or ate
                                // the neighbor. We don't track the union-find
                                // logic we would need to answer that question,
                                // but both 2007 algorithm implementations I've
                                // seen deal with this by tracking DFS counter
                                // intervals/subtree sizes, and deciding that
                                // the last thin on our path visited no later
                                // than the neighbor, and exited no earlier
                                // than the neighbor (i.e. the last ancestor of
                                // the neighbor on our path) should be our
                                // replacement neighbor.
                                
                                // This makes sense because if the neighbor
                                // merged into anything, it's an ancestor of
                                // the neighbor. So we go looking for it.
                                
                                // TODO: let absorb_all_along_path do this instead?
                                
                                // Start out with ourselves as the replacement neighbor ancestor.
                                size_t replacement_neighbor_number = frame.current;
                                // Consider the next candidate
                                size_t candidate = nodes[replacement_neighbor_number].path_tail;
                                while (candidate != -1 &&
                                       nodes[candidate].dfs_counter <= neighbor.dfs_counter &&
                                       nodes[candidate].dfs_exit >= neighbor.dfs_exit) {
                                    
                                    // This candidate is a lower ancestor of the neighbor, so adopt it.
                                    replacement_neighbor_number = candidate;
                                    candidate = nodes[replacement_neighbor_number].path_tail;
                                }
                                
                                auto& replacement_neighbor = nodes[replacement_neighbor_number];
                                
                                // Absorb along our path from ourselves to the
                                // replacement neighbor, inclusive.
                                // Ignores trivial paths.
                                absorb_all_along_path(-1, frame.current, replacement_neighbor.path_tail);
                                
                                // We also have to (or at least can) adopt the
                                // path of the replacement neighbor as our own
                                // path now. That's basically the rest of the
                                // path that we didn't merge.
                                // This isn't mentioned in the paper either,
                                // but I've seen the official implementation do
                                // it, and if we don't do it our path is going
                                // to go through a bunch of stuff we already
                                // merged, and waste time when we merge again.
                                
                                // If we ever merge us down our path again,
                                // continue with the part we didn't already
                                // eat.
                                node.path_tail = replacement_neighbor.path_tail;
                            }
                            else {
                                // The other possibility is the neighbor is just
                                // us. Officially self loops aren't allowed, so
                                // we censor the edge.
                                
                                node.effective_degree--;
                            }
                            
                            // Clean up the neighbor from the to do list; we
                            // finished it without recursing.
                            frame.neighbors.pop_back();
                            
                            // Kick back to the work loop to do the next
                            // neighbor, if any.
                        }
                    }
                    else {
                        // We have returned from a recursive call on this neighbor.
                        
                        // Support bridge edges: detect if we are returning
                        // across a bridge edge and censor it. Norouzi and Tsin
                        // 2014 as written in the paper assumes no bridge
                        // edges, and what we're about to do relies on all
                        // neighbors connecting back somewhere.
                        if (neighbor.low_point == neighbor.dfs_counter) {
                            // It has no back-edges out of its own subtree, so it must be across a bridge.
                            
                            // Hide the edge we just took from degree calculations.
                            neighbor.effective_degree--;
                            node.effective_degree--;
                            
                            // Don't do anything else with the edge
                        }
                        else {
                            // Wasn't a bridge edge, so we care about more than just traversing that part of the graph.
                            
                            // Do steps 1.1.1 and 1.1.2 of the algorithm as described in the paper.
                            if (neighbor.effective_degree == 2) {
                                // This neighbor gets absorbed and possibly ejected.
                                
                                // Take it off of its own path.
                                neighbor.is_on_path = false;
                            }
                            
                            // Because we hid the bridge edges, degree 1 nodes should never happen
                            assert(neighbor.effective_degree != 1);
                            
                            if (node.low_point <= neighbor.low_point) {
                                
                                // Absorb all along the path starting with here and
                                // continuing with this neighbor's path, to the
                                // end.
                                absorb_all_along_path(frame.current, neighbor_number, -1);
                            }
                            else {
                                
                                // Lower our low point to that of the neighbor
                                node.low_point = neighbor.low_point;
                                
                                // Absorb all along our own path
                                absorb_all_along_path(-1, frame.current, -1);
                                // Adjust our path to be us and then our neighbor's path
                                node.is_on_path = true;
                                node.path_tail = neighbor_number;
                            }
                        }
                        
                        // Say we aren't coming back from a recursive call
                        // anymore.
                        frame.recursing = false;
                        
                        // Clean up the neighbor,
                        frame.neighbors.pop_back();
                        
                        // Kick back to the work loop to do the next neighbor,
                        // if any.
                    }
                    
                }
                else {
                    // All the neighbors left to do for this node are done.
                    
                    // This node is done.
                    
                    // Remember when we exited it
                    node.dfs_exit = dfs_counter;
                    
                    // Clean up the stack frame.
                    stack.pop_back();
                }
            }
        }
    }
    
    return connected_components(merge_graph);
}



}

#endif /* centrolign_three_connected_components_hpp */

################################################################################
#                                                                              #
# - GranMapache: GRAphs-and-Networks MAPping Applications                      #
#   with Cython and HEuristics                                                 #
#                                                                              #
# - Module: partial_maps                                                       #
#                                                                              #
# - Description: analysis of properties of partial maps, like maximum          #
#   connected extensions, overlaps, consistency, and others.                   #
#                                                                              #
################################################################################



# dependencies #################################################################



# already in python ------------------------------------------------------------
from sys import getrecursionlimit, setrecursionlimit



# not in python ----------------------------------------------------------------
import networkx as nx



# cython specifics -------------------------------------------------------------
import cython
from libcpp cimport bool as cpp_bool
from libcpp.map cimport map as cpp_map
from libcpp.pair cimport pair as cpp_pair
from libcpp.vector cimport vector as cpp_vector
cdef extern from "<algorithm>" namespace "std":
    # find element in vector
    Iter find[Iter, Const](Iter first, Iter last, Const value)



# custom dependencies ----------------------------------------------------------
from .integerization import encode_graphs, decode_graphs, encode_match, decode_match



# algorithms ###################################################################



# functions - maximum connected extensions - wrapper ###########################



# function: callable wrapper for the maximum connected extensions --------------
def maximum_connected_extensions(G = nx.Graph(),       # can also receive a DiGraph
                                 H = nx.Graph(),       # can also receive a DiGraph
                                 input_anchor = []):   # should be non-empty list
    # description
    """
    > description: receives two graphs G and H, and a match between them (here called
    anchor), and uses a VF2-like approach to obtain the maximum extensions of the anchor
    producing connected common subgraphs (not necessarily maximum themselves). The anchor
    alone also produces a subgraph, which may not be an induced common subgraph, but
    the subgraph produced by any extension after removing the achor is always induced.

    > input:
    * G - first networkx (di)graph being matched.
    * H - second networkx (di)graph being matched.
    * input_anchor - inyective map as a non-empty list of 2-tuples (x, y) of nodes x
    from G and y from H. An exception is raised if the anchor is empty.

    > output:
    * extensions - list of injective maps each as a list of 2-tuples (x, y) of nodes x
    from G and y from H representing the maximum connected extensions of the anchor (each
    extension contains the anchor as a sublist).
    * good_anchor - boolean value indicating if the extensions cover all nodes of G and
    of H, i.e., if they are bijections between G and H. If so, the anchor is what we
    have refered to as a "good partial atom map", and equivalenteĺy the match obtained
    when removing the anchhor from any extension is a graph-isomorphism between the
    "remainder" graphs it induces from G and H.

    > calls:
    * .integerization.encode_graphs
    * .integerization.decode_graphs
    * .integerization.encode_match
    * .integerization.decode_match
    * undirected_maximum_connected_extensions
    * directed_maximum_connected_extensions
    """
    # exception handling and input correctness
    test_list = [0, 0]
    test_tuple = (0, 0)
    test_undir = nx.Graph()
    test_dir = nx.DiGraph()
    if(type(G) not in [type(test_undir), type(test_dir)]):
        raise(ValueError("first argument must be a networkx graph or digraph."))
    if(type(H) not in [type(test_undir), type(test_dir)]):
        raise(ValueError("second argument must be a networkx graph or digraph."))
    if((nx.is_directed(G)) and (not nx.is_directed(H))):
        raise(ValueError("input graphs must be both directed or both undirected."))
    if((not nx.is_directed(G)) and (nx.is_directed(H))):
        raise(ValueError("input graphs must be both directed or both undirected."))
    if(not type(input_anchor) in [type(test_list)]):
        raise(ValueError("third argument must be a non-empty list of 2-tuples."))
    if(len(input_anchor) == 0):
        raise(ValueError("third argument must be a non-empty list of 2-tuples."))
    for test_entry in input_anchor:
        if(not type(test_entry) in [type(test_tuple)]):
            raise(ValueError("all elements in input list must be tuples."))
        if(not len(test_entry) == 2):
            raise(ValueError("all tuples in input list must be of lenght 2."))
        if(test_entry[0] not in list(G.nodes())):
            raise(ValueError("the input list is matching a vertex not present in the first graph."))
        if(test_entry[1] not in list(H.nodes())):
            raise(ValueError("the input list is matching a vertex not present in the second graph."))
    if(not len(list(set([x for (x, y) in input_anchor]))) == len(input_anchor)):
        raise(ValueError("the input list must be an injective map and without repeated elements."))
    if(not len(list(set([y for (x, y) in input_anchor]))) == len(input_anchor)):
        raise(ValueError("the input list must be an injective map and without repeated elements."))
    # output holders
    cdef list extensions = []
    good_anchor = False
    # local variables (cython)
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int counter = 0
    cdef int current_limit = 0
    cdef int required_limit = 0
    cdef float scalation_value = 0
    cdef long unsigned int expected_order = 0
    cdef cpp_vector[cpp_pair[int, int]] encoded_anchor
    cdef cpp_vector[cpp_vector[cpp_pair[int, int]]] encoded_extensions
    cdef cpp_map[int, int] total_order
    cdef cpp_map[int, int] nodes_G
    cdef cpp_map[int, int] nodes_H
    cdef cpp_map[int, cpp_vector[int]] neigh_G        # only used if undirected
    cdef cpp_map[int, cpp_vector[int]] neigh_H        # only used if undirected
    cdef cpp_map[int, cpp_vector[int]] in_neigh_G     # only used if directed
    cdef cpp_map[int, cpp_vector[int]] in_neigh_H     # only used if directed
    cdef cpp_map[int, cpp_vector[int]] out_neigh_G    # only used if directed
    cdef cpp_map[int, cpp_vector[int]] out_neigh_H    # only used if directed
    cdef cpp_map[cpp_pair[int, int], int] edges_G
    cdef cpp_map[cpp_pair[int, int], int] edges_H
    # local variables (python)
    cdef list encoded_graphs = []
    cdef list inside_anchor_H = []
    cdef list outside_anchor_H = []
    cdef dict info = dict()
    cdef dict encoded_node_names = dict()
    cdef dict encoded_node_label = dict()
    cdef dict encoded_edge_label = dict()
    # encode graphs
    encoded_graphs, encoded_node_names, encoded_node_label, encoded_edge_label = encode_graphs([G, H])
    # encode match
    encoded_anchor = encode_match(input_anchor, encoded_node_names)
    # prepare nodes
    nodes_G = {node:info["GMNL"] for (node, info) in encoded_graphs[0].nodes(data = True)}
    nodes_H = {node:info["GMNL"] for (node, info) in encoded_graphs[1].nodes(data = True)}
    # prepare edges
    if(nx.is_directed(G)):
        edges_G = {(node1, node2):info["GMEL"] for (node1, node2, info) in encoded_graphs[0].edges(data = True)}
        edges_H = {(node1, node2):info["GMEL"] for (node1, node2, info) in encoded_graphs[1].edges(data = True)}
    else:
        edges_G = {tuple(sorted([node1, node2])):info["GMEL"] for (node1, node2, info) in encoded_graphs[0].edges(data = True)}
        edges_H = {tuple(sorted([node1, node2])):info["GMEL"] for (node1, node2, info) in encoded_graphs[1].edges(data = True)}
    # prepare neighbors
    if(nx.is_directed(G)):
        in_neigh_G = {node:list(encoded_graphs[0].predecessors(node)) for node in list(encoded_graphs[0].nodes())}
        in_neigh_H = {node:list(encoded_graphs[1].predecessors(node)) for node in list(encoded_graphs[1].nodes())}
        out_neigh_G = {node:list(encoded_graphs[0].neighbors(node)) for node in list(encoded_graphs[0].nodes())}
        out_neigh_H = {node:list(encoded_graphs[1].neighbors(node)) for node in list(encoded_graphs[1].nodes())}
    else:
        neigh_G = {node:list(encoded_graphs[0].neighbors(node)) for node in list(encoded_graphs[0].nodes())}
        neigh_H = {node:list(encoded_graphs[1].neighbors(node)) for node in list(encoded_graphs[1].nodes())}
    # get total order for VF2-like analysis
    inside_anchor_H = [node2 for (node1, node2) in encoded_anchor]
    outside_anchor_H = [node for node in list(encoded_graphs[1].nodes()) if(node not in inside_anchor_H)]
    for node in inside_anchor_H:
        counter = counter + 1
        total_order[node] = counter
    for node in outside_anchor_H:
        counter = counter + 1
        total_order[node] = counter
    # get expected order
    expected_order = min([len(nodes_G), len(nodes_H)])
    # set recursion limit
    scalation_value = 1.5
    required_limit = max([len(nodes_G), len(nodes_H)])
    current_limit = getrecursionlimit()
    if(current_limit < (scalation_value * required_limit)):
        setrecursionlimit(int(scalation_value * required_limit))
    # get maximum extensions
    if(nx.is_directed(G)):
        directed_maximum_connected_extensions(nodes_G,
                                              nodes_H,
                                              in_neigh_G,
                                              in_neigh_H,
                                              out_neigh_G,
                                              out_neigh_H,
                                              edges_G,
                                              edges_H,
                                              total_order,
                                              expected_order,
                                              encoded_anchor,
                                              encoded_extensions)
    else:
        undirected_maximum_connected_extensions(nodes_G,
                                                nodes_H,
                                                neigh_G,
                                                neigh_H,
                                                edges_G,
                                                edges_H,
                                                total_order,
                                                expected_order,
                                                encoded_anchor,
                                                encoded_extensions)
    # decode maximum extensions
    for each_extension in encoded_extensions:
        extensions.append(decode_match(list(each_extension), encoded_node_names))
    # check if the anchor was a good partial map
    if((len(extensions[0]) == len(nodes_G)) and (len(extensions[0]) == len(nodes_H))):
        good_anchor = True
    # end of function
    return(extensions, good_anchor)



# functions - maximum connected extensions - undirected ########################



# function: core routine of VF2-like undirected approach -----------------------
cdef void undirected_maximum_connected_extensions(cpp_map[int, int] & nodes_G,
                                                  cpp_map[int, int] & nodes_H,
                                                  cpp_map[int, cpp_vector[int]] & neigh_G,
                                                  cpp_map[int, cpp_vector[int]] & neigh_H,
                                                  cpp_map[cpp_pair[int, int], int] & edges_G,
                                                  cpp_map[cpp_pair[int, int], int] & edges_H,
                                                  cpp_map[int, int] & total_order,
                                                  long unsigned int & expected_order,
                                                  cpp_vector[cpp_pair[int, int]] current_match,
                                                  cpp_vector[cpp_vector[cpp_pair[int, int]]] & all_matches) noexcept:
    # local variables (cython)
    cdef long unsigned int new_score = 0
    cdef long unsigned int old_score = 0
    cdef cpp_bool semantic_feasibility = False
    cdef cpp_bool syntactic_feasibility = False
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[int] current_match_G
    cdef cpp_vector[int] current_match_H
    cdef cpp_vector[cpp_pair[int, int]] new_match
    cdef cpp_vector[cpp_pair[int, int]] candidates
    cdef cpp_map[int, int] forward_match
    # test initial match and improve if possible
    if(all_matches.empty()):
        if(not current_match.empty()):
            all_matches.push_back(current_match)
    else:
        # test improvement in matches
        new_score = current_match.size()
        old_score = all_matches[0].size()
        # save match if it has the same score
        if(new_score == old_score):
            all_matches.push_back(current_match)
        # overwrite everything with new match if it improves score
        if(new_score > old_score):
            all_matches.clear()
            all_matches.push_back(current_match)
    # if not optimal yet then obtain available pairs
    if(current_match.size() < expected_order):
        # generate auxiliary structures
        for each_pair in current_match:
            current_match_G.push_back(each_pair.first)
            current_match_H.push_back(each_pair.second)
            forward_match[each_pair.first] = each_pair.second
        # get candidate pairs
        candidates = undirected_candidates(current_match_G,
                                           current_match_H,
                                           neigh_G,
                                           neigh_H,
                                           total_order)
        # evaluate candidates
        for each_pair in candidates:
            # evaluate sintactic feasibility
            syntactic_feasibility = undirected_syntactic_feasibility(each_pair.first,
                                                                     each_pair.second,
                                                                     current_match_G,
                                                                     current_match_H,
                                                                     forward_match,
                                                                     neigh_G,
                                                                     neigh_H)
            if(syntactic_feasibility):
                # evaluate semantic feasibility
                semantic_feasibility = undirected_semantic_feasibility(each_pair.first,
                                                                       each_pair.second,
                                                                       current_match_G,
                                                                       forward_match,
                                                                       nodes_G,
                                                                       nodes_H,
                                                                       neigh_G,
                                                                       edges_G,
                                                                       edges_H)
                if(semantic_feasibility):
                    # build new match
                    new_match = current_match
                    new_match.push_back(each_pair)
                    # extend match
                    undirected_maximum_connected_extensions(nodes_G,
                                                            nodes_H,
                                                            neigh_G,
                                                            neigh_H,
                                                            edges_G,
                                                            edges_H,
                                                            total_order,
                                                            expected_order,
                                                            new_match,
                                                            all_matches)
    # end of function



# function: get candidate pairs for undirected extension search ----------------
cdef cpp_vector[cpp_pair[int, int]] undirected_candidates(cpp_vector[int] & current_match_G,
                                                          cpp_vector[int] & current_match_H,
                                                          cpp_map[int, cpp_vector[int]] & neigh_G,
                                                          cpp_map[int, cpp_vector[int]] & neigh_H,
                                                          cpp_map[int, int] & total_order) noexcept:
    # local variables
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int reference_maximum = 0
    cdef cpp_pair[int, int] temp_pair
    cdef cpp_vector[int] valid_G
    cdef cpp_vector[int] valid_H
    cdef cpp_vector[cpp_pair[int, int]] candidate_pairs
    # get maximum value of total order in match
    for node in current_match_H:
        if(total_order[node] > reference_maximum):
            reference_maximum = total_order[node]
    # get valid sets
    for node1 in current_match_G:
        for node2 in neigh_G[node1]:
            # if not already valid in G
            if(find(valid_G.begin(), valid_G.end(), node2) == valid_G.end()):
                # if not yet in match
                if(find(current_match_G.begin(), current_match_G.end(), node2) == current_match_G.end()):
                    valid_G.push_back(node2)
    for node1 in current_match_H:
        for node2 in neigh_H[node1]:
            # if not already valid in H
            if(find(valid_H.begin(), valid_H.end(), node2) == valid_H.end()):
                # if not yet in match
                if(find(current_match_H.begin(), current_match_H.end(), node2) == current_match_H.end()):
                    # if total order greater than in match
                    if(total_order[node2] > reference_maximum):
                        valid_H.push_back(node2)
    # get candidates
    if((not valid_G.empty()) and (not valid_H.empty())):
        for node1 in valid_G:
            for node2 in valid_H:
                temp_pair.first = node1
                temp_pair.second = node2
                candidate_pairs.push_back(temp_pair)
    # end of function
    return(candidate_pairs)



# function: evaluate syntactic feasability for undirected extension ------------
cdef cpp_bool undirected_syntactic_feasibility(int & node1,
                                               int & node2,
                                               cpp_vector[int] & current_match_G,
                                               cpp_vector[int] & current_match_H,
                                               cpp_map[int, int] & forward_match,
                                               cpp_map[int, cpp_vector[int]] & neigh_G,
                                               cpp_map[int, cpp_vector[int]] & neigh_H) noexcept:
    # local variables
    cdef int node = 0
    cdef int mapped = 0
    cdef cpp_vector[int] neighbors_match_G
    cdef cpp_vector[int] neighbors_match_H
    # loop-consistency-test
    if(find(neigh_G[node1].begin(), neigh_G[node1].end(), node1) != neigh_G[node1].end()):
        if(find(neigh_H[node2].begin(), neigh_H[node2].end(), node2) == neigh_H[node2].end()):
            # node1 has a loop in G but node2 has no loop in H
            return(False)
    if(find(neigh_G[node1].begin(), neigh_G[node1].end(), node1) == neigh_G[node1].end()):
        if(find(neigh_H[node2].begin(), neigh_H[node2].end(), node2) != neigh_H[node2].end()):
            # node1 has no loop in G but node2 has a loop in H
            return(False)
    # look ahead 0: consistency of neighbors in match
    for node in neigh_G[node1]:
        if(find(current_match_G.begin(), current_match_G.end(), node) != current_match_G.end()):
            mapped = forward_match[node]
            neighbors_match_G.push_back(mapped)
    for node in neigh_H[node2]:
        if(find(current_match_H.begin(), current_match_H.end(), node) != current_match_H.end()):
            neighbors_match_H.push_back(node)
    if(neighbors_match_G.size() != neighbors_match_H.size()):
        # one node has more neighbors in the match than the other
        return(False)
    else:
        for mapped in neighbors_match_G:
            if(find(neighbors_match_H.begin(), neighbors_match_H.end(), mapped) == neighbors_match_H.end()):
                # the neighbors dont respect the match
                return(False)
    # end of function
    return(True)



# function: evaluate semantic feasability for undirected extension -------------
cdef cpp_bool undirected_semantic_feasibility(int & node1,
                                              int & node2,
                                              cpp_vector[int] & current_match_G,
                                              cpp_map[int, int] & forward_match,
                                              cpp_map[int, int] & nodes_G,
                                              cpp_map[int, int] & nodes_H,
                                              cpp_map[int, cpp_vector[int]] & neigh_G,
                                              cpp_map[cpp_pair[int, int], int] & edges_G,
                                              cpp_map[cpp_pair[int, int], int] & edges_H) noexcept:
    # local variables
    cdef int node = 0
    cdef int mapped = 0
    cdef cpp_pair[int, int] labeled_edge_G
    cdef cpp_pair[int, int] labeled_edge_H
    cdef cpp_vector[int] neighbors_match_G
    # compare vertex-labels
    if(nodes_G[node1] != nodes_H[node2]):
        return(False)
    # compare loop-labels
    if(find(neigh_G[node1].begin(), neigh_G[node1].end(), node1) != neigh_G[node1].end()):
        labeled_edge_G.first = node1
        labeled_edge_G.second = node1
        labeled_edge_H.first = node2
        labeled_edge_H.second = node2
        if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
            return(False)
    # compare non-loop edge-labels
    for node in neigh_G[node1]:
        if(find(current_match_G.begin(), current_match_G.end(), node) != current_match_G.end()):
            neighbors_match_G.push_back(node)
    for node in neighbors_match_G:
        # edge in G with only one end in match
        if(node1 < node):
            labeled_edge_G.first = node1
            labeled_edge_G.second = node
        else:
            labeled_edge_G.first = node
            labeled_edge_G.second = node1
        # edge in H with only one end in match
        mapped = forward_match[node]
        if(node2 < mapped):
            labeled_edge_H.first = node2
            labeled_edge_H.second = mapped
        else:
            labeled_edge_H.first = mapped
            labeled_edge_H.second = node2
        if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
            return(False)
    # end of function
    return(True)



# functions - maximum connected extensions - directed ##########################



# function: core routine of VF2-like directed approach -------------------------
cdef void directed_maximum_connected_extensions(cpp_map[int, int] & nodes_G,
                                                cpp_map[int, int] & nodes_H,
                                                cpp_map[int, cpp_vector[int]] & in_neigh_G,
                                                cpp_map[int, cpp_vector[int]] & in_neigh_H,
                                                cpp_map[int, cpp_vector[int]] & out_neigh_G,
                                                cpp_map[int, cpp_vector[int]] & out_neigh_H,
                                                cpp_map[cpp_pair[int, int], int] & edges_G,
                                                cpp_map[cpp_pair[int, int], int] & edges_H,
                                                cpp_map[int, int] & total_order,
                                                long unsigned int & expected_order,
                                                cpp_vector[cpp_pair[int, int]] current_match,
                                                cpp_vector[cpp_vector[cpp_pair[int, int]]] & all_matches) noexcept:
    # local variables
    cdef long unsigned int new_score = 0
    cdef long unsigned int old_score = 0
    cdef cpp_bool semantic_feasibility = False
    cdef cpp_bool syntactic_feasibility = False
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[int] current_match_G
    cdef cpp_vector[int] current_match_H
    cdef cpp_vector[cpp_pair[int, int]] new_match
    cdef cpp_vector[cpp_pair[int, int]] candidates
    cdef cpp_map[int, int] forward_match
    # test initial match and consecutive matches
    if(all_matches.empty()):
        if(not current_match.empty()):
            all_matches.push_back(current_match)
    else:
        # test improvement in matches
        new_score = current_match.size()
        old_score = all_matches[0].size()
        # save match if it has the same score
        if(new_score == old_score):
            all_matches.push_back(current_match)
        # overwrite with new match if it improves score
        if(new_score > old_score):
            all_matches.clear()
            all_matches.push_back(current_match)
    # if not optimal yet then obtain available pairs
    if(current_match.size() < expected_order):
        # generate auxiliary structures
        for each_pair in current_match:
            current_match_G.push_back(each_pair.first)
            current_match_H.push_back(each_pair.second)
            forward_match[each_pair.first] = each_pair.second
        # get candidate pairs
        candidates = directed_candidates(current_match_G,
                                         current_match_H,
                                         in_neigh_G,
                                         in_neigh_H,
                                         out_neigh_G,
                                         out_neigh_H,
                                         total_order)
        # evaluate candidates
        for each_pair in candidates:
            # evaluate sintactic feasibility
            syntactic_feasibility = directed_syntactic_feasibility(each_pair.first,
                                                                   each_pair.second,
                                                                   current_match_G,
                                                                   current_match_H,
                                                                   forward_match,
                                                                   in_neigh_G,
                                                                   in_neigh_H,
                                                                   out_neigh_G,
                                                                   out_neigh_H)
            if(syntactic_feasibility):
                # evaluate semantic feasibility
                semantic_feasibility = directed_semantic_feasibility(each_pair.first,
                                                                     each_pair.second,
                                                                     current_match_G,
                                                                     forward_match,
                                                                     nodes_G,
                                                                     nodes_H,
                                                                     in_neigh_G,
                                                                     out_neigh_G,
                                                                     edges_G,
                                                                     edges_H)
                if(semantic_feasibility):
                    # build new match
                    new_match = current_match
                    new_match.push_back(each_pair)
                    # extend match
                    directed_maximum_connected_extensions(nodes_G,
                                                          nodes_H,
                                                          in_neigh_G,
                                                          in_neigh_H,
                                                          out_neigh_G,
                                                          out_neigh_H,
                                                          edges_G,
                                                          edges_H,
                                                          total_order,
                                                          expected_order,
                                                          new_match,
                                                          all_matches)
    # end of function



# function: get candidate pairs for directed extension search ------------------
cdef cpp_vector[cpp_pair[int, int]] directed_candidates(cpp_vector[int] & current_match_G,
                                                        cpp_vector[int] & current_match_H,
                                                        cpp_map[int, cpp_vector[int]] & in_neigh_G,
                                                        cpp_map[int, cpp_vector[int]] & in_neigh_H,
                                                        cpp_map[int, cpp_vector[int]] & out_neigh_G,
                                                        cpp_map[int, cpp_vector[int]] & out_neigh_H,
                                                        cpp_map[int, int] & total_order) noexcept:
    # local variables
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int reference_maximum = 0
    cdef cpp_pair[int, int] temp_pair
    cdef cpp_vector[int] valid_G
    cdef cpp_vector[int] valid_H
    cdef cpp_vector[int] temp_vector
    cdef cpp_vector[int] all_neighbors_G
    cdef cpp_vector[int] all_neighbors_H
    cdef cpp_vector[cpp_pair[int, int]] candidate_pairs
    # get maximum value of total order in match
    for node in current_match_H:
        if(total_order[node] > reference_maximum):
            reference_maximum = total_order[node]
    # get valid sets
    for node1 in current_match_G:
        all_neighbors_G = in_neigh_G[node1]
        temp_vector = out_neigh_G[node1]
        all_neighbors_G.insert(all_neighbors_G.end(), temp_vector.begin(), temp_vector.end())
        for node2 in all_neighbors_G:
            # if not already valid in G
            if(find(valid_G.begin(), valid_G.end(), node2) == valid_G.end()):
                # if not yet in match
                if(find(current_match_G.begin(), current_match_G.end(), node2) == current_match_G.end()):
                    valid_G.push_back(node2)
    for node1 in current_match_H:
        all_neighbors_H = in_neigh_H[node1]
        temp_vector = out_neigh_H[node1]
        all_neighbors_H.insert(all_neighbors_H.end(), temp_vector.begin(), temp_vector.end())
        for node2 in all_neighbors_H:
            # if not already valid in H
            if(find(valid_H.begin(), valid_H.end(), node2) == valid_H.end()):
                # if not yet in match
                if(find(current_match_H.begin(), current_match_H.end(), node2) == current_match_H.end()):
                    # if total order greater than in match
                    if(total_order[node2] > reference_maximum):
                        valid_H.push_back(node2)
    # get candidates
    if((not valid_G.empty()) and (not valid_H.empty())):
        for node1 in valid_G:
            for node2 in valid_H:
                temp_pair.first = node1
                temp_pair.second = node2
                candidate_pairs.push_back(temp_pair)
    # end of function
    return(candidate_pairs)



# function: evaluate syntactic feasability for directed extension --------------
cdef cpp_bool directed_syntactic_feasibility(int & node1,
                                             int & node2,
                                             cpp_vector[int] & current_match_G,
                                             cpp_vector[int] & current_match_H,
                                             cpp_map[int, int] & forward_match,
                                             cpp_map[int, cpp_vector[int]] & in_neigh_G,
                                             cpp_map[int, cpp_vector[int]] & in_neigh_H,
                                             cpp_map[int, cpp_vector[int]] & out_neigh_G,
                                             cpp_map[int, cpp_vector[int]] & out_neigh_H) noexcept:
    # local variables
    cdef int node = 0
    cdef int mapped = 0
    cdef cpp_vector[int] in_neighbors_match_G
    cdef cpp_vector[int] in_neighbors_match_H
    cdef cpp_vector[int] out_neighbors_match_G
    cdef cpp_vector[int] out_neighbors_match_H
    # loop-consistency-test
    if(find(in_neigh_G[node1].begin(), in_neigh_G[node1].end(), node1) != in_neigh_G[node1].end()):
        if(find(in_neigh_H[node2].begin(), in_neigh_H[node2].end(), node2) == in_neigh_H[node2].end()):
            # node1 has a loop in G but node2 has no loop in H
            return(False)
    if(find(in_neigh_G[node1].begin(), in_neigh_G[node1].end(), node1) == in_neigh_G[node1].end()):
        if(find(in_neigh_H[node2].begin(), in_neigh_H[node2].end(), node2) != in_neigh_H[node2].end()):
            # node1 has no loop in G but node2 has a loop in H
            return(False)
    # look ahead 0: consistency of in-neighbors in match
    for node in in_neigh_G[node1]:
        if(find(current_match_G.begin(), current_match_G.end(), node) != current_match_G.end()):
            mapped = forward_match[node]
            in_neighbors_match_G.push_back(mapped)
    for node in in_neigh_H[node2]:
        if(find(current_match_H.begin(), current_match_H.end(), node) != current_match_H.end()):
            in_neighbors_match_H.push_back(node)
    if(in_neighbors_match_G.size() != in_neighbors_match_H.size()):
        # one node has more in-neighbors in the match than the other
        return(False)
    else:
        for mapped in in_neighbors_match_G:
            if(find(in_neighbors_match_H.begin(), in_neighbors_match_H.end(), mapped) == in_neighbors_match_H.end()):
                # the in-neighbors dont respect the match
                return(False)
    # look ahead 0: consistency of out-neighbors in match
    for node in out_neigh_G[node1]:
        if(find(current_match_G.begin(), current_match_G.end(), node) != current_match_G.end()):
            mapped = forward_match[node]
            out_neighbors_match_G.push_back(mapped)
    for node in out_neigh_H[node2]:
        if(find(current_match_H.begin(), current_match_H.end(), node) != current_match_H.end()):
            out_neighbors_match_H.push_back(node)
    if(out_neighbors_match_G.size() != out_neighbors_match_H.size()):
        # one node has more out-neighbors in the match than the other
        return(False)
    else:
        for mapped in out_neighbors_match_G:
            if(find(out_neighbors_match_H.begin(), out_neighbors_match_H.end(), mapped) == out_neighbors_match_H.end()):
                # the out-neighbors dont respect the match
                return(False)
    # end of function
    return(True)



# function: evaluate semantic feasability for directed extension ---------------
cdef cpp_bool directed_semantic_feasibility(int & node1,
                                            int & node2,
                                            cpp_vector[int] & current_match_G,
                                            cpp_map[int, int] & forward_match,
                                            cpp_map[int, int] & nodes_G,
                                            cpp_map[int, int] & nodes_H,
                                            cpp_map[int, cpp_vector[int]] & in_neigh_G,
                                            cpp_map[int, cpp_vector[int]] & out_neigh_G,
                                            cpp_map[cpp_pair[int, int], int] & edges_G,
                                            cpp_map[cpp_pair[int, int], int] & edges_H) noexcept:
    # local variables
    cdef int node = 0
    cdef int mapped = 0
    cdef cpp_pair[int, int] labeled_edge_G
    cdef cpp_pair[int, int] labeled_edge_H
    cdef cpp_vector[int] in_neighbors_match_G
    cdef cpp_vector[int] out_neighbors_match_G
    # compare vertex-labels
    if(nodes_G[node1] != nodes_H[node2]):
        return(False)
    # compare loop-labels
    if(find(in_neigh_G[node1].begin(), in_neigh_G[node1].end(), node1) != in_neigh_G[node1].end()):
        labeled_edge_G.first = node1
        labeled_edge_G.second = node1
        labeled_edge_H.first = node2
        labeled_edge_H.second = node2
        if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
            return(False)
    # compare non-loop in-edge-labels
    for node in in_neigh_G[node1]:
        if(find(current_match_G.begin(), current_match_G.end(), node) != current_match_G.end()):
            in_neighbors_match_G.push_back(node)
    for node in in_neighbors_match_G:
        # edge in G with only one end in match
        labeled_edge_G.first = node
        labeled_edge_G.second = node1
        labeled_edge_H.first = forward_match[node]
        labeled_edge_H.second = node2
        if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
            return(False)
    # compare non-loop out-edge-labels
    for node in out_neigh_G[node1]:
        if(find(current_match_G.begin(), current_match_G.end(), node) != current_match_G.end()):
            out_neighbors_match_G.push_back(node)
    for node in out_neighbors_match_G:
        # edge in G with only one end in match
        labeled_edge_G.first = node1
        labeled_edge_G.second = node
        labeled_edge_H.first = node2
        labeled_edge_H.second = forward_match[node]
        if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
            return(False)
    # end of function
    return(True)



################################################################################
################################################################################

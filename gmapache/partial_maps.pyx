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
from copy import deepcopy
from sys import getrecursionlimit, setrecursionlimit




# not in python ----------------------------------------------------------------
import networkx as nx




# cython specifics -------------------------------------------------------------
import cython
from libcpp cimport bool as cpp_bool
from libcpp.map cimport map as cpp_map
from libcpp.pair cimport pair as cpp_pair
from libcpp.stack cimport stack as cpp_stack
from libcpp.vector cimport vector as cpp_vector
from libcpp.unordered_map cimport unordered_map as cpp_unordered_map
cdef extern from "<algorithm>" namespace "std":
    # find element in vector
    Iter find[Iter, Const](Iter first, Iter last, Const value)




# custom dependencies ----------------------------------------------------------
from .integerization import encode_graphs, decode_graphs, encode_match, decode_match




# C/C++ structs ################################################################




# struct: undirected graph ---------------------------------------------------
cdef struct partial_maps_undirected_graph:
    cpp_unordered_map[int, int] nodes
    cpp_map[cpp_pair[int, int], int] edges
    cpp_unordered_map[int, cpp_vector[int]] neighbors




# struct: directed graph -------------------------------------------------------
cdef struct partial_maps_directed_graph:
    cpp_unordered_map[int, int] nodes
    cpp_map[cpp_pair[int, int], int] edges
    cpp_unordered_map[int, cpp_vector[int]] in_neighbors
    cpp_unordered_map[int, cpp_vector[int]] out_neighbors




# algorithms ###################################################################




# functions - maximum connected extensions - wrapper ###########################




# function: callable wrapper for the maximum connected extensions --------------
def maximum_connected_extensions(nx_G = nx.Graph(),          # can be nx.DiGraph
                                 nx_H = nx.Graph(),          # can be nx.DiGraph
                                 input_anchor = [],          # should be non-empty list
                                 all_extensions = False,     # by default stops when finding one complete extension (if any)
                                 iterative_search = True):   # by default an iterative search is used, otherwise a recursive version is called
    # description
    """
    > description: receives two networkx graphs G and H, and a match between them (here
    called anchor), and uses a VF2-like approach to obtain the maximum extensions of the
    anchor producing connected common subgraphs (not necessarily maximum themselves). The
    anchor alone also produces a subgraph, which may not be an induced common subgraph,
    but the subgraph produced by any extension after removing the achor is always induced.

    > input:
    * nx_G - first networkx (di)graph being matched.
    * nx_H - second networkx (di)graph being matched.
    * input_anchor - inyective map as a non-empty list of 2-tuples (x, y) of nodes x
    from G and y from H. An exception is raised if the anchor is empty.
    * all_extensions - boolean indicating if the function should stop as soon as one complete
    extension is found (if any) - this is the default behavior- or if it should search for all
    possible (complete) extensions. NOTE: mathematically speaking, the complete extension of
    a FIXED and GOOD anchor is unique up to equivalence of bijections, that is equivalence
    of atom maps or correspondingly isomorphism of their ITS graphs. Nontheless it should be
    noted that changing the anchor may produce non-equivalent (complete) extensions. In other
    words, each call to this function can (mathematically) produce only one complete extension,
    but calls with different anchors can produce non-equivalent extensions, even if the anchors
    produce isomorphic partial ITS graphs.
    * iterative_search - boolean indicating if the iterative version of this algorithm should
    used (the default), or if a recursive version of it should be used instead.

    > output:
    * extensions - list of injective maps each as a list of 2-tuples (x, y) of nodes x
    from G and y from H representing the maximum connected extensions of the anchor (each
    extension contains the anchor as a sublist).
    * good_anchor - boolean value indicating if the extensions are complete, i.e., cover
    all nodes of G and of H, and thus if they are bijections between G and H. If so, the
    anchor is what we have refered to as a "good partial atom map", and equivalenteĺy the
    match obtained when removing the anchhor from any extension is a graph-isomorphism
    between the "remainder" graphs it induces from G and H.

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
    test_bool = False
    if(type(nx_G) not in [type(test_undir), type(test_dir)]):
        raise(ValueError("gmapache: first argument must be a networkx graph or digraph."))
    if(type(nx_H) not in [type(test_undir), type(test_dir)]):
        raise(ValueError("gmapache: second argument must be a networkx graph or digraph."))
    if((nx.is_directed(nx_G)) and (not nx.is_directed(nx_H))):
        raise(ValueError("gmapache: input graphs must be both directed or both undirected."))
    if((not nx.is_directed(nx_G)) and (nx.is_directed(nx_H))):
        raise(ValueError("gmapache: input graphs must be both directed or both undirected."))
    if(type(all_extensions) not in [type(test_bool)]):
        raise(ValueError("gmapache: fourth argument must be a boolean variable."))
    if(type(iterative_search) not in [type(test_bool)]):
        raise(ValueError("gmapache: fifth argument must be a boolean variable."))
    if(not type(input_anchor) in [type(test_list)]):
        raise(ValueError("gmapache: third argument must be a non-empty list of 2-tuples."))
    if(len(input_anchor) == 0):
        raise(ValueError("gmapache: third argument must be a non-empty list of 2-tuples."))
    for test_entry in input_anchor:
        if(not type(test_entry) in [type(test_tuple)]):
            raise(ValueError("gmapache: all elements in input list must be tuples."))
        if(not len(test_entry) == 2):
            raise(ValueError("gmapache: all tuples in input list must be of lenght 2."))
        if(test_entry[0] not in list(nx_G.nodes())):
            raise(ValueError("gmapache: the input list is matching a vertex not present in the first graph."))
        if(test_entry[1] not in list(nx_H.nodes())):
            raise(ValueError("gmapache: the input list is matching a vertex not present in the second graph."))
    if(not len(list(set([x for (x, y) in input_anchor]))) == len(input_anchor)):
        raise(ValueError("gmapache: the input list must be an injective map and without repeated elements."))
    if(not len(list(set([y for (x, y) in input_anchor]))) == len(input_anchor)):
        raise(ValueError("gmapache: the input list must be an injective map and without repeated elements."))

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
    cdef size_t expected_order = 0
    cdef cpp_pair[int, int] temp_pair
    cdef cpp_pair[int, cpp_vector[int]] each_pair
    cdef cpp_vector[int] next_level
    cdef cpp_vector[int] current_level
    cdef cpp_vector[cpp_pair[int, int]] encoded_anchor
    cdef cpp_vector[cpp_pair[int, int]] each_extension
    cdef cpp_vector[cpp_vector[cpp_pair[int, int]]] encoded_extensions
    cdef cpp_unordered_map[int, int] total_order
    cdef cpp_unordered_map[int, cpp_bool] visited
    cdef cpp_unordered_map[int, cpp_vector[int]] connectivity_neighbors
    cdef partial_maps_undirected_graph undirected_G
    cdef partial_maps_undirected_graph undirected_H
    cdef partial_maps_directed_graph directed_G
    cdef partial_maps_directed_graph directed_H

    # local variables (python)
    cdef list encoded_graphs = []
    cdef dict info = dict()
    cdef dict encoded_node_names = dict()
    cdef dict encoded_node_label = dict()
    cdef dict encoded_edge_label = dict()
    undirected_copy_H = None

    # encode graphs
    encoded_graphs, encoded_node_names, encoded_node_label, encoded_edge_label = encode_graphs([nx_G, nx_H])

    # encode match
    encoded_anchor = encode_match(input_anchor, encoded_node_names)

    # prepare nodes
    if(nx.is_directed(nx_G)):
        directed_G.nodes = {node:info["GMNL"] for (node, info) in encoded_graphs[0].nodes(data = True)}
        directed_H.nodes = {node:info["GMNL"] for (node, info) in encoded_graphs[1].nodes(data = True)}
    else:
        undirected_G.nodes = {node:info["GMNL"] for (node, info) in encoded_graphs[0].nodes(data = True)}
        undirected_H.nodes = {node:info["GMNL"] for (node, info) in encoded_graphs[1].nodes(data = True)}

    # prepare edges
    if(nx.is_directed(nx_G)):
        directed_G.edges = {(node1, node2):info["GMEL"] for (node1, node2, info) in encoded_graphs[0].edges(data = True)}
        directed_H.edges = {(node1, node2):info["GMEL"] for (node1, node2, info) in encoded_graphs[1].edges(data = True)}
    else:
        # prepare edges of undirected G
        for (node1, node2, info) in encoded_graphs[0].edges(data = True):
            if(node1 == node2):
                temp_pair.first = node1
                temp_pair.second = node1
                undirected_G.edges[temp_pair] = info["GMEL"]
            else:
                # save the two label edges to simplify future access
                temp_pair.first = node1
                temp_pair.second = node2
                undirected_G.edges[temp_pair] = info["GMEL"]
                temp_pair.first = node2
                temp_pair.second = node1
                undirected_G.edges[temp_pair] = info["GMEL"]
        # prepare edges of undirected H
        for (node1, node2, info) in encoded_graphs[1].edges(data = True):
            if(node1 == node2):
                temp_pair.first = node1
                temp_pair.second = node1
                undirected_H.edges[temp_pair] = info["GMEL"]
            else:
                # save the two label edges to simplify future access
                temp_pair.first = node1
                temp_pair.second = node2
                undirected_H.edges[temp_pair] = info["GMEL"]
                temp_pair.first = node2
                temp_pair.second = node1
                undirected_H.edges[temp_pair] = info["GMEL"]

    # prepare neighbors
    if(nx.is_directed(nx_G)):
        directed_G.in_neighbors = {node:list(encoded_graphs[0].predecessors(node)) for node in list(encoded_graphs[0].nodes())}
        directed_H.in_neighbors = {node:list(encoded_graphs[1].predecessors(node)) for node in list(encoded_graphs[1].nodes())}
        directed_G.out_neighbors = {node:list(encoded_graphs[0].neighbors(node)) for node in list(encoded_graphs[0].nodes())}
        directed_H.out_neighbors = {node:list(encoded_graphs[1].neighbors(node)) for node in list(encoded_graphs[1].nodes())}
    else:
        undirected_G.neighbors = {node:list(encoded_graphs[0].neighbors(node)) for node in list(encoded_graphs[0].nodes())}
        undirected_H.neighbors = {node:list(encoded_graphs[1].neighbors(node)) for node in list(encoded_graphs[1].nodes())}

    # get total order for VF2-like analysis
    # NOTE: for connected extensions this should be given by concentric neighborhoods around the anchor with a multisource BFS
    undirected_copy_H = deepcopy(encoded_graphs[1])
    if(nx.is_directed(nx_H)):
        undirected_copy_H = undirected_copy_H.to_undirected()
    connectivity_neighbors = {node:list(undirected_copy_H.neighbors(node)) for node in list(undirected_copy_H.nodes())}
    current_level = [node2 for (node1, node2) in encoded_anchor]
    for each_pair in connectivity_neighbors:
        node = each_pair.first
        visited[node] = False

    for node in current_level:
        counter = counter + 1
        total_order[node] = counter
        visited[node] = True

    while(not current_level.empty()):
        # reinitialize next_level
        next_level.clear()
        # iterate getting immediate neighbors not already ordered
        for node1 in current_level:
            for node2 in connectivity_neighbors[node1]:
                # only assign oreder of not visited yet
                if(not visited[node2]):
                    # increase and assign counter
                    counter = counter + 1
                    total_order[node2] = counter
                    visited[node2] = True
                    # level management
                    next_level.push_back(node2)
        # update nodes to be ordered
        current_level.clear()
        current_level = next_level

    for each_pair in connectivity_neighbors:
        node = each_pair.first
        if(not visited[node]):
            counter = counter + 1
            total_order[node] = counter
            visited[node] = True

    # get expected order
    expected_order = min([nx_G.order(), nx_H.order()])

    # set recursion limit if recursive version was requested
    if(not iterative_search):
        scalation_value = 1.5
        required_limit = max([nx_G.order(), nx_H.order()])
        current_limit = getrecursionlimit()
        if(current_limit < (scalation_value * required_limit)):
            setrecursionlimit(int(scalation_value * required_limit))

    # get maximum extensions
    if(nx.is_directed(nx_G)):
        directed_maximum_connected_extensions(all_extensions,
                                              expected_order,
                                              encoded_anchor,
                                              directed_G,
                                              directed_H,
                                              total_order,
                                              encoded_extensions)
    else:
        if(iterative_search):
            undirected_maximum_connected_extensions_iterative(all_extensions,
                                                              expected_order,
                                                              encoded_anchor,
                                                              total_order,
                                                              undirected_G,
                                                              undirected_H,
                                                              encoded_extensions)
        else:
            undirected_maximum_connected_extensions_recursive(all_extensions,
                                                              expected_order,
                                                              encoded_anchor,
                                                              total_order,
                                                              undirected_G,
                                                              undirected_H,
                                                              encoded_extensions)
    # decode maximum extensions
    for each_extension in encoded_extensions:
        extensions.append(decode_match(list(each_extension), encoded_node_names))

    # check if the anchor was a good partial map
    if((len(extensions[0]) == nx_G.order()) and (len(extensions[0]) == nx_H.order())):
        good_anchor = True

    # end of function
    return(extensions, good_anchor)




# functions - maximum connected extensions - undirected ########################




# function: core routine of VF2-like undirected approach - iterative -----------
# NOTE: an iterative DFS version of this algorithm can be implemented without a "visited"
# list, since the total order given to the vertices of the second graph guarantess that
# the search space is actually a search tree, and thus cannot have repeated states.
cdef void undirected_maximum_connected_extensions_iterative(cpp_bool all_extensions,
                                                            size_t expected_order,
                                                            cpp_vector[cpp_pair[int, int]] input_anchor,
                                                            cpp_unordered_map[int, int] & total_order,
                                                            partial_maps_undirected_graph & G,
                                                            partial_maps_undirected_graph & H,
                                                            cpp_vector[cpp_vector[cpp_pair[int, int]]] & all_matches) noexcept:

    # local variables (cython)
    cdef size_t new_score = 0
    cdef size_t old_score = 0
    cdef cpp_bool semantic_feasibility = False
    cdef cpp_bool syntactic_feasibility = False
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[int] current_match_G
    cdef cpp_vector[int] current_match_H
    cdef cpp_vector[cpp_pair[int, int]] new_match
    cdef cpp_vector[cpp_pair[int, int]] candidates
    cdef cpp_vector[cpp_pair[int, int]] current_match
    cdef cpp_unordered_map[int, int] forward_match
    cdef cpp_stack[cpp_vector[cpp_pair[int, int]]] dfs_stack

    # initialize stack with (non-empty) input anchor match
    dfs_stack.push(input_anchor)

    # iterative DFS search
    while(not dfs_stack.empty()):
        # get current state and pop it
        current_match.clear()
        current_match = dfs_stack.top()
        dfs_stack.pop()

        # test initial match and improve if possible
        if(all_matches.empty()):
            if(not current_match.empty()):
                all_matches.push_back(current_match)
                # if complete match and only one then return
                if(G.nodes.size() == H.nodes.size()):
                    if(current_match.size() == G.nodes.size()):
                        if(not all_extensions):
                            return
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
                # if complete match and only one then return
                if(G.nodes.size() == H.nodes.size()):
                    if(current_match.size() == G.nodes.size()):
                        if(not all_extensions):
                            return

        # if not optimal yet then obtain available pairs
        if(current_match.size() < expected_order):
            # reinitialize variables used for filtering candidates
            current_match_G.clear()
            current_match_H.clear()
            forward_match.clear()
            candidates.clear()

            # generate auxiliary structures
            for each_pair in current_match:
                current_match_G.push_back(each_pair.first)
                current_match_H.push_back(each_pair.second)
                forward_match[each_pair.first] = each_pair.second

            # get candidate pairs
            candidates = undirected_candidates(current_match,
                                               current_match_G,
                                               current_match_H,
                                               G.neighbors,
                                               H.neighbors,
                                               total_order)

            # evaluate candidates
            for each_pair in candidates:
                # evaluate sintactic feasibility
                syntactic_feasibility = undirected_syntactic_feasibility(each_pair.first,
                                                                         each_pair.second,
                                                                         current_match_G,
                                                                         current_match_H,
                                                                         forward_match,
                                                                         G.neighbors,
                                                                         H.neighbors)

                if(syntactic_feasibility):
                    # evaluate semantic feasibility
                    semantic_feasibility = undirected_semantic_feasibility(each_pair.first,
                                                                           each_pair.second,
                                                                           current_match_G,
                                                                           forward_match,
                                                                           G.nodes,
                                                                           H.nodes,
                                                                           G.neighbors,
                                                                           G.edges,
                                                                           H.edges)

                    if(semantic_feasibility):
                        # build new match
                        new_match.clear()
                        new_match = current_match
                        new_match.push_back(each_pair)
                        # add new valid candidate states
                        dfs_stack.push(new_match)
    # end of function




# function: core routine of VF2-like undirected approach - recursive -----------
cdef void undirected_maximum_connected_extensions_recursive(cpp_bool all_extensions,
                                                            size_t expected_order,
                                                            cpp_vector[cpp_pair[int, int]] current_match,
                                                            cpp_unordered_map[int, int] & total_order,
                                                            partial_maps_undirected_graph & G,
                                                            partial_maps_undirected_graph & H,
                                                            cpp_vector[cpp_vector[cpp_pair[int, int]]] & all_matches) noexcept:

    # local variables (cython)
    cdef size_t new_score = 0
    cdef size_t old_score = 0
    cdef cpp_bool semantic_feasibility = False
    cdef cpp_bool syntactic_feasibility = False
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[int] current_match_G
    cdef cpp_vector[int] current_match_H
    cdef cpp_vector[cpp_pair[int, int]] new_match
    cdef cpp_vector[cpp_pair[int, int]] candidates
    cdef cpp_unordered_map[int, int] forward_match

    # test initial match and improve if possible
    if(all_matches.empty()):
        if(not current_match.empty()):
            all_matches.push_back(current_match)
            # if complete match and only one then return
            if(G.nodes.size() == H.nodes.size()):
                if(current_match.size() == G.nodes.size()):
                    if(not all_extensions):
                        return
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
            # if complete match and only one then return
            if(G.nodes.size() == H.nodes.size()):
                if(current_match.size() == G.nodes.size()):
                    if(not all_extensions):
                        return

    # if not optimal yet then obtain available pairs
    if(current_match.size() < expected_order):
        # generate auxiliary structures
        for each_pair in current_match:
            current_match_G.push_back(each_pair.first)
            current_match_H.push_back(each_pair.second)
            forward_match[each_pair.first] = each_pair.second

        # get candidate pairs
        candidates = undirected_candidates(current_match,
                                           current_match_G,
                                           current_match_H,
                                           G.neighbors,
                                           H.neighbors,
                                           total_order)

        # evaluate candidates
        for each_pair in candidates:
            # evaluate sintactic feasibility
            syntactic_feasibility = undirected_syntactic_feasibility(each_pair.first,
                                                                     each_pair.second,
                                                                     current_match_G,
                                                                     current_match_H,
                                                                     forward_match,
                                                                     G.neighbors,
                                                                     H.neighbors)
            if(syntactic_feasibility):
                # evaluate semantic feasibility
                semantic_feasibility = undirected_semantic_feasibility(each_pair.first,
                                                                       each_pair.second,
                                                                       current_match_G,
                                                                       forward_match,
                                                                       G.nodes,
                                                                       H.nodes,
                                                                       G.neighbors,
                                                                       G.edges,
                                                                       H.edges)
                if(semantic_feasibility):
                    # build new match
                    new_match.clear()
                    new_match = current_match
                    new_match.push_back(each_pair)
                    # extend match
                    undirected_maximum_connected_extensions_recursive(all_extensions,
                                                                      expected_order,
                                                                      new_match,
                                                                      total_order,
                                                                      G,
                                                                      H,
                                                                      all_matches)
                    # finish if only one complete extension was requested and it was already found
                    if(G.nodes.size() == H.nodes.size()):
                        if(not all_extensions):
                            # anchor is always present in vector of all matches at this point
                            if(all_matches[0].size() == G.nodes.size()):
                                return
    # end of function




# function: get candidate pairs for undirected extension search ----------------
cdef cpp_vector[cpp_pair[int, int]] undirected_candidates(cpp_vector[cpp_pair[int, int]] & current_match,
                                                          cpp_vector[int] & current_match_G,
                                                          cpp_vector[int] & current_match_H,
                                                          cpp_unordered_map[int, cpp_vector[int]] & neigh_G,
                                                          cpp_unordered_map[int, cpp_vector[int]] & neigh_H,
                                                          cpp_unordered_map[int, int] & total_order) noexcept:

    # output holders
    cdef cpp_vector[cpp_pair[int, int]] candidate_pairs

    # local variables
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int reference_maximum = 0
    cdef cpp_pair[int, int] each_pair
    cdef cpp_pair[int, int] temp_pair
    cdef cpp_vector[int] valid_G
    cdef cpp_vector[int] valid_H

    # get maximum value of total order in match
    for node in current_match_H:
        if(total_order[node] > reference_maximum):
            reference_maximum = total_order[node]

    # get candidate pairs based on valid sets
    for each_pair in current_match:
        # reinitialize valid neighbors
        valid_G.clear()
        valid_H.clear()

        # get valid neighbors in G
        for node in neigh_G[each_pair.first]:
            # if not yet in match
            if(find(current_match_G.begin(), current_match_G.end(), node) == current_match_G.end()):
                valid_G.push_back(node)

        # get valid neighbors in H
        for node in neigh_H[each_pair.second]:
            # if total order greater than in match
            if(total_order[node] > reference_maximum):
                # if not yet in match
                if(find(current_match_H.begin(), current_match_H.end(), node) == current_match_H.end()):
                    valid_H.push_back(node)

        # make product of valid neighbors
        if((not valid_G.empty()) and (not valid_H.empty())):
            for node1 in valid_G:
                for node2 in valid_H:
                    temp_pair.first = node1
                    temp_pair.second = node2
                    if(find(candidate_pairs.begin(), candidate_pairs.end(), temp_pair) == candidate_pairs.end()):
                        candidate_pairs.push_back(temp_pair)

    # end of function
    return(candidate_pairs)




# function: evaluate syntactic feasability for undirected extension ------------
cdef cpp_bool undirected_syntactic_feasibility(int node1,
                                               int node2,
                                               cpp_vector[int] & current_match_G,
                                               cpp_vector[int] & current_match_H,
                                               cpp_unordered_map[int, int] & forward_match,
                                               cpp_unordered_map[int, cpp_vector[int]] & neigh_G,
                                               cpp_unordered_map[int, cpp_vector[int]] & neigh_H) noexcept:

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
cdef cpp_bool undirected_semantic_feasibility(int node1,
                                              int node2,
                                              cpp_vector[int] & current_match_G,
                                              cpp_unordered_map[int, int] & forward_match,
                                              cpp_unordered_map[int, int] & nodes_G,
                                              cpp_unordered_map[int, int] & nodes_H,
                                              cpp_unordered_map[int, cpp_vector[int]] & neigh_G,
                                              cpp_map[cpp_pair[int, int], int] & edges_G,
                                              cpp_map[cpp_pair[int, int], int] & edges_H) noexcept:

    # local variables
    cdef int node = 0
    cdef cpp_pair[int, int] labeled_edge_G
    cdef cpp_pair[int, int] labeled_edge_H

    # compare vertex-labels
    if(nodes_G[node1] != nodes_H[node2]):
        return(False)

    # compare loop-labels
    if(find(neigh_G[node1].begin(), neigh_G[node1].end(), node1) != neigh_G[node1].end()):
        # loop in G
        labeled_edge_G.first = node1
        labeled_edge_G.second = node1
        # loop in H
        labeled_edge_H.first = node2
        labeled_edge_H.second = node2
        # compare edge labels
        if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
            return(False)

    # compare non-loop edge-labels
    for node in neigh_G[node1]:
        if(find(current_match_G.begin(), current_match_G.end(), node) != current_match_G.end()):
            # edge in G with only one end in match
            labeled_edge_G.first = node1
            labeled_edge_G.second = node
            # edge in H with only one end in match
            labeled_edge_H.first = node2
            labeled_edge_H.second = forward_match[node]
            # compare edge labels
            if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
                return(False)

    # end of function
    return(True)




# functions - maximum connected extensions - directed ##########################




# function: core routine of VF2-like directed approach -------------------------
cdef void directed_maximum_connected_extensions(cpp_bool all_extensions,
                                                size_t expected_order,
                                                cpp_vector[cpp_pair[int, int]] current_match,
                                                partial_maps_directed_graph & G,
                                                partial_maps_directed_graph & H,
                                                cpp_unordered_map[int, int] & total_order,
                                                cpp_vector[cpp_vector[cpp_pair[int, int]]] & all_matches) noexcept:

    # local variables
    cdef size_t new_score = 0
    cdef size_t old_score = 0
    cdef cpp_bool semantic_feasibility = False
    cdef cpp_bool syntactic_feasibility = False
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[int] current_match_G
    cdef cpp_vector[int] current_match_H
    cdef cpp_vector[cpp_pair[int, int]] new_match
    cdef cpp_vector[cpp_pair[int, int]] candidates
    cdef cpp_unordered_map[int, int] forward_match

    # test initial match and consecutive matches
    if(all_matches.empty()):
        if(not current_match.empty()):
            all_matches.push_back(current_match)
            # if complete match and only one then return
            if(G.nodes.size() == H.nodes.size()):
                if(current_match.size() == G.nodes.size()):
                    if(not all_extensions):
                        return
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
            # if complete match and only one then return
            if(G.nodes.size() == H.nodes.size()):
                if(current_match.size() == G.nodes.size()):
                    if(not all_extensions):
                        return

    # if not optimal yet then obtain available pairs
    if(current_match.size() < expected_order):
        # generate auxiliary structures
        for each_pair in current_match:
            current_match_G.push_back(each_pair.first)
            current_match_H.push_back(each_pair.second)
            forward_match[each_pair.first] = each_pair.second

        # get candidate pairs
        candidates = directed_candidates(current_match,
                                         current_match_G,
                                         current_match_H,
                                         G.in_neighbors,
                                         H.in_neighbors,
                                         G.out_neighbors,
                                         H.out_neighbors,
                                         total_order)

        # evaluate candidates
        for each_pair in candidates:
            # evaluate sintactic feasibility
            syntactic_feasibility = directed_syntactic_feasibility(each_pair.first,
                                                                   each_pair.second,
                                                                   current_match_G,
                                                                   current_match_H,
                                                                   forward_match,
                                                                   G.in_neighbors,
                                                                   H.in_neighbors,
                                                                   G.out_neighbors,
                                                                   H.out_neighbors)
            if(syntactic_feasibility):
                # evaluate semantic feasibility
                semantic_feasibility = directed_semantic_feasibility(each_pair.first,
                                                                     each_pair.second,
                                                                     current_match_G,
                                                                     forward_match,
                                                                     G.nodes,
                                                                     H.nodes,
                                                                     G.in_neighbors,
                                                                     G.out_neighbors,
                                                                     G.edges,
                                                                     H.edges)
                if(semantic_feasibility):
                    # build new match
                    new_match.clear()
                    new_match = current_match
                    new_match.push_back(each_pair)
                    # extend match
                    directed_maximum_connected_extensions(all_extensions,
                                                          expected_order,
                                                          new_match,
                                                          G,
                                                          H,
                                                          total_order,
                                                          all_matches)
                    # finish if only one complete extension was requested and it was already found
                    if(G.nodes.size() == H.nodes.size()):
                        if(not all_extensions):
                            # anchor is always present in vector of all matches at this point
                            if(all_matches[0].size() == G.nodes.size()):
                                return
    # end of function




# function: get candidate pairs for directed extension search ------------------
cdef cpp_vector[cpp_pair[int, int]] directed_candidates(cpp_vector[cpp_pair[int, int]] & current_match,
                                                        cpp_vector[int] & current_match_G,
                                                        cpp_vector[int] & current_match_H,
                                                        cpp_unordered_map[int, cpp_vector[int]] & in_neigh_G,
                                                        cpp_unordered_map[int, cpp_vector[int]] & in_neigh_H,
                                                        cpp_unordered_map[int, cpp_vector[int]] & out_neigh_G,
                                                        cpp_unordered_map[int, cpp_vector[int]] & out_neigh_H,
                                                        cpp_unordered_map[int, int] & total_order) noexcept:
    # output holders
    cdef cpp_vector[cpp_pair[int, int]] candidate_pairs

    # local variables
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int reference_maximum = 0
    cdef cpp_pair[int, int] temp_pair
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[int] valid_G
    cdef cpp_vector[int] valid_H

    # get maximum value of total order in match
    for node in current_match_H:
        if(total_order[node] > reference_maximum):
            reference_maximum = total_order[node]

    # get candidates based on valid sets
    for each_pair in current_match:
        # reinitialize valid neighbors
        valid_G.clear()
        valid_H.clear()

        # get valid in-neighbors in G
        for node in in_neigh_G[each_pair.first]:
            # if not yet in match
            if(find(current_match_G.begin(), current_match_G.end(), node) == current_match_G.end()):
                valid_G.push_back(node)

        # get valid in-neighbors in H
        for node in in_neigh_H[each_pair.second]:
            # if total order greater than in match
            if(total_order[node] > reference_maximum):
                # if not yet in match
                if(find(current_match_H.begin(), current_match_H.end(), node) == current_match_H.end()):
                    valid_H.push_back(node)

        # make product of valid in-neighbors
        if((not valid_G.empty()) and (not valid_H.empty())):
            for node1 in valid_G:
                for node2 in valid_H:
                    temp_pair.first = node1
                    temp_pair.second = node2
                    if(find(candidate_pairs.begin(), candidate_pairs.end(), temp_pair) == candidate_pairs.end()):
                        candidate_pairs.push_back(temp_pair)

        # reinitialize valid neighbors
        valid_G.clear()
        valid_H.clear()

        # get valid out-neighbors in G
        for node in out_neigh_G[each_pair.first]:
            # if not yet in match
            if(find(current_match_G.begin(), current_match_G.end(), node) == current_match_G.end()):
                valid_G.push_back(node)

        # get valid out-neighbors in H
        for node in out_neigh_H[each_pair.second]:
            # if total order greater than in match
            if(total_order[node] > reference_maximum):
                # if not yet in match
                if(find(current_match_H.begin(), current_match_H.end(), node) == current_match_H.end()):
                    valid_H.push_back(node)

        # make product of valid out-neighbors
        if((not valid_G.empty()) and (not valid_H.empty())):
            for node1 in valid_G:
                for node2 in valid_H:
                    temp_pair.first = node1
                    temp_pair.second = node2
                    if(find(candidate_pairs.begin(), candidate_pairs.end(), temp_pair) == candidate_pairs.end()):
                        candidate_pairs.push_back(temp_pair)

    # end of function
    return(candidate_pairs)




# function: evaluate syntactic feasability for directed extension --------------
cdef cpp_bool directed_syntactic_feasibility(int node1,
                                             int node2,
                                             cpp_vector[int] & current_match_G,
                                             cpp_vector[int] & current_match_H,
                                             cpp_unordered_map[int, int] & forward_match,
                                             cpp_unordered_map[int, cpp_vector[int]] & in_neigh_G,
                                             cpp_unordered_map[int, cpp_vector[int]] & in_neigh_H,
                                             cpp_unordered_map[int, cpp_vector[int]] & out_neigh_G,
                                             cpp_unordered_map[int, cpp_vector[int]] & out_neigh_H) noexcept:

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
cdef cpp_bool directed_semantic_feasibility(int node1,
                                            int node2,
                                            cpp_vector[int] & current_match_G,
                                            cpp_unordered_map[int, int] & forward_match,
                                            cpp_unordered_map[int, int] & nodes_G,
                                            cpp_unordered_map[int, int] & nodes_H,
                                            cpp_unordered_map[int, cpp_vector[int]] & in_neigh_G,
                                            cpp_unordered_map[int, cpp_vector[int]] & out_neigh_G,
                                            cpp_map[cpp_pair[int, int], int] & edges_G,
                                            cpp_map[cpp_pair[int, int], int] & edges_H) noexcept:

    # local variables
    cdef int node = 0
    cdef cpp_pair[int, int] labeled_edge_G
    cdef cpp_pair[int, int] labeled_edge_H

    # compare vertex-labels
    if(nodes_G[node1] != nodes_H[node2]):
        return(False)

    # compare loop-labels
    if(find(in_neigh_G[node1].begin(), in_neigh_G[node1].end(), node1) != in_neigh_G[node1].end()):
        # loop in G
        labeled_edge_G.first = node1
        labeled_edge_G.second = node1
        # loop in H
        labeled_edge_H.first = node2
        labeled_edge_H.second = node2
        # compare edge labels
        if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
            return(False)

    # compare non-loop in-edge-labels
    for node in in_neigh_G[node1]:
        if(find(current_match_G.begin(), current_match_G.end(), node) != current_match_G.end()):
            # edge in G with only one end in match
            labeled_edge_G.first = node
            labeled_edge_G.second = node1
            # edge in H with only one end in match
            labeled_edge_H.first = forward_match[node]
            labeled_edge_H.second = node2
            # compare edge labels
            if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
                return(False)

    # compare non-loop out-edge-labels
    for node in out_neigh_G[node1]:
        if(find(current_match_G.begin(), current_match_G.end(), node) != current_match_G.end()):
            # edge in G with only one end in match
            labeled_edge_G.first = node1
            labeled_edge_G.second = node
            # edge in H with only one end in match
            labeled_edge_H.first = node2
            labeled_edge_H.second = forward_match[node]
            # compare edge labels
            if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
                return(False)

    # end of function
    return(True)




################################################################################
################################################################################

################################################################################
#                                                                              #
# - GranMapache: GRAphs-and-Networks MAPping Applications                      #
#   with Cython and HEuristics                                                 #
#                                                                              #
# - Module: partial_maps                                                       #
#                                                                              #
# - Description: analysis of properties of partial maps, like maximum induced  #
#   connected extension(s), overlaps, consistency, and others.                 #
#                                                                              #
################################################################################




# C++ string encoding ##########################################################
# cython: c_string_type=unicode, c_string_encoding=utf8




# dependencies #################################################################




# already in python ------------------------------------------------------------
from copy import deepcopy
from sys import getrecursionlimit, setrecursionlimit




# not in python ----------------------------------------------------------------
import networkx as nx




# cython specifics -------------------------------------------------------------
import cython
from libcpp cimport bool as cpp_bool
from libcpp.pair cimport pair as cpp_pair
from libcpp.stack cimport stack as cpp_stack
from libcpp.vector cimport vector as cpp_vector
from libcpp.string cimport string as cpp_string
from libcpp.unordered_set cimport unordered_set as cpp_unordered_set
from libcpp.unordered_map cimport unordered_map as cpp_unordered_map
cdef extern from "<algorithm>" namespace "std":
    # find element in vector
    cpp_string to_string(int value)




# custom dependencies ----------------------------------------------------------
from .integerization import encode_graphs, decode_graphs, encode_match, decode_match




# C/C++ structs ################################################################




# struct: undirected graph ---------------------------------------------------
# using unordered maps and unordered sets for scalability
cdef struct partial_maps_undirected_graph:
    cpp_unordered_map[int, int] nodes
    cpp_unordered_map[cpp_string, int] edges
    cpp_unordered_map[int, cpp_unordered_set[int]] neighbors




# struct: directed graph -------------------------------------------------------
# using unordered maps and unordered sets for scalability
cdef struct partial_maps_directed_graph:
    cpp_unordered_map[int, int] nodes
    cpp_unordered_map[cpp_string, int] edges
    cpp_unordered_map[int, cpp_unordered_set[int]] in_neighbors
    cpp_unordered_map[int, cpp_unordered_set[int]] out_neighbors




# structure: container for information of candidate matches undirected ---------
# using unordered maps and unordered sets for scalability
cdef struct partial_maps_candidates_struct_undirected:
    cpp_vector[cpp_pair[int, int]] candidates
    cpp_unordered_set[int] ring_G
    cpp_unordered_set[int] ring_H




# structure: container for information of candidate matches directed -----------
# using unordered maps and unordered sets for scalability
cdef struct partial_maps_candidates_struct_directed:
    cpp_vector[cpp_pair[int, int]] candidates
    cpp_unordered_set[int] in_ring_G
    cpp_unordered_set[int] in_ring_H
    cpp_unordered_set[int] out_ring_G
    cpp_unordered_set[int] out_ring_H




# algorithms ###################################################################




# functions - induced connected extensions - wrapper ###########################




# function: callable wrapper for the induced connected extensions --------------
def induced_connected_extensions(nx_G = nx.Graph(),          # can be nx.DiGraph
                                 nx_H = nx.Graph(),          # can be nx.DiGraph
                                 input_anchor = [],          # should be non-empty list
                                 node_labels = True,         # consider node labels when evaluating the extensions
                                 edge_labels = True,         # consider edge labels when evaluating the extensions
                                 all_extensions = False,     # by default stops when finding one complete extension (if any)
                                 iterative_search = True):   # by default an iterative search is used, otherwise a recursive version is called
    # description
    """
    > description: receives two networkx graphs G and H of the same order, and a match
    between them (here called anchor), and uses a variant of the VF2 algorithm by doing
    an isomorphism-like search to obtain the induced complete connected extensions of the
    anchor, if any. This exists if and only if (1) the input graphs are "balanced" (with the
    same number of nodes of each label), (2) if the provided partial map is a good atom map,
    i.e., if it already covers all the changing edges, and (3) if the ITS that this map
    induces is connected. Such extension, if it exists, is unique up to equivalence of atom
    maps, and by default the function will stop when finding it. This can be changed with
    the boolean parameter all_extensions, but it should be noted that such search can be more
    time consuming depending on the number of automotphisms of the input graphs.

    > input:
    * nx_G - first networkx (di)graph being matched.
    * nx_H - second networkx (di)graph being matched.
    * input_anchor - inyective map as a non-empty list of 2-tuples (x, y) of nodes x
    from G and y from H. An exception is raised if the anchor is empty.
    * node_labels - boolean indicating if node labels should be considered for the search,
    which is the default behavior, or if they should be ignored.
    * edge_labels - boolean indicating if edge labels should be considered for the search,
    which is the default behavior, or if they should be ignored.
    * all_extensions - boolean indicating if the function should stop as soon as one complete
    extension is found (if any) - this is the default behavior- or if it should search for all
    possible (complete) extensions. NOTE: mathematically speaking, the complete extension of
    a FIXED and GOOD anchor is unique up to equivalence of bijections, that is, equivalence
    of atom maps or correspondingly isomorphism of their ITS graphs. Nontheless it should be
    noted that changing the anchor may produce a non-equivalent (complete) extension. In other
    words, each call to this function can (mathematically) produce only one complete extension,
    but calls with different anchors can produce non-equivalent extensions, even if the anchors
    themselves produce isomorphic partial ITS graphs.
    * iterative_search - boolean indicating if the iterative version of this algorithm should
    be used (the default), or if a recursive version of it should be used instead.

    > output:
    * extensions - possibly empty list of injective maps each as a list of 2-tuples (x, y) of
    nodes x from G and y from H, each representing the induced complete connected extensions
    of the anchor (each extension contains the anchor as a sublist), if it exists.
    * good_anchor - boolean value indicating if the extensions where found, i.e., they cover
    all nodes of G and of H, and thus if they are bijections between G and H. If so, the anchor
    is what we refered to as a "good partial atom map", and equivalenteĺy the match obtained when
    removing the anchhor from any such extension is a graph-isomorphism between the "remainder"
    graphs it induces from G and H.

    > calls:
    * .integerization.encode_graphs
    * .integerization.decode_graphs
    * .integerization.encode_match
    * .integerization.decode_match
    * undirected_induced_connected_extensions_iterative
    * undirected_induced_connected_extensions_recursive
    * directed_induced_connected_extensions_iterative
    * directed_induced_connected_extensions_recursive
    """

    # exception handling and input correctness
    test_list = [0, 0]
    test_tuple = (0, 0)
    test_undir = nx.Graph()
    test_dir = nx.DiGraph()
    test_bool = False
    # check that first argument is networkx graph or digraph
    if(type(nx_G) not in [type(test_undir), type(test_dir)]):
        raise(ValueError("gmapache: first argument must be a networkx graph or digraph."))
    # check that second argument is networkx graph or digraph
    if(type(nx_H) not in [type(test_undir), type(test_dir)]):
        raise(ValueError("gmapache: second argument must be a networkx graph or digraph."))
    # check that the input graphs have the same type
    if((nx.is_directed(nx_G)) and (not nx.is_directed(nx_H))):
        raise(ValueError("gmapache: input graphs must be both directed or both undirected."))
    if((not nx.is_directed(nx_G)) and (nx.is_directed(nx_H))):
        raise(ValueError("gmapache: input graphs must be both directed or both undirected."))
    # check that input graphs are not null graphs
    if((nx_G.order() == 0) or (nx_H.order() == 0)):
        raise(ValueError("gmapache: input graphs must have at least one node each."))
    # check that input graphs have the same number of vertices
    if(not nx_G.order() == nx_H.order()):
        raise(ValueError("gmapache: input graphs must have the same number of vertices."))
    # check that fourth argument is a boolean variable
    if(type(node_labels) not in [type(test_bool)]):
        raise(ValueError("gmapache: fourth argument must be a boolean variable."))
    # check that fifth argument is a boolean variable
    if(type(edge_labels) not in [type(test_bool)]):
        raise(ValueError("gmapache: fifth argument must be a boolean variable."))
    # check that sixth argument is a boolean variable
    if(type(all_extensions) not in [type(test_bool)]):
        raise(ValueError("gmapache: sixth argument must be a boolean variable."))
    # check that seventh argument is a boolean variable
    if(type(iterative_search) not in [type(test_bool)]):
        raise(ValueError("gmapache: seventh argument must be a boolean variable."))
    # check that third argument is a list
    if(not type(input_anchor) in [type(test_list)]):
        raise(ValueError("gmapache: third argument must be a non-empty list of 2-tuples."))
    if(len(input_anchor) == 0):
        raise(ValueError("gmapache: third argument must be a non-empty list of 2-tuples."))
    # check correctness of entries in third argument
    for test_entry in input_anchor:
        if(not type(test_entry) in [type(test_tuple)]):
            raise(ValueError("gmapache: all elements in input list must be tuples."))
        if(not len(test_entry) == 2):
            raise(ValueError("gmapache: all tuples in input list must be of lenght 2."))
        if(test_entry[0] not in list(nx_G.nodes())):
            raise(ValueError("gmapache: the input list is matching a vertex not present in the first graph."))
        if(test_entry[1] not in list(nx_H.nodes())):
            raise(ValueError("gmapache: the input list is matching a vertex not present in the second graph."))
    # check amount of entries in third argument
    if(not len(list(set([x for (x, y) in input_anchor]))) == len(input_anchor)):
        raise(ValueError("gmapache: the input list must be an injective map and without repeated elements."))
    if(not len(list(set([y for (x, y) in input_anchor]))) == len(input_anchor)):
        raise(ValueError("gmapache: the input list must be an injective map and without repeated elements."))

    # output holders
    cdef list extensions = []
    good_anchor = False

    # local variables (cython)
    cdef int node = 0
    cdef int label = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int counter = 0
    cdef int current_limit = 0
    cdef int required_limit = 0
    cdef float scalation_value = 0
    cdef size_t expected_order = 0
    cdef cpp_string comma
    comma.push_back(44)
    cdef cpp_string temp_str
    cdef cpp_pair[int, int] node_and_label
    cdef cpp_pair[int, int] label_and_count
    cdef cpp_pair[int, cpp_unordered_set[int]] each_pair
    cdef cpp_vector[int] next_level
    cdef cpp_vector[int] current_level
    cdef cpp_vector[cpp_pair[int, int]] encoded_anchor
    cdef cpp_vector[cpp_pair[int, int]] each_extension
    cdef cpp_vector[cpp_vector[cpp_pair[int, int]]] encoded_extensions
    cdef cpp_unordered_map[int, int] total_order
    cdef cpp_unordered_map[int, cpp_bool] visited
    cdef cpp_unordered_map[int, int] count_node_labels_G
    cdef cpp_unordered_map[int, int] count_node_labels_H
    cdef cpp_unordered_map[int, cpp_unordered_set[int]] connectivity_neighbors
    cdef partial_maps_directed_graph directed_G
    cdef partial_maps_directed_graph directed_H
    cdef partial_maps_undirected_graph undirected_G
    cdef partial_maps_undirected_graph undirected_H

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

    # check that the graphs share the same number of nodes with each label (only if the search should preserve node labels)
    # NOTE: the direction is not important, this is just choosing the object that was initialize and is now being used
    if(node_labels):
        if(nx.is_directed(nx_G)):
            # get label count on directed G
            for node_and_label in directed_G.nodes:
                # get node label
                label = node_and_label.second
                # count node label
                if(count_node_labels_G.find(label) != count_node_labels_G.end()):
                    count_node_labels_G[label] = count_node_labels_G[label] + 1
                else:
                    count_node_labels_G[label] = 1
            # get label count on directed H
            for node_and_label in directed_H.nodes:
                # get node label
                label = node_and_label.second
                # count node label
                if(count_node_labels_H.find(label) != count_node_labels_H.end()):
                    count_node_labels_H[label] = count_node_labels_H[label] + 1
                else:
                    count_node_labels_H[label] = 1
        else:
            # get label count on undirected G
            for node_and_label in undirected_G.nodes:
                # get node label
                label = node_and_label.second
                # count node label
                if(count_node_labels_G.find(label) != count_node_labels_G.end()):
                    count_node_labels_G[label] = count_node_labels_G[label] + 1
                else:
                    count_node_labels_G[label] = 1
            # get label count on undirected H
            for node_and_label in undirected_H.nodes:
                # get node label
                label = node_and_label.second
                # count node label
                if(count_node_labels_H.find(label) != count_node_labels_H.end()):
                    count_node_labels_H[label] = count_node_labels_H[label] + 1
                else:
                    count_node_labels_H[label] = 1

        # compare node counts in both graphs per label
        if(count_node_labels_G.size() != count_node_labels_H.size()):
            raise(ValueError("gmapache: requested preservation of node labels but input graphs have different sets of node labels."))
        for label_and_count in count_node_labels_G:
            if(count_node_labels_H.find(label_and_count.first) == count_node_labels_H.end()):
                raise(ValueError("gmapache: requested preservation of node labels but input graphs have different sets of node labels."))
            else:
                if(label_and_count.second != count_node_labels_H[label_and_count.first]):
                    raise(ValueError("gmapache: requested preservation of node labels but input graphs have different numbers of nodes for some labels."))

    # prepare edges
    if(nx.is_directed(nx_G)):
        directed_G.edges = {str(node1)+","+str(node2):info["GMEL"] for (node1, node2, info) in encoded_graphs[0].edges(data = True)}
        directed_H.edges = {str(node1)+","+str(node2):info["GMEL"] for (node1, node2, info) in encoded_graphs[1].edges(data = True)}
    else:
        # prepare edges of undirected G
        for (node1, node2, info) in encoded_graphs[0].edges(data = True):
            if(node1 == node2):
                temp_str = to_string(node1) + comma + to_string(node1)
                undirected_G.edges[temp_str] = info["GMEL"]
            else:
                # save the two label edges to simplify future access
                temp_str = to_string(node1) + comma + to_string(node2)
                undirected_G.edges[temp_str] = info["GMEL"]
                temp_str = to_string(node2) + comma + to_string(node1)
                undirected_G.edges[temp_str] = info["GMEL"]
        # prepare edges of undirected H
        for (node1, node2, info) in encoded_graphs[1].edges(data = True):
            if(node1 == node2):
                temp_str = to_string(node1) + comma + to_string(node1)
                undirected_H.edges[temp_str] = info["GMEL"]
            else:
                # save the two label edges to simplify future access
                temp_str = to_string(node1) + comma + to_string(node2)
                undirected_H.edges[temp_str] = info["GMEL"]
                temp_str = to_string(node2) + comma + to_string(node1)
                undirected_H.edges[temp_str] = info["GMEL"]

    # prepare neighbors
    if(nx.is_directed(nx_G)):
        directed_G.in_neighbors = {node:set(encoded_graphs[0].predecessors(node)) for node in list(encoded_graphs[0].nodes())}
        directed_H.in_neighbors = {node:set(encoded_graphs[1].predecessors(node)) for node in list(encoded_graphs[1].nodes())}
        directed_G.out_neighbors = {node:set(encoded_graphs[0].neighbors(node)) for node in list(encoded_graphs[0].nodes())}
        directed_H.out_neighbors = {node:set(encoded_graphs[1].neighbors(node)) for node in list(encoded_graphs[1].nodes())}
    else:
        undirected_G.neighbors = {node:set(encoded_graphs[0].neighbors(node)) for node in list(encoded_graphs[0].nodes())}
        undirected_H.neighbors = {node:set(encoded_graphs[1].neighbors(node)) for node in list(encoded_graphs[1].nodes())}

    # get total order for VF2-like analysis
    # NOTE: for connected extensions this should be given by concentric neighborhoods around the anchor with a multisource BFS
    undirected_copy_H = deepcopy(encoded_graphs[1])
    if(nx.is_directed(nx_H)):
        undirected_copy_H = undirected_copy_H.to_undirected()
    connectivity_neighbors = {node:set(undirected_copy_H.neighbors(node)) for node in list(undirected_copy_H.nodes())}
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
    expected_order = nx_G.order()

    # set recursion limit if recursive version was requested
    if(not iterative_search):
        scalation_value = 1.5
        required_limit = max([nx_G.order(), nx_H.order()])
        current_limit = getrecursionlimit()
        if(current_limit < (scalation_value * required_limit)):
            setrecursionlimit(int(scalation_value * required_limit))

    # get induced extensions
    if(nx.is_directed(nx_G)):
        if(iterative_search):
            directed_induced_connected_extensions_iterative(node_labels,
                                                            edge_labels,
                                                            all_extensions,
                                                            expected_order,
                                                            encoded_anchor,
                                                            total_order,
                                                            directed_G,
                                                            directed_H,
                                                            encoded_extensions)
        else:
            directed_induced_connected_extensions_recursive(node_labels,
                                                            edge_labels,
                                                            all_extensions,
                                                            expected_order,
                                                            encoded_anchor,
                                                            total_order,
                                                            directed_G,
                                                            directed_H,
                                                            encoded_extensions)
    else:
        if(iterative_search):
            undirected_induced_connected_extensions_iterative(node_labels,
                                                              edge_labels,
                                                              all_extensions,
                                                              expected_order,
                                                              encoded_anchor,
                                                              total_order,
                                                              undirected_G,
                                                              undirected_H,
                                                              encoded_extensions)
        else:
            undirected_induced_connected_extensions_recursive(node_labels,
                                                              edge_labels,
                                                              all_extensions,
                                                              expected_order,
                                                              encoded_anchor,
                                                              total_order,
                                                              undirected_G,
                                                              undirected_H,
                                                              encoded_extensions)

    # decode induced extensions
    for each_extension in encoded_extensions:
        extensions.append(decode_match(list(each_extension), encoded_node_names))

    # check if the anchor was a good partial map
    if(len(extensions[0]) == nx_G.order()):
        good_anchor = True

    # end of function
    return(extensions, good_anchor)




# functions - induced connected extensions - undirected ########################




# function: core routine of VF2-like undirected approach - iterative -----------
# NOTE: an iterative DFS version of this algorithm can be implemented without a "visited"
# list, since the total order given to the vertices of the second graph guarantees that
# the search space is actually a search tree, and thus cannot have repeated states.
cdef void undirected_induced_connected_extensions_iterative(cpp_bool node_labels,
                                                            cpp_bool edge_labels,
                                                            cpp_bool all_extensions,
                                                            size_t expected_order,
                                                            cpp_vector[cpp_pair[int, int]] input_anchor,
                                                            cpp_unordered_map[int, int] & total_order,
                                                            partial_maps_undirected_graph & G,
                                                            partial_maps_undirected_graph & H,
                                                            cpp_vector[cpp_vector[cpp_pair[int, int]]] & all_matches) noexcept:

    # local variables (cython)
    cdef size_t new_score = 0
    cdef size_t old_score = 0
    cdef cpp_bool semantic_feasibility_res = True
    cdef cpp_bool syntactic_feasibility_res = True
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[cpp_pair[int, int]] new_match
    cdef cpp_vector[cpp_pair[int, int]] current_match
    cdef cpp_unordered_set[int] current_match_G
    cdef cpp_unordered_set[int] current_match_H
    cdef cpp_unordered_map[int, int] forward_match
    cdef cpp_stack[cpp_vector[cpp_pair[int, int]]] dfs_stack
    cdef partial_maps_candidates_struct_undirected candidates_struct

    # initialize stack with (non-empty) input anchor match
    dfs_stack.push(input_anchor)

    # iterative DFS search
    while(not dfs_stack.empty()):
        # get current state and pop it
        current_match.clear()
        current_match = dfs_stack.top()
        dfs_stack.pop()

        # save if complete extension was reached
        if(current_match.size() == expected_order):
            all_matches.push_back(current_match)
            if(not all_extensions):
                return

        # if not optimal yet then obtain available pairs
        if(current_match.size() < expected_order):
            # reinitialize variables used for filtering candidates
            current_match_G.clear()
            current_match_H.clear()
            forward_match.clear()

            # generate auxiliary structures
            for each_pair in current_match:
                current_match_G.insert(each_pair.first)
                current_match_H.insert(each_pair.second)
                forward_match[each_pair.first] = each_pair.second

            # get candidate pairs
            candidates_struct = undirected_candidates(expected_order,
                                                      current_match,
                                                      current_match_G,
                                                      current_match_H,
                                                      G.neighbors,
                                                      H.neighbors,
                                                      total_order)

            # evaluate candidates
            for each_pair in candidates_struct.candidates:
                # evaluate syntactic feasibility
                syntactic_feasibility_res = syntactic_feasibility(each_pair.first,
                                                                  each_pair.second,
                                                                  candidates_struct.ring_G,
                                                                  candidates_struct.ring_H,
                                                                  current_match_G,
                                                                  current_match_H,
                                                                  forward_match,
                                                                  G.neighbors,
                                                                  H.neighbors)

                if(syntactic_feasibility_res):
                    # evaluate semantic feasibility
                    if(node_labels or edge_labels):
                        semantic_feasibility_res = semantic_feasibility(node_labels,
                                                                        edge_labels,
                                                                        each_pair.first,
                                                                        each_pair.second,
                                                                        candidates_struct.ring_G,
                                                                        candidates_struct.ring_H,
                                                                        current_match_G,
                                                                        current_match_H,
                                                                        forward_match,
                                                                        G.nodes,
                                                                        H.nodes,
                                                                        G.neighbors,
                                                                        H.neighbors,
                                                                        G.edges,
                                                                        H.edges)

                    # push to stack if valid
                    if(semantic_feasibility_res):
                        # build new match
                        new_match.clear()
                        new_match = current_match
                        new_match.push_back(each_pair)
                        # add new valid candidate states
                        dfs_stack.push(new_match)
    # end of function




# function: core routine of VF2-like undirected approach - recursive -----------
cdef void undirected_induced_connected_extensions_recursive(cpp_bool node_labels,
                                                            cpp_bool edge_labels,
                                                            cpp_bool all_extensions,
                                                            size_t expected_order,
                                                            cpp_vector[cpp_pair[int, int]] current_match,
                                                            cpp_unordered_map[int, int] & total_order,
                                                            partial_maps_undirected_graph & G,
                                                            partial_maps_undirected_graph & H,
                                                            cpp_vector[cpp_vector[cpp_pair[int, int]]] & all_matches) noexcept:

    # local variables (cython)
    cdef size_t new_score = 0
    cdef size_t old_score = 0
    cdef cpp_bool semantic_feasibility_res = True
    cdef cpp_bool syntactic_feasibility_res = True
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[cpp_pair[int, int]] new_match
    cdef cpp_unordered_set[int] current_match_G
    cdef cpp_unordered_set[int] current_match_H
    cdef cpp_unordered_map[int, int] forward_match
    cdef partial_maps_candidates_struct_undirected candidates_struct

    # save if complete extension was reached
    if(current_match.size() == expected_order):
        all_matches.push_back(current_match)
        if(not all_extensions):
            return

    # if not optimal yet then obtain available pairs
    if(current_match.size() < expected_order):
        # generate auxiliary structures
        for each_pair in current_match:
            current_match_G.insert(each_pair.first)
            current_match_H.insert(each_pair.second)
            forward_match[each_pair.first] = each_pair.second

        # get candidate pairs
        candidates_struct = undirected_candidates(expected_order,
                                                  current_match,
                                                  current_match_G,
                                                  current_match_H,
                                                  G.neighbors,
                                                  H.neighbors,
                                                  total_order)

        # evaluate candidates
        for each_pair in candidates_struct.candidates:
            # evaluate syntactic feasibility
            syntactic_feasibility_res = syntactic_feasibility(each_pair.first,
                                                              each_pair.second,
                                                              candidates_struct.ring_G,
                                                              candidates_struct.ring_H,
                                                              current_match_G,
                                                              current_match_H,
                                                              forward_match,
                                                              G.neighbors,
                                                              H.neighbors)

            if(syntactic_feasibility_res):
                # evaluate semantic feasibility
                if(node_labels or edge_labels):
                    semantic_feasibility_res = semantic_feasibility(node_labels,
                                                                    edge_labels,
                                                                    each_pair.first,
                                                                    each_pair.second,
                                                                    candidates_struct.ring_G,
                                                                    candidates_struct.ring_H,
                                                                    current_match_G,
                                                                    current_match_H,
                                                                    forward_match,
                                                                    G.nodes,
                                                                    H.nodes,
                                                                    G.neighbors,
                                                                    H.neighbors,
                                                                    G.edges,
                                                                    H.edges)

                # push to stack if valid
                if(semantic_feasibility_res):
                    # build new match
                    new_match.clear()
                    new_match = current_match
                    new_match.push_back(each_pair)
                    # extend match
                    undirected_induced_connected_extensions_recursive(node_labels,
                                                                      edge_labels,
                                                                      all_extensions,
                                                                      expected_order,
                                                                      new_match,
                                                                      total_order,
                                                                      G,
                                                                      H,
                                                                      all_matches)

                    # finish if only one extension was requested and it was already found
                    if(not all_matches.empty()):
                        if(not all_extensions):
                            return
    # end of function




# function: get candidate pairs for undirected extension search ----------------
cdef partial_maps_candidates_struct_undirected undirected_candidates(size_t expected_order,
                                                                     cpp_vector[cpp_pair[int, int]] & current_match,
                                                                     cpp_unordered_set[int] & current_match_G,
                                                                     cpp_unordered_set[int] & current_match_H,
                                                                     cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_G,
                                                                     cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_H,
                                                                     cpp_unordered_map[int, int] & total_order) noexcept:

    # output holders
    cdef partial_maps_candidates_struct_undirected candidates_struct

    # local variables
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int reference_minimum = 0
    cdef cpp_string comma
    comma.push_back(44)
    cdef cpp_string temp_string
    cdef cpp_pair[int, int] temp_pair
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[int] valid_G
    cdef cpp_vector[int] valid_H
    cdef cpp_vector[cpp_pair[int, int]] candidate_pairs
    cdef cpp_unordered_set[int] ring_G
    cdef cpp_unordered_set[int] ring_H
    cdef cpp_unordered_set[cpp_string] candidate_pairs_member

    # initialize reference minimum
    reference_minimum = 1 + 2*(<int>expected_order)

    # get valid sets
    for each_pair in current_match:
        # reinitialize valid neighbors
        valid_G.clear()
        valid_H.clear()

        # get valid neighbors in G
        for node in neigh_G[each_pair.first]:
            # if not yet in match
            if(current_match_G.find(node) == current_match_G.end()):
                valid_G.push_back(node)

        # get valid neighbors in H
        for node in neigh_H[each_pair.second]:
            # if not yet in match
            if(current_match_H.find(node) == current_match_H.end()):
                valid_H.push_back(node)

        # make product of valid neighbors
        if((not valid_G.empty()) and (not valid_H.empty())):
            for node1 in valid_G:
                for node2 in valid_H:
                    # check if pair satisfies minimum total order condition
                    if(total_order[node2] <= reference_minimum):
                        temp_string = to_string(node1) + comma + to_string(node2)
                        if(candidate_pairs_member.find(temp_string) == candidate_pairs_member.end()):
                            temp_pair.first = node1
                            temp_pair.second = node2
                            if(total_order[node2] == reference_minimum):
                                # add proper pair
                                candidate_pairs.push_back(temp_pair)
                                # add string version for constant look ups
                                candidate_pairs_member.insert(temp_string)
                            if(total_order[node2] < reference_minimum):
                                # update control values
                                reference_minimum = total_order[node2]
                                candidate_pairs.clear()
                                candidate_pairs_member.clear()
                                # add proper pair
                                candidate_pairs.push_back(temp_pair)
                                # add string version for constant look ups
                                candidate_pairs_member.insert(temp_string)

        # additionally save neighbors in ring around match in G
        for node in valid_G:
            # save if unrepeated
            if(ring_G.find(node) == ring_G.end()):
                ring_G.insert(node)

        # additionally save neighbors in ring around match in H
        for node in valid_H:
            # save if unrepeated
            if(ring_H.find(node) == ring_H.end()):
                ring_H.insert(node)

    # pack return structure
    candidates_struct.candidates = candidate_pairs
    candidates_struct.ring_G = ring_G
    candidates_struct.ring_H = ring_H
    # end of function
    return(candidates_struct)




# functions - induced connected extensions - undirected ########################




# function: core routine of VF2-like directed approach - iterative -------------
# NOTE: an iterative DFS version of this algorithm can be implemented without a "visited"
# list, since the total order given to the vertices of the second graph guarantees that
# the search space is actually a search tree, and thus cannot have repeated states.
cdef void directed_induced_connected_extensions_iterative(cpp_bool node_labels,
                                                          cpp_bool edge_labels,
                                                          cpp_bool all_extensions,
                                                          size_t expected_order,
                                                          cpp_vector[cpp_pair[int, int]] input_anchor,
                                                          cpp_unordered_map[int, int] & total_order,
                                                          partial_maps_directed_graph & G,
                                                          partial_maps_directed_graph & H,
                                                          cpp_vector[cpp_vector[cpp_pair[int, int]]] & all_matches) noexcept:

    # local variables (cython)
    cdef size_t new_score = 0
    cdef size_t old_score = 0
    cdef cpp_bool in_semantic_feasibility = True
    cdef cpp_bool in_syntactic_feasibility = True
    cdef cpp_bool out_semantic_feasibility = True
    cdef cpp_bool out_syntactic_feasibility = True
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[cpp_pair[int, int]] new_match
    cdef cpp_vector[cpp_pair[int, int]] current_match
    cdef cpp_unordered_set[int] current_match_G
    cdef cpp_unordered_set[int] current_match_H
    cdef cpp_unordered_map[int, int] forward_match
    cdef cpp_stack[cpp_vector[cpp_pair[int, int]]] dfs_stack
    cdef partial_maps_candidates_struct_directed candidates_struct

    # initialize stack with (non-empty) input anchor match
    dfs_stack.push(input_anchor)

    # iterative DFS search
    while(not dfs_stack.empty()):
        # get current state and pop it
        current_match.clear()
        current_match = dfs_stack.top()
        dfs_stack.pop()

        # save if complete extension was reached
        if(current_match.size() == expected_order):
            all_matches.push_back(current_match)
            if(not all_extensions):
                return

        # if not optimal yet then obtain available pairs
        if(current_match.size() < expected_order):
            # reinitialize variables used for filtering candidates
            current_match_G.clear()
            current_match_H.clear()
            forward_match.clear()

            # generate auxiliary structures
            for each_pair in current_match:
                current_match_G.insert(each_pair.first)
                current_match_H.insert(each_pair.second)
                forward_match[each_pair.first] = each_pair.second

            # get candidate pairs
            candidates_struct = directed_candidates(expected_order,
                                                    current_match,
                                                    current_match_G,
                                                    current_match_H,
                                                    G.in_neighbors,
                                                    H.in_neighbors,
                                                    G.out_neighbors,
                                                    H.out_neighbors,
                                                    total_order)

            # evaluate candidates
            for each_pair in candidates_struct.candidates:
                # evaluate syntactic feasibility of in-neighbors
                in_syntactic_feasibility = syntactic_feasibility(each_pair.first,
                                                                 each_pair.second,
                                                                 candidates_struct.in_ring_G,
                                                                 candidates_struct.in_ring_H,
                                                                 current_match_G,
                                                                 current_match_H,
                                                                 forward_match,
                                                                 G.in_neighbors,
                                                                 H.in_neighbors)

                if(in_syntactic_feasibility):
                    # evaluate syntactic feasibility of out-neighbors
                    out_syntactic_feasibility = syntactic_feasibility(each_pair.first,
                                                                      each_pair.second,
                                                                      candidates_struct.out_ring_G,
                                                                      candidates_struct.out_ring_H,
                                                                      current_match_G,
                                                                      current_match_H,
                                                                      forward_match,
                                                                      G.out_neighbors,
                                                                      H.out_neighbors)

                    if(out_syntactic_feasibility):
                        # evaluate semantic feasibility of in-neighors
                        if(node_labels or edge_labels):
                            in_semantic_feasibility = semantic_feasibility(node_labels,
                                                                           edge_labels,
                                                                           each_pair.first,
                                                                           each_pair.second,
                                                                           candidates_struct.in_ring_G,
                                                                           candidates_struct.in_ring_H,
                                                                           current_match_G,
                                                                           current_match_H,
                                                                           forward_match,
                                                                           G.nodes,
                                                                           H.nodes,
                                                                           G.in_neighbors,
                                                                           H.in_neighbors,
                                                                           G.edges,
                                                                           H.edges)

                        if(in_semantic_feasibility):
                            # evaluate semantic feasibility of out-neighors
                            if(node_labels or edge_labels):
                                out_semantic_feasibility = semantic_feasibility(node_labels,
                                                                                edge_labels,
                                                                                each_pair.first,
                                                                                each_pair.second,
                                                                                candidates_struct.out_ring_G,
                                                                                candidates_struct.out_ring_H,
                                                                                current_match_G,
                                                                                current_match_H,
                                                                                forward_match,
                                                                                G.nodes,
                                                                                H.nodes,
                                                                                G.out_neighbors,
                                                                                H.out_neighbors,
                                                                                G.edges,
                                                                                H.edges)

                            # push to stack if valid
                            if(out_semantic_feasibility):
                                # build new match
                                new_match.clear()
                                new_match = current_match
                                new_match.push_back(each_pair)
                                # add new valid candidate states
                                dfs_stack.push(new_match)
    # end of function




# function: core routine of VF2-like directed approach - recursive -------------
cdef void directed_induced_connected_extensions_recursive(cpp_bool node_labels,
                                                          cpp_bool edge_labels,
                                                          cpp_bool all_extensions,
                                                          size_t expected_order,
                                                          cpp_vector[cpp_pair[int, int]] current_match,
                                                          cpp_unordered_map[int, int] & total_order,
                                                          partial_maps_directed_graph & G,
                                                          partial_maps_directed_graph & H,
                                                          cpp_vector[cpp_vector[cpp_pair[int, int]]] & all_matches) noexcept:

    # local variables (cython)
    cdef size_t new_score = 0
    cdef size_t old_score = 0
    cdef cpp_bool in_semantic_feasibility = True
    cdef cpp_bool in_syntactic_feasibility = True
    cdef cpp_bool out_semantic_feasibility = True
    cdef cpp_bool out_syntactic_feasibility = True
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[cpp_pair[int, int]] new_match
    cdef cpp_unordered_set[int] current_match_G
    cdef cpp_unordered_set[int] current_match_H
    cdef cpp_unordered_map[int, int] forward_match
    cdef partial_maps_candidates_struct_directed candidates_struct

    # save if complete extension was reached
    if(current_match.size() == expected_order):
        all_matches.push_back(current_match)
        if(not all_extensions):
            return

    # if not optimal yet then obtain available pairs
    if(current_match.size() < expected_order):
        # generate auxiliary structures
        for each_pair in current_match:
            current_match_G.insert(each_pair.first)
            current_match_H.insert(each_pair.second)
            forward_match[each_pair.first] = each_pair.second

        # get candidate pairs
        candidates_struct = directed_candidates(expected_order,
                                                current_match,
                                                current_match_G,
                                                current_match_H,
                                                G.in_neighbors,
                                                H.in_neighbors,
                                                G.out_neighbors,
                                                H.out_neighbors,
                                                total_order)

        # evaluate candidates
        for each_pair in candidates_struct.candidates:
            # evaluate syntactic feasibility of in-neighbors
            in_syntactic_feasibility = syntactic_feasibility(each_pair.first,
                                                             each_pair.second,
                                                             candidates_struct.in_ring_G,
                                                             candidates_struct.in_ring_H,
                                                             current_match_G,
                                                             current_match_H,
                                                             forward_match,
                                                             G.in_neighbors,
                                                             H.in_neighbors)

            if(in_syntactic_feasibility):
                # evaluate syntactic feasibility of out-neighbors
                out_syntactic_feasibility = syntactic_feasibility(each_pair.first,
                                                                  each_pair.second,
                                                                  candidates_struct.out_ring_G,
                                                                  candidates_struct.out_ring_H,
                                                                  current_match_G,
                                                                  current_match_H,
                                                                  forward_match,
                                                                  G.out_neighbors,
                                                                  H.out_neighbors)

                if(out_syntactic_feasibility):
                    # evaluate semantic feasibility of in-neighors
                    if(node_labels or edge_labels):
                        in_semantic_feasibility = semantic_feasibility(node_labels,
                                                                       edge_labels,
                                                                       each_pair.first,
                                                                       each_pair.second,
                                                                       candidates_struct.in_ring_G,
                                                                       candidates_struct.in_ring_H,
                                                                       current_match_G,
                                                                       current_match_H,
                                                                       forward_match,
                                                                       G.nodes,
                                                                       H.nodes,
                                                                       G.in_neighbors,
                                                                       H.in_neighbors,
                                                                       G.edges,
                                                                       H.edges)

                    if(in_semantic_feasibility):
                        # evaluate semantic feasibility of out-neighors
                        if(node_labels or edge_labels):
                            out_semantic_feasibility = semantic_feasibility(node_labels,
                                                                            edge_labels,
                                                                            each_pair.first,
                                                                            each_pair.second,
                                                                            candidates_struct.out_ring_G,
                                                                            candidates_struct.out_ring_H,
                                                                            current_match_G,
                                                                            current_match_H,
                                                                            forward_match,
                                                                            G.nodes,
                                                                            H.nodes,
                                                                            G.out_neighbors,
                                                                            H.out_neighbors,
                                                                            G.edges,
                                                                            H.edges)

                        # push to stack if valid
                        if(out_semantic_feasibility):
                            # build new match
                            new_match.clear()
                            new_match = current_match
                            new_match.push_back(each_pair)
                            # extend match
                            directed_induced_connected_extensions_recursive(node_labels,
                                                                            edge_labels,
                                                                            all_extensions,
                                                                            expected_order,
                                                                            new_match,
                                                                            total_order,
                                                                            G,
                                                                            H,
                                                                            all_matches)

                            # finish if only one complete extension was requested and it was already found
                            if(not all_matches.empty()):
                                if(not all_extensions):
                                    return
    # end of function




# function: get candidate pairs for directed extension search ------------------
cdef partial_maps_candidates_struct_directed directed_candidates(size_t expected_order,
                                                                 cpp_vector[cpp_pair[int, int]] & current_match,
                                                                 cpp_unordered_set[int] & current_match_G,
                                                                 cpp_unordered_set[int] & current_match_H,
                                                                 cpp_unordered_map[int, cpp_unordered_set[int]] & in_neigh_G,
                                                                 cpp_unordered_map[int, cpp_unordered_set[int]] & in_neigh_H,
                                                                 cpp_unordered_map[int, cpp_unordered_set[int]] & out_neigh_G,
                                                                 cpp_unordered_map[int, cpp_unordered_set[int]] & out_neigh_H,
                                                                 cpp_unordered_map[int, int] & total_order) noexcept:

    # output holders
    cdef partial_maps_candidates_struct_directed candidates_struct

    # local variables
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int reference_minimum = 0
    cdef cpp_string comma
    comma.push_back(44)
    cdef cpp_string temp_string
    cdef cpp_pair[int, int] temp_pair
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[int] valid_G
    cdef cpp_vector[int] valid_H
    cdef cpp_vector[cpp_pair[int, int]] candidate_pairs
    cdef cpp_unordered_set[int] in_ring_G
    cdef cpp_unordered_set[int] in_ring_H
    cdef cpp_unordered_set[int] out_ring_G
    cdef cpp_unordered_set[int] out_ring_H
    cdef cpp_unordered_set[cpp_string] candidate_pairs_member

    # initialize reference minimum
    reference_minimum = 1 + 2*(<int>expected_order)

    # get candidate pairs from in-neighbors
    for each_pair in current_match:
        # reinitialize valid neighbors
        valid_G.clear()
        valid_H.clear()

        # get valid in-neighbors in G
        for node in in_neigh_G[each_pair.first]:
            # if not yet in match
            if(current_match_G.find(node) == current_match_G.end()):
                valid_G.push_back(node)

        # get valid in-neighbors in H
        for node in in_neigh_H[each_pair.second]:
            # if not yet in match
            if(current_match_H.find(node) == current_match_H.end()):
                valid_H.push_back(node)

        # make product of valid in-neighbors
        if((not valid_G.empty()) and (not valid_H.empty())):
            for node1 in valid_G:
                for node2 in valid_H:
                    # check if pair satisfies minimum total order condition
                    if(total_order[node2] <= reference_minimum):
                        temp_string = to_string(node1) + comma + to_string(node2)
                        if(candidate_pairs_member.find(temp_string) == candidate_pairs_member.end()):
                            temp_pair.first = node1
                            temp_pair.second = node2
                            if(total_order[node2] == reference_minimum):
                                # add proper pair
                                candidate_pairs.push_back(temp_pair)
                                # add string version for constant look ups
                                candidate_pairs_member.insert(temp_string)
                            if(total_order[node2] < reference_minimum):
                                # update control values
                                reference_minimum = total_order[node2]
                                candidate_pairs.clear()
                                candidate_pairs_member.clear()
                                # add proper pair
                                candidate_pairs.push_back(temp_pair)
                                # add string version for constant look ups
                                candidate_pairs_member.insert(temp_string)

        # additionally save neighbors in in-ring around match in G
        for node in valid_G:
            # save if unrepeated
            if(in_ring_G.find(node) == in_ring_G.end()):
                in_ring_G.insert(node)

        # additionally save neighbors in in-ring around match in H
        for node in valid_H:
            # save if unrepeated
            if(in_ring_H.find(node) == in_ring_H.end()):
                in_ring_H.insert(node)

        # additionally save valid out-neighbors in out-ring around match in G
        for node in out_neigh_G[each_pair.first]:
            # if not yet in match
            if(current_match_G.find(node) == current_match_G.end()):
                # save if unrepeated
                if(out_ring_G.find(node) == out_ring_G.end()):
                    out_ring_G.insert(node)

        # additionally save valid out-neighbors in out-ring around match in H
        for node in out_neigh_H[each_pair.second]:
            # if not yet in match
            if(current_match_H.find(node) == current_match_H.end()):
                # save if unrepeated
                if(out_ring_H.find(node) == out_ring_H.end()):
                    out_ring_H.insert(node)

    # alternatively get candidate pairs from out-neighbors
    if(candidate_pairs.empty()):
        # get candidate pairs from out-neighbors
        for each_pair in current_match:
            # reinitialize valid neighbors
            valid_G.clear()
            valid_H.clear()

            # get valid out-neighbors in G
            for node in out_neigh_G[each_pair.first]:
                # if not yet in match
                if(current_match_G.find(node) == current_match_G.end()):
                    valid_G.push_back(node)

            # get valid out-neighbors in H
            for node in out_neigh_H[each_pair.second]:
                # if not yet in match
                if(current_match_H.find(node) == current_match_H.end()):
                    valid_H.push_back(node)

            # make product of valid out-neighbors
            if((not valid_G.empty()) and (not valid_H.empty())):
                for node1 in valid_G:
                    for node2 in valid_H:
                        # check if pair satisfies minimum total order condition
                        if(total_order[node2] <= reference_minimum):
                            temp_string = to_string(node1) + comma + to_string(node2)
                            if(candidate_pairs_member.find(temp_string) == candidate_pairs_member.end()):
                                temp_pair.first = node1
                                temp_pair.second = node2
                                if(total_order[node2] == reference_minimum):
                                    # add proper pair
                                    candidate_pairs.push_back(temp_pair)
                                    # add string version for constant look ups
                                    candidate_pairs_member.insert(temp_string)
                                if(total_order[node2] < reference_minimum):
                                    # update control values
                                    reference_minimum = total_order[node2]
                                    candidate_pairs.clear()
                                    candidate_pairs_member.clear()
                                    # add proper pair
                                    candidate_pairs.push_back(temp_pair)
                                    # add string version for constant look ups
                                    candidate_pairs_member.insert(temp_string)

    # pack return structure
    candidates_struct.candidates = candidate_pairs
    candidates_struct.in_ring_G = in_ring_G
    candidates_struct.in_ring_H = in_ring_H
    candidates_struct.out_ring_G = out_ring_G
    candidates_struct.out_ring_H = out_ring_H
    # end of function
    return(candidates_struct)




# functions - feasability of matches - undirected and directed #################




# function: evaluate syntactic feasability for extension search ----------------
cdef cpp_bool syntactic_feasibility(int node1,
                                    int node2,
                                    cpp_unordered_set[int] & ring_G,
                                    cpp_unordered_set[int] & ring_H,
                                    cpp_unordered_set[int] & current_match_G,
                                    cpp_unordered_set[int] & current_match_H,
                                    cpp_unordered_map[int, int] & forward_match,
                                    cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_G,
                                    cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_H) noexcept:

    # local variables
    cdef int node = 0
    cdef int mapped = 0
    cdef cpp_unordered_set[int] neighbors_ring_G
    cdef cpp_unordered_set[int] neighbors_ring_H
    cdef cpp_unordered_set[int] neighbors_match_G
    cdef cpp_unordered_set[int] neighbors_match_H
    cdef cpp_unordered_set[int] neighbors_extern_G
    cdef cpp_unordered_set[int] neighbors_extern_H

    # loop-consistency-test
    if(neigh_G[node1].find(node1) != neigh_G[node1].end()):
        if(neigh_H[node2].find(node2) == neigh_H[node2].end()):
            # node1 has a loop in G but node2 has no loop in H
            return(False)
    if(neigh_G[node1].find(node1) == neigh_G[node1].end()):
        if(neigh_H[node2].find(node2) != neigh_H[node2].end()):
            # node1 has no loop in G but node2 has a loop in H
            return(False)

    # obtain tripartition of neighbors
    for node in neigh_G[node1]:
        if(current_match_G.find(node) != current_match_G.end()):
            # save map of neighbor to be compared later
            mapped = forward_match[node]
            neighbors_match_G.insert(mapped)
        else:
            if(ring_G.find(node) != ring_G.end()):
                # save neighbor since we are just comparing numbers later
                neighbors_ring_G.insert(node)
            else:
                # save neighbor since we are just comparing numbers later
                neighbors_extern_G.insert(node)
    for node in neigh_H[node2]:
        if(current_match_H.find(node) != current_match_H.end()):
            # save node of H to be compared to map of nodes from G
            neighbors_match_H.insert(node)
        else:
            if(ring_H.find(node) != ring_H.end()):
                # save neighbor since we are just comparing numbers later
                neighbors_ring_H.insert(node)
            else:
                # save neighbor since we are just comparing numbers later
                neighbors_extern_H.insert(node)

    # look ahead 0: consistency of neighbors in match
    if(neighbors_match_G.size() != neighbors_match_H.size()):
        # one node has more neighbors in the match than the other
        return(False)
    else:
        for mapped in neighbors_match_G:
            if(neighbors_match_H.find(mapped) == neighbors_match_H.end()):
                # the neighbors dont respect the match
                return(False)

    # look ahead 1: consistency of neighbors in ring (not in match but adjacent to match)
    if(neighbors_ring_G.size() != neighbors_ring_H.size()):
        return(False)

    # look ahead 2: consistency of extern neighbors (neither in match nor adjacent to match)
    if(neighbors_extern_G.size() != neighbors_extern_H.size()):
        return(False)

    # end of function
    return(True)




# function: evaluate semantic feasability for extension search -----------------
cdef cpp_bool semantic_feasibility(cpp_bool node_labels,
                                   cpp_bool edge_labels,
                                   int node1,
                                   int node2,
                                   cpp_unordered_set[int] & ring_G,
                                   cpp_unordered_set[int] & ring_H,
                                   cpp_unordered_set[int] & current_match_G,
                                   cpp_unordered_set[int] & current_match_H,
                                   cpp_unordered_map[int, int] & forward_match,
                                   cpp_unordered_map[int, int] & nodes_G,
                                   cpp_unordered_map[int, int] & nodes_H,
                                   cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_G,
                                   cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_H,
                                   cpp_unordered_map[cpp_string, int] & edges_G,
                                   cpp_unordered_map[cpp_string, int] & edges_H) noexcept:

    # local variables
    cdef int node = 0
    cdef int mapped = 0
    cdef int node_label = 0
    cdef int edge_label = 0
    cdef cpp_string comma
    comma.push_back(44)
    cdef cpp_string labeled_edge_G
    cdef cpp_string labeled_edge_H
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[int] neighbors_ring_G
    cdef cpp_vector[int] neighbors_ring_H
    cdef cpp_vector[int] neighbors_match_G
    cdef cpp_vector[int] neighbors_match_H
    cdef cpp_vector[int] neighbors_extern_G
    cdef cpp_vector[int] neighbors_extern_H
    cdef cpp_unordered_map[int, int] count_node_ring_G
    cdef cpp_unordered_map[int, int] count_node_ring_H
    cdef cpp_unordered_map[int, int] count_edge_ring_G
    cdef cpp_unordered_map[int, int] count_edge_ring_H
    cdef cpp_unordered_map[int, int] count_node_extern_G
    cdef cpp_unordered_map[int, int] count_node_extern_H
    cdef cpp_unordered_map[int, int] count_edge_extern_G
    cdef cpp_unordered_map[int, int] count_edge_extern_H

    if(node_labels):
        # compare vertex-labels
        if(nodes_G[node1] != nodes_H[node2]):
            return(False)

    if(edge_labels):
        # compare loop-labels
        if(neigh_G[node1].find(node1) != neigh_G[node1].end()):
            labeled_edge_G = to_string(node1) + comma + to_string(node1)
            labeled_edge_H = to_string(node2) + comma + to_string(node2)
            # compare labeled edges
            if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
                return(False)

    # obtain tripartition of neighbors in G
    for node in neigh_G[node1]:
        if(current_match_G.find(node) != current_match_G.end()):
            neighbors_match_G.push_back(node)
        else:
            if(ring_G.find(node) != ring_G.end()):
                # save neighbor since we are just comparing numbers later
                neighbors_ring_G.push_back(node)
            else:
                # save neighbor since we are just comparing numbers later
                neighbors_extern_G.push_back(node)

    if(edge_labels):
        # label look ahead 0: compare non-loop edge-labels in match
        for node in neighbors_match_G:
            # edge in G with only one end in match
            labeled_edge_G = to_string(node1) + comma + to_string(node)
            # edge in H with only one end in match
            labeled_edge_H = to_string(node2) + comma + to_string(forward_match[node])
            # compare labeled edges
            if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
                return(False)

    # obtain tripartition of neighbors in H
    for node in neigh_H[node2]:
        if(current_match_H.find(node) != current_match_H.end()):
            neighbors_match_H.push_back(node)
        else:
            if(ring_H.find(node) != ring_H.end()):
                neighbors_ring_H.push_back(node)
            else:
                neighbors_extern_H.push_back(node)

    # label look ahead 1: compare labels of neighbors in ring (not in match but adjacent to match)
    if(not neighbors_ring_G.empty()):
        # count in G
        for node in neighbors_ring_G:
            if(node_labels):
                # get node label
                node_label = nodes_G[node]
                # count node label
                if(count_node_ring_G.find(node_label) != count_node_ring_G.end()):
                    count_node_ring_G[node_label] = count_node_ring_G[node_label] + 1
                else:
                    count_node_ring_G[node_label] = 1
            if(edge_labels):
                # get edge label (possibly loop label)
                labeled_edge_G = to_string(node1) + comma + to_string(node)
                edge_label = edges_G[labeled_edge_G]
                # count edge label
                if(count_edge_ring_G.find(edge_label) != count_edge_ring_G.end()):
                    count_edge_ring_G[edge_label] = count_edge_ring_G[edge_label] + 1
                else:
                    count_edge_ring_G[edge_label] = 1

        # count in H
        for node in neighbors_ring_H:
            if(node_labels):
                # get node label
                node_label = nodes_H[node]
                # count node label
                if(count_node_ring_H.find(node_label) != count_node_ring_H.end()):
                    count_node_ring_H[node_label] = count_node_ring_H[node_label] + 1
                else:
                    count_node_ring_H[node_label] = 1
            if(edge_labels):
                # get edge label (possibly loop label)
                labeled_edge_H = to_string(node2) + comma + to_string(node)
                edge_label = edges_H[labeled_edge_H]
                # count edge label
                if(count_edge_ring_H.find(edge_label) != count_edge_ring_H.end()):
                    count_edge_ring_H[edge_label] = count_edge_ring_H[edge_label] + 1
                else:
                    count_edge_ring_H[edge_label] = 1

        # compare number of types of adjacent nodes
        if(node_labels):
            if(count_node_ring_G.size() != count_node_ring_H.size()):
                return(False)

        # compare number of types of incident edges
        if(edge_labels):
            if(count_edge_ring_G.size() != count_edge_ring_H.size()):
                return(False)

        # compare types of adjacent nodes
        if(node_labels):
            for each_pair in count_node_ring_G:
                if(count_node_ring_H.find(each_pair.first) == count_node_ring_H.end()):
                    return(False)
                else:
                    if(each_pair.second != count_node_ring_H[each_pair.first]):
                        return(False)

        # compare types of incident edges
        if(edge_labels):
            for each_pair in count_edge_ring_G:
                if(count_edge_ring_H.find(each_pair.first) == count_edge_ring_H.end()):
                    return(False)
                else:
                    if(each_pair.second != count_edge_ring_H[each_pair.first]):
                        return(False)

    # label look ahead 2: compare labels of extern neighbors (neither in match nor adjacent to match)
    if(not neighbors_extern_G.empty()):
        for node in neighbors_extern_G:
            if(node_labels):
                # get node label
                node_label = nodes_G[node]
                # count node label
                if(count_node_extern_G.find(node_label) != count_node_extern_G.end()):
                    count_node_extern_G[node_label] = count_node_extern_G[node_label] + 1
                else:
                    count_node_extern_G[node_label] = 1
            if(edge_labels):
                # get edge label (possibly loop label)
                labeled_edge_G = to_string(node1) + comma + to_string(node)
                edge_label = edges_G[labeled_edge_G]
                # count edge label
                if(count_edge_extern_G.find(edge_label) != count_edge_extern_G.end()):
                    count_edge_extern_G[edge_label] = count_edge_extern_G[edge_label] + 1
                else:
                    count_edge_extern_G[edge_label] = 1

        for node in neighbors_extern_H:
            if(node_labels):
                # get node label
                node_label = nodes_H[node]
                # count node label
                if(count_node_extern_H.find(node_label) != count_node_extern_H.end()):
                    count_node_extern_H[node_label] = count_node_extern_H[node_label] + 1
                else:
                    count_node_extern_H[node_label] = 1
            if(edge_labels):
                # get edge label (possibly loop label)
                labeled_edge_H = to_string(node2) + comma + to_string(node)
                edge_label = edges_H[labeled_edge_H]
                # count edge label
                if(count_edge_extern_H.find(edge_label) != count_edge_extern_H.end()):
                    count_edge_extern_H[edge_label] = count_edge_extern_H[edge_label] + 1
                else:
                    count_edge_extern_H[edge_label] = 1

        # compare number of types of adjacent nodes
        if(node_labels):
            if(count_node_extern_G.size() != count_node_extern_H.size()):
                return(False)

        # compare number of types of incident edges
        if(edge_labels):
            if(count_edge_extern_G.size() != count_edge_extern_H.size()):
                return(False)

        # compare types of adjacent nodes
        if(node_labels):
            for each_pair in count_node_extern_G:
                if(count_node_extern_H.find(each_pair.first) == count_node_extern_H.end()):
                    return(False)
                else:
                    if(each_pair.second != count_node_extern_H[each_pair.first]):
                        return(False)

        # compare types of incident edges
        if(edge_labels):
            for each_pair in count_edge_extern_G:
                if(count_edge_extern_H.find(each_pair.first) == count_edge_extern_H.end()):
                    return(False)
                else:
                    if(each_pair.second != count_edge_extern_H[each_pair.first]):
                        return(False)

    # end of function
    return(True)




################################################################################
################################################################################

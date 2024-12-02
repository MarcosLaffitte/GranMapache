################################################################################
#                                                                              #
# - GranMapache: GRAphs-and-Networks MAPping Applications                      #
#   with Cython and HEuristics                                                 #
#                                                                              #
# - Module: common_subgraphs                                                   #
#                                                                              #
# - Description: methods for obtaining maximum common subgraphs of different   #
#   kinds between graphs and study their properties.                           #
#                                                                              #
################################################################################




# C++ string encoding ##########################################################
# cython: c_string_type=unicode, c_string_encoding=utf8




# dependencies #################################################################




# already in python ------------------------------------------------------------
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
cdef struct common_subgraphs_undirected_graph:
    cpp_unordered_map[int, int] nodes
    cpp_unordered_map[cpp_string, int] edges
    cpp_unordered_map[int, cpp_unordered_set[int]] neighbors




# struct: directed graph -------------------------------------------------------
# using unordered maps and unordered sets for scalability
cdef struct common_subgraphs_directed_graph:
    cpp_unordered_map[int, int] nodes
    cpp_unordered_map[cpp_string, int] edges
    cpp_unordered_map[int, cpp_unordered_set[int]] in_neighbors
    cpp_unordered_map[int, cpp_unordered_set[int]] out_neighbors




# algorithms ###################################################################




# functions - maximum common induced subgraphs - wrapper #######################




# function: callable wrapper for the MCIS search -------------------------------
def maximum_common_induced_subgraphs(nx_G = nx.Graph(),          # can be nx.DiGraph
                                     nx_H = nx.Graph(),          # can be nx.DiGraph
                                     input_anchor = [],          # should be non-empty list
                                     node_labels = True,         # consider node labels when evaluating the common subgraphs
                                     edge_labels = True,         # consider edge labels when evaluating the common subgraphs
                                     stop_if_covering = True,    # stops when finding a maximum common subgraph that covers one of the input graphs
                                     iterative_search = False):  # by default a recursive search is used, otherwise an iterative version is called
    # description
    """
    > description: receives two networkx graphs G and H, and a possibly-empty match
    between them (here called anchor), and uses a VF2-like approach to obtain all maximum
    common induced subgraphs that extend or cover the anchor. If the anchor is empty then
    these are all the general maximum common induced subgraphs between G and H.

    > input:
    * nx_G - first networkx (di)graph being matched.
    * nx_H - second networkx (di)graph being matched.
    * input_anchor - inyective map as a non-empty list of 2-tuples (x, y) of nodes x
    from G and y from H.
    * node_labels - boolean indicating if node labels should be considered for the search,
    which is the default behavior, or if they should be ignored.
    * edge_labels - boolean indicating if edge labels should be considered for the search,
    which is the default behavior, or if they should be ignored.
    * stop_if_covering - boolean indicating if the function should stop as soon as finding
    one common induced subgraphs extending the anchor and completely covering all the vertices
    of one of the input graphs - this is the default behavior. It should be noticed that such
    common subgraph will be a subgraph isomorphism if and only if the input anchor already
    induces a common subgraph, or if it is empty. Setting this to False might result into a
    more time consuming call to this function.
    * iterative_search - boolean indicating if the iterative version of this algorithm should
    be used instead of the recursive version (default).

    > output:
    * mci_subgraphs - list of injective maps each as a list of 2-tuples (x, y) of nodes x
    from G and y from H representing each a different maximum common induced subgraph containing
    the anchor. If the anchor is empty then these are general common induced subgraphs.
    * covering - boolean value indicating if the common subgraph that is being returned covers
    all the vertices of the smallest of the input graphs.

    > calls:
    * .integerization.encode_graphs
    * .integerization.decode_graphs
    * .integerization.encode_match
    * .integerization.decode_match
    * undirected_maximum_common_induced_subgraphs_iterative
    * undirected_maximum_common_induced_subgraphs_recursive
    * directed_maximum_common_induced_subgraphs_iterative
    * directed_maximum_common_induced_subgraphs_recursive
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
    # check that input graphs are not null graphs
    if((nx_G.order() == 0) or (nx_H.order() == 0)):
        raise(ValueError("gmapache: input graphs must have at least one node each."))
    # check that the input graphs have the same type
    if((nx.is_directed(nx_G)) and (not nx.is_directed(nx_H))):
        raise(ValueError("gmapache: input graphs must be both directed or both undirected."))
    if((not nx.is_directed(nx_G)) and (nx.is_directed(nx_H))):
        raise(ValueError("gmapache: input graphs must be both directed or both undirected."))
    # check that fourth argument is a boolean variable
    if(type(node_labels) not in [type(test_bool)]):
        raise(ValueError("gmapache: fourth argument must be a boolean variable."))
    # check that fifth argument is a boolean variable
    if(type(edge_labels) not in [type(test_bool)]):
        raise(ValueError("gmapache: fifth argument must be a boolean variable."))
    # check that sixth argument is a boolean variable
    if(type(stop_if_covering) not in [type(test_bool)]):
        raise(ValueError("gmapache: sixth argument must be a boolean variable."))
    # check that seventh argument is a boolean variable
    if(type(iterative_search) not in [type(test_bool)]):
        raise(ValueError("gmapache: seventh argument must be a boolean variable."))
    # check that third argument is a list
    if(not type(input_anchor) in [type(test_list)]):
        raise(ValueError("gmapache: third argument must be a list of 2-tuples."))
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
    if(len(input_anchor) > 0):
        if(not len(list(set([x for (x, y) in input_anchor]))) == len(input_anchor)):
            raise(ValueError("gmapache: the input list must be an injective map and without repeated elements."))
        if(not len(list(set([y for (x, y) in input_anchor]))) == len(input_anchor)):
            raise(ValueError("gmapache: the input list must be an injective map and without repeated elements."))

    # output holders
    cdef list mci_subgraphs = []
    covering = False

    # local variables (cython)
    cdef int node = 0
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
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[cpp_pair[int, int]] encoded_anchor
    cdef cpp_vector[cpp_pair[int, int]] each_subgraph
    cdef cpp_vector[cpp_vector[cpp_pair[int, int]]] encoded_subgraphs
    cdef cpp_unordered_map[int, int] total_order
    cdef partial_maps_directed_graph directed_G
    cdef partial_maps_directed_graph directed_H
    cdef partial_maps_undirected_graph undirected_G
    cdef partial_maps_undirected_graph undirected_H

    # local variables (python)
    cdef set node_labels_G = set()
    cdef set node_labels_H = set()
    cdef list node_labels_intersection = set()
    cdef list encoded_graphs = []
    cdef dict info = dict()
    cdef dict encoded_node_names = dict()
    cdef dict encoded_node_label = dict()
    cdef dict encoded_edge_label = dict()

    # encode graphs
    encoded_graphs, encoded_node_names, encoded_node_label, encoded_edge_label = encode_graphs([nx_G, nx_H])

    # encode match
    encoded_anchor = encode_match(input_anchor, encoded_node_names)

    # check that the graphs share at least one vertex with the same label (only if the search should preserve node labels)
    if(node_labels):
        node_labels_G = set([info["GMNL"] for (node, info) in encoded_graphs[0].nodes(data = True)])
        node_labels_H = set([info["GMNL"] for (node, info) in encoded_graphs[1].nodes(data = True)])
        node_labels_intersection = list(node_labels_G.intersection(node_labels_H))
        if(len(node_labels_intersection) == 0):
            raise(ValueError("gmapache: requested preservation of node labels but the input graphs share no nodes with common label."))

    # prepare nodes
    if(nx.is_directed(nx_G)):
        directed_G.nodes = {node:info["GMNL"] for (node, info) in encoded_graphs[0].nodes(data = True)}
        directed_H.nodes = {node:info["GMNL"] for (node, info) in encoded_graphs[1].nodes(data = True)}
    else:
        undirected_G.nodes = {node:info["GMNL"] for (node, info) in encoded_graphs[0].nodes(data = True)}
        undirected_H.nodes = {node:info["GMNL"] for (node, info) in encoded_graphs[1].nodes(data = True)}

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
    # NOTE: really doesn't depend on directedness of the graph, it's just selecting the structure that was initialized
    if(nx.is_directed(nx_H)):
        for each_pair in directed_H.nodes:
            node = each_pair.first
            counter = counter + 1
            total_order[node] = counter
    else:
        for each_pair in undirected_H.nodes:
            node = each_pair.first
            counter = counter + 1
            total_order[node] = counter

    # get expected order
    expected_order = min([nx_G.order(), nx_H.order()])

    # set recursion limit if recursive version was requested
    if(iterative_search):
        scalation_value = 1.5
        required_limit = max([nx_G.order(), nx_H.order()])
        current_limit = getrecursionlimit()
        if(current_limit < (scalation_value * required_limit)):
            setrecursionlimit(int(scalation_value * required_limit))

    # get maximum common induced subgraphs
    if(nx.is_directed(nx_G)):
        if(iterative_search):
            directed_maximum_common_induced_subgraphs_iterative(node_labels,
                                                                edge_labels,
                                                                stop_if_covering,
                                                                expected_order,
                                                                encoded_anchor,
                                                                total_order,
                                                                directed_G,
                                                                directed_H,
                                                                encoded_subgraphs)
        else:
            directed_maximum_common_induced_subgraphs_recursive(node_labels,
                                                                edge_labels,
                                                                stop_if_covering,
                                                                expected_order,
                                                                encoded_anchor,
                                                                total_order,
                                                                directed_G,
                                                                directed_H,
                                                                encoded_subgraphs)
    else:
        if(iterative_search):
            undirected_maximum_common_induced_subgraphs_iterative(node_labels,
                                                                  edge_labels,
                                                                  stop_if_covering,
                                                                  expected_order,
                                                                  encoded_anchor,
                                                                  total_order,
                                                                  undirected_G,
                                                                  undirected_H,
                                                                  encoded_subgraphs)
        else:
            undirected_maximum_common_induced_subgraphs_recursive(node_labels,
                                                                  edge_labels,
                                                                  stop_if_covering,
                                                                  expected_order,
                                                                  encoded_anchor,
                                                                  total_order,
                                                                  undirected_G,
                                                                  undirected_H,
                                                                  encoded_subgraphs)

    # decode maximum common induced subgraphs
    for each_subgraph in encoded_subgraphs:
        mci_subgraphs.append(decode_match(list(each_subgraph), encoded_node_names))

    # check if the returned subgraph covers all the vertices of the smallest graph
    # NOTE: at this point there should be at least one common subgraph with at least one node
    if((len(mci_subgraphs[0]) == nx_G.order()) or (len(mci_subgraphs[0]) == nx_H.order())):
        covering = True

    # end of function
    return(mci_subgraphs, covering)




# functions - maximum common induced subgraphs - undirected ####################




# function: core routine of VF2-like undirected approach - iterative -----------
# NOTE: an iterative DFS version of this algorithm can be implemented without a "visited"
# list, since the total order given to the vertices of the second graph guarantees that
# the search space is actually a search tree, and thus cannot have repeated states.
cdef void undirected_maximum_common_induced_subgraphs_iterative(cpp_bool node_labels,
                                                                cpp_bool edge_labels,
                                                                cpp_bool stop_if_covering,
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
    cdef cpp_vector[cpp_pair[int, int]] candidates
    cdef cpp_vector[cpp_pair[int, int]] current_match
    cdef cpp_unordered_set[int] current_match_G
    cdef cpp_unordered_set[int] current_match_H
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
                if(current_match.size() == expected_order):
                    if(stop_if_covering):
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
                if(current_match.size() == expected_order):
                    if(stop_if_covering):
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
                current_match_G.insert(each_pair.first)
                current_match_H.insert(each_pair.second)
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
                syntactic_feasibility_res = syntactic_feasibility(each_pair.first,
                                                                  each_pair.second,
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
                                                                        current_match_G,
                                                                        forward_match,
                                                                        G.nodes,
                                                                        H.nodes,
                                                                        G.neighbors,
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
cdef void undirected_maximum_common_induced_subgraphs_recursive(cpp_bool node_labels,
                                                                cpp_bool edge_labels,
                                                                cpp_bool stop_if_covering,
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
    cdef cpp_vector[cpp_pair[int, int]] candidates
    cdef cpp_unordered_set[int] current_match_G
    cdef cpp_unordered_set[int] current_match_H
    cdef cpp_unordered_map[int, int] forward_match

    # test initial match and improve if possible
    if(all_matches.empty()):
        if(not current_match.empty()):
            all_matches.push_back(current_match)
            # if complete match and only one then return
            if(current_match.size() == expected_order):
                if(stop_if_covering):
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
            if(current_match.size() == expected_order):
                if(stop_if_covering):
                    return

    # if not optimal yet then obtain available pairs
    if(current_match.size() < expected_order):
        # generate auxiliary structures
        for each_pair in current_match:
            current_match_G.insert(each_pair.first)
            current_match_H.insert(each_pair.second)
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
            syntactic_feasibility_res = syntactic_feasibility(each_pair.first,
                                                              each_pair.second,
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
                                                                    current_match_G,
                                                                    forward_match,
                                                                    G.nodes,
                                                                    H.nodes,
                                                                    G.neighbors,
                                                                    G.edges,
                                                                    H.edges)

                if(semantic_feasibility_res):
                    # build new match
                    new_match.clear()
                    new_match = current_match
                    new_match.push_back(each_pair)
                    # extend match
                    undirected_maximum_common_induced_subgraphs_recursive(node_labels,
                                                                          edge_labels,
                                                                          stop_if_covering,
                                                                          expected_order,
                                                                          new_match,
                                                                          total_order,
                                                                          G,
                                                                          H,
                                                                          all_matches)

                    # finish if stop when covering was requested
                    if(not all_matches.empty()):
                        if(all_matches[0].size() == expected_order):
                            if(stop_if_covering):
                                return
    # end of function




# function: get candidate pairs for undirected extension search ----------------
cdef cpp_vector[cpp_pair[int, int]] undirected_candidates(cpp_vector[cpp_pair[int, int]] & current_match,
                                                          cpp_unordered_set[int] & current_match_G,
                                                          cpp_unordered_set[int] & current_match_H,
                                                          cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_G,
                                                          cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_H,
                                                          cpp_unordered_map[int, int] & total_order) noexcept:

    # output holders
    cdef cpp_vector[cpp_pair[int, int]] candidate_pairs

    # local variables (cython)
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int reference_maximum = 0
    cdef cpp_string comma
    comma.push_back(44)
    cdef cpp_string temp_string
    cdef cpp_pair[int, int] temp_pair
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[int] valid_G
    cdef cpp_vector[int] valid_H
    cdef cpp_unordered_set[cpp_string] candidate_pairs_member

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
            if(current_match_G.find(node) == current_match_G.end()):
                valid_G.push_back(node)

        # get valid neighbors in H
        for node in neigh_H[each_pair.second]:
            # if total order greater than in match
            if(total_order[node] > reference_maximum):
                # if not yet in match
                if(current_match_H.find(node) == current_match_H.end()):
                    valid_H.push_back(node)

        # make product of valid neighbors
        if((not valid_G.empty()) and (not valid_H.empty())):
            for node1 in valid_G:
                for node2 in valid_H:
                    temp_string = to_string(node1) + comma + to_string(node2)
                    if(candidate_pairs_member.find(temp_string) == candidate_pairs_member.end()):
                        # add proper pair
                        temp_pair.first = node1
                        temp_pair.second = node2
                        candidate_pairs.push_back(temp_pair)
                        # add string version for constant look ups
                        candidate_pairs_member.insert(temp_string)

    # end of function
    return(candidate_pairs)




# functions - maximum connected extensions - directed ##########################




# function: core routine of VF2-like directed approach - iterative -------------
# NOTE: an iterative DFS version of this algorithm can be implemented without a "visited"
# list, since the total order given to the vertices of the second graph guarantees that
# the search space is actually a search tree, and thus cannot have repeated states.
cdef void directed_maximum_common_induced_subgraphs_iterative(cpp_bool node_labels,
                                                              cpp_bool edge_labels,
                                                              cpp_bool stop_if_covering,
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
    cdef cpp_vector[cpp_pair[int, int]] candidates
    cdef cpp_vector[cpp_pair[int, int]] current_match
    cdef cpp_unordered_set[int] current_match_G
    cdef cpp_unordered_set[int] current_match_H
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
                if(current_match.size() == expected_order):
                    if(stop_if_covering):
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
                if(current_match.size() == expected_order):
                    if(stop_if_covering):
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
                current_match_G.insert(each_pair.first)
                current_match_H.insert(each_pair.second)
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
                # evaluate syntactic feasibility with in-neighbors
                in_syntactic_feasibility = syntactic_feasibility(each_pair.first,
                                                                 each_pair.second,
                                                                 current_match_G,
                                                                 current_match_H,
                                                                 forward_match,
                                                                 G.in_neighbors,
                                                                 H.in_neighbors)

                if(in_syntactic_feasibility):
                    # evaluate syntactic feasibility with out-neighbors
                    out_syntactic_feasibility = syntactic_feasibility(each_pair.first,
                                                                      each_pair.second,
                                                                      current_match_G,
                                                                      current_match_H,
                                                                      forward_match,
                                                                      G.out_neighbors,
                                                                      H.out_neighbors)

                    if(out_syntactic_feasibility):
                        # evaluate semantic feasibility with in-neighbors
                        if(node_labels or edge_labels):
                            in_semantic_feasibility = semantic_feasibility(node_labels,
                                                                           edge_labels,
                                                                           each_pair.first,
                                                                           each_pair.second,
                                                                           current_match_G,
                                                                           forward_match,
                                                                           G.nodes,
                                                                           H.nodes,
                                                                           G.in_neighbors,
                                                                           G.edges,
                                                                           H.edges)

                        if(in_semantic_feasibility):
                            # evaluate semantic feasibility with out-neighbors
                            if(node_labels or edge_labels):
                                out_semantic_feasibility = semantic_feasibility(node_labels,
                                                                                edge_labels,
                                                                                each_pair.first,
                                                                                each_pair.second,
                                                                                current_match_G,
                                                                                forward_match,
                                                                                G.nodes,
                                                                                H.nodes,
                                                                                G.out_neighbors,
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
cdef void directed_maximum_common_induced_subgraphs_recursive(cpp_bool node_labels,
                                                              cpp_bool edge_labels,
                                                              cpp_bool stop_if_covering,
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
    cdef cpp_vector[cpp_pair[int, int]] candidates
    cdef cpp_unordered_set[int] current_match_G
    cdef cpp_unordered_set[int] current_match_H
    cdef cpp_unordered_map[int, int] forward_match

    # test initial match and consecutive matches
    if(all_matches.empty()):
        if(not current_match.empty()):
            all_matches.push_back(current_match)
            # if complete match and only one then return
            if(current_match.size() == expected_order):
                if(stop_if_covering):
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
            if(current_match.size() == expected_order):
                if(stop_if_covering):
                    return

    # if not optimal yet then obtain available pairs
    if(current_match.size() < expected_order):
        # generate auxiliary structures
        for each_pair in current_match:
            current_match_G.insert(each_pair.first)
            current_match_H.insert(each_pair.second)
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
            # evaluate syntactic feasibility with in-neighbors
            in_syntactic_feasibility = syntactic_feasibility(each_pair.first,
                                                             each_pair.second,
                                                             current_match_G,
                                                             current_match_H,
                                                             forward_match,
                                                             G.in_neighbors,
                                                             H.in_neighbors)

            if(in_syntactic_feasibility):
                # evaluate syntactic feasibility with out-neighbors
                out_syntactic_feasibility = syntactic_feasibility(each_pair.first,
                                                                  each_pair.second,
                                                                  current_match_G,
                                                                  current_match_H,
                                                                  forward_match,
                                                                  G.out_neighbors,
                                                                  H.out_neighbors)

                if(out_syntactic_feasibility):
                    # evaluate semantic feasibility with in-neighbors
                    if(node_labels or edge_labels):
                        in_semantic_feasibility = semantic_feasibility(node_labels,
                                                                       edge_labels,
                                                                       each_pair.first,
                                                                       each_pair.second,
                                                                       current_match_G,
                                                                       forward_match,
                                                                       G.nodes,
                                                                       H.nodes,
                                                                       G.in_neighbors,
                                                                       G.edges,
                                                                       H.edges)

                    if(in_semantic_feasibility):
                        # evaluate semantic feasibility with out-neighbors
                        if(node_labels or edge_labels):
                            out_semantic_feasibility = semantic_feasibility(node_labels,
                                                                            edge_labels,
                                                                            each_pair.first,
                                                                            each_pair.second,
                                                                            current_match_G,
                                                                            forward_match,
                                                                            G.nodes,
                                                                            H.nodes,
                                                                            G.out_neighbors,
                                                                            G.edges,
                                                                            H.edges)

                        # push to stack if valid
                        if(out_semantic_feasibility):
                            # build new match
                            new_match.clear()
                            new_match = current_match
                            new_match.push_back(each_pair)
                            # extend match
                            directed_maximum_common_induced_subgraphs_recursive(node_labels,
                                                                                edge_labels,
                                                                                stop_if_covering,
                                                                                expected_order,
                                                                                new_match,
                                                                                total_order,
                                                                                G,
                                                                                H,
                                                                                all_matches)
                            # finish if stop when covering was requested
                            if(not all_matches.empty()):
                                if(all_matches[0].size() == expected_order):
                                    if(stop_if_covering):
                                        return
    # end of function




# function: get candidate pairs for directed extension search ------------------
cdef cpp_vector[cpp_pair[int, int]] directed_candidates(cpp_vector[cpp_pair[int, int]] & current_match,
                                                        cpp_unordered_set[int] & current_match_G,
                                                        cpp_unordered_set[int] & current_match_H,
                                                        cpp_unordered_map[int, cpp_unordered_set[int]] & in_neigh_G,
                                                        cpp_unordered_map[int, cpp_unordered_set[int]] & in_neigh_H,
                                                        cpp_unordered_map[int, cpp_unordered_set[int]] & out_neigh_G,
                                                        cpp_unordered_map[int, cpp_unordered_set[int]] & out_neigh_H,
                                                        cpp_unordered_map[int, int] & total_order) noexcept:

    # output holders
    cdef cpp_vector[cpp_pair[int, int]] candidate_pairs

    # local variables (cython)
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int reference_maximum = 0
    cdef cpp_string comma
    comma.push_back(44)
    cdef cpp_string temp_string
    cdef cpp_pair[int, int] temp_pair
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[int] valid_G
    cdef cpp_vector[int] valid_H
    cdef cpp_unordered_set[cpp_string] candidate_pairs_member

    # get maximum value of total order in match
    for node in current_match_H:
        if(total_order[node] > reference_maximum):
            reference_maximum = total_order[node]

    # get candidates based on in-neighbors
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
            # if total order greater than in match
            if(total_order[node] > reference_maximum):
                # if not yet in match
                if(current_match_H.find(node) == current_match_H.end()):
                    valid_H.push_back(node)

        # make product of valid in-neighbors
        if((not valid_G.empty()) and (not valid_H.empty())):
            for node1 in valid_G:
                for node2 in valid_H:
                    temp_string = to_string(node1) + comma + to_string(node2)
                    if(candidate_pairs_member.find(temp_string) == candidate_pairs_member.end()):
                        # add proper pair
                        temp_pair.first = node1
                        temp_pair.second = node2
                        candidate_pairs.push_back(temp_pair)
                        # add string version for constant look ups
                        candidate_pairs_member.insert(temp_string)

    # alternatively get candidate pairs from out-neighbors
    if(candidate_pairs.empty()):
        # get candidates based on valid sets
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
                # if total order greater than in match
                if(total_order[node] > reference_maximum):
                    # if not yet in match
                    if(current_match_H.find(node) == current_match_H.end()):
                        valid_H.push_back(node)

            # make product of valid out-neighbors
            if((not valid_G.empty()) and (not valid_H.empty())):
                for node1 in valid_G:
                    for node2 in valid_H:
                        temp_string = to_string(node1) + comma + to_string(node2)
                        if(candidate_pairs_member.find(temp_string) == candidate_pairs_member.end()):
                            # add proper pair
                            temp_pair.first = node1
                            temp_pair.second = node2
                            candidate_pairs.push_back(temp_pair)
                            # add string version for constant look ups
                            candidate_pairs_member.insert(temp_string)

    # end of function
    return(candidate_pairs)




# functions - feasability of matches - undirected and directed #################




# function: evaluate syntactic feasability for connected extension -------------
cdef cpp_bool syntactic_feasibility(int node1,
                                    int node2,
                                    cpp_unordered_set[int] & current_match_G,
                                    cpp_unordered_set[int] & current_match_H,
                                    cpp_unordered_map[int, int] & forward_match,
                                    cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_G,
                                    cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_H) noexcept:

    # local variables (cython)
    cdef int node = 0
    cdef int mapped = 0
    cdef cpp_unordered_set[int] neighbors_match_G
    cdef cpp_unordered_set[int] neighbors_match_H

    # loop-consistency-test
    if(neigh_G[node1].find(node1) != neigh_G[node1].end()):
        if(neigh_H[node2].find(node2) == neigh_H[node2].end()):
            # node1 has a loop in G but node2 has no loop in H
            return(False)
    if(neigh_G[node1].find(node1) == neigh_G[node1].end()):
        if(neigh_H[node2].find(node2) != neigh_H[node2].end()):
            # node1 has no loop in G but node2 has a loop in H
            return(False)

    # look ahead 0: consistency of neighbors in match
    for node in neigh_G[node1]:
        if(current_match_G.find(node) != current_match_G.end()):
            mapped = forward_match[node]
            neighbors_match_G.insert(mapped)

    for node in neigh_H[node2]:
        if(current_match_H.find(node) != current_match_H.end()):
            neighbors_match_H.insert(node)

    if(neighbors_match_G.size() != neighbors_match_H.size()):
        # one node has more neighbors in the match than the other
        return(False)
    else:
        for mapped in neighbors_match_G:
            if(neighbors_match_H.find(mapped) == neighbors_match_H.end()):
                # the neighbors dont respect the match
                return(False)

    # end of function
    return(True)




# function: evaluate semantic feasability for connected extension --------------
cdef cpp_bool semantic_feasibility(cpp_bool node_labels,
                                   cpp_bool edge_labels,
                                   int node1,
                                   int node2,
                                   cpp_unordered_set[int] & current_match_G,
                                   cpp_unordered_map[int, int] & forward_match,
                                   cpp_unordered_map[int, int] & nodes_G,
                                   cpp_unordered_map[int, int] & nodes_H,
                                   cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_G,
                                   cpp_unordered_map[cpp_string, int] & edges_G,
                                   cpp_unordered_map[cpp_string, int] & edges_H) noexcept:

    # local variables (cython)
    cdef int node = 0
    cdef cpp_string comma
    comma.push_back(44)
    cdef cpp_string labeled_edge_G
    cdef cpp_string labeled_edge_H

    if(node_labels):
        # compare vertex-labels
        if(nodes_G[node1] != nodes_H[node2]):
            return(False)

    if(edge_labels):
        # compare loop-labels
        if(neigh_G[node1].find(node1) != neigh_G[node1].end()):
            # loop in G
            labeled_edge_G = to_string(node1) + comma + to_string(node1)
            # loop in H
            labeled_edge_H = to_string(node2) + comma + to_string(node2)
            # compare edge labels
            if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
                return(False)

        # compare non-loop edge-labels
        for node in neigh_G[node1]:
            if(current_match_G.find(node) != current_match_G.end()):
                # edge in G with only one end in match
                labeled_edge_G = to_string(node1) + comma + to_string(node)
                # edge in H with only one end in match
                labeled_edge_H = to_string(node2) + comma + to_string(forward_match[node])
                # compare edge labels
                if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
                    return(False)

    # end of function
    return(True)




################################################################################
################################################################################




# CONSTRUCTION ZONE ############################################################




# TRANSOFRMING THIS INTO ANCHORED_COMMON_CONNECTED_INDUCED_SUBGRAPH

# functions - maximum connected extensions - wrapper ###########################




# # function: callable wrapper for the maximum connected extensions --------------
# def maximum_connected_extensions(nx_G = nx.Graph(),          # can be nx.DiGraph
#                                  nx_H = nx.Graph(),          # can be nx.DiGraph
#                                  input_anchor = [],          # should be non-empty list
#                                  node_labels = True,         # consider node labels when evaluating the extensions
#                                  edge_labels = True,         # consider edge labels when evaluating the extensions
#                                  all_extensions = False,     # by default stops when finding one complete extension (if any)
#                                  iterative_search = True):   # by default an iterative search is used, otherwise a recursive version is called
#     # description
#     """
#     > description: receives two networkx graphs G and H, and a match between them (here
#     called anchor), and uses a VF2-like approach to obtain the maximum extensions of the
#     anchor producing connected common subgraphs (not necessarily maximum themselves). The
#     anchor alone also produces a subgraph, which may not be an induced common subgraph,
#     but the subgraph produced by any extension after removing the achor is always induced.
#     An additional boolean variable determines if the search should be done iteratively
#     (by default), or recursively.

#     > input:
#     * nx_G - first networkx (di)graph being matched.
#     * nx_H - second networkx (di)graph being matched.
#     * input_anchor - inyective map as a non-empty list of 2-tuples (x, y) of nodes x
#     from G and y from H. An exception is raised if the anchor is empty.
#     * node_labels - boolean indicating if node labels should be considered for the search,
#     which is the default behavior, or if they should be ignored.
#     * edge_labels - boolean indicating if edge labels should be considered for the search,
#     which is the default behavior, or if they should be ignored.
#     * all_extensions - boolean indicating if the function should stop as soon as one complete
#     extension is found (if any) - this is the default behavior- or if it should search for all
#     possible (complete) extensions. NOTE: mathematically speaking, the complete extension of
#     a FIXED and GOOD anchor is unique up to equivalence of bijections, that is, equivalence
#     of atom maps or correspondingly isomorphism of their ITS graphs. Nontheless it should be
#     noted that changing the anchor may produce a non-equivalent (complete) extension. In other
#     words, each call to this function can (mathematically) produce only one complete extension,
#     but calls with different anchors can produce non-equivalent extensions, even if the anchors
#     themselves produce isomorphic partial ITS graphs.
#     * iterative_search - boolean indicating if the iterative version of this algorithm should
#     be used (the default), or if a recursive version of it should be used instead.

#     > output:
#     * extensions - list of injective maps each as a list of 2-tuples (x, y) of nodes x
#     from G and y from H representing the maximum connected extensions of the anchor (each
#     extension contains the anchor as a sublist).
#     * good_anchor - boolean value indicating if the extensions are complete, i.e., cover
#     all nodes of G and of H, and thus if they are bijections between G and H. If so, the
#     anchor is what we have refered to as a "good partial atom map", and equivalenteĺy the
#     match obtained when removing the anchhor from any extension is a graph-isomorphism
#     between the "remainder" graphs it induces from G and H.

#     > calls:
#     * .integerization.encode_graphs
#     * .integerization.decode_graphs
#     * .integerization.encode_match
#     * .integerization.decode_match
#     * undirected_maximum_connected_extensions_iterative
#     * undirected_maximum_connected_extensions_recursive
#     * directed_maximum_connected_extensions_iterative
#     * directed_maximum_connected_extensions_recursive
#     """

#     # exception handling and input correctness
#     test_list = [0, 0]
#     test_tuple = (0, 0)
#     test_undir = nx.Graph()
#     test_dir = nx.DiGraph()
#     test_bool = False
#     # check that first argument is networkx graph or digraph
#     if(type(nx_G) not in [type(test_undir), type(test_dir)]):
#         raise(ValueError("gmapache: first argument must be a networkx graph or digraph."))
#     # check that second argument is networkx graph or digraph
#     if(type(nx_H) not in [type(test_undir), type(test_dir)]):
#         raise(ValueError("gmapache: second argument must be a networkx graph or digraph."))
#     # check that input graphs are not null graphs
#     if((nx_G.order() == 0) or (nx_H.order() == 0)):
#         raise(ValueError("gmapache: input graphs must have at least one node each."))
#     # check that the input graphs have the same type
#     if((nx.is_directed(nx_G)) and (not nx.is_directed(nx_H))):
#         raise(ValueError("gmapache: input graphs must be both directed or both undirected."))
#     if((not nx.is_directed(nx_G)) and (nx.is_directed(nx_H))):
#         raise(ValueError("gmapache: input graphs must be both directed or both undirected."))
#     # check that fourth argument is a boolean variable
#     if(type(node_labels) not in [type(test_bool)]):
#         raise(ValueError("gmapache: fourth argument must be a boolean variable."))
#     # check that fifth argument is a boolean variable
#     if(type(edge_labels) not in [type(test_bool)]):
#         raise(ValueError("gmapache: fifth argument must be a boolean variable."))
#     # check that sixth argument is a boolean variable
#     if(type(all_extensions) not in [type(test_bool)]):
#         raise(ValueError("gmapache: sixth argument must be a boolean variable."))
#     # check that seventh argument is a boolean variable
#     if(type(iterative_search) not in [type(test_bool)]):
#         raise(ValueError("gmapache: seventh argument must be a boolean variable."))
#     # check that third argument is a list
#     if(not type(input_anchor) in [type(test_list)]):
#         raise(ValueError("gmapache: third argument must be a non-empty list of 2-tuples."))
#     if(len(input_anchor) == 0):
#         raise(ValueError("gmapache: third argument must be a non-empty list of 2-tuples."))
#     # check correctness of entries in third argument
#     for test_entry in input_anchor:
#         if(not type(test_entry) in [type(test_tuple)]):
#             raise(ValueError("gmapache: all elements in input list must be tuples."))
#         if(not len(test_entry) == 2):
#             raise(ValueError("gmapache: all tuples in input list must be of lenght 2."))
#         if(test_entry[0] not in list(nx_G.nodes())):
#             raise(ValueError("gmapache: the input list is matching a vertex not present in the first graph."))
#         if(test_entry[1] not in list(nx_H.nodes())):
#             raise(ValueError("gmapache: the input list is matching a vertex not present in the second graph."))
#     # check amount of entries in third argument
#     if(not len(list(set([x for (x, y) in input_anchor]))) == len(input_anchor)):
#         raise(ValueError("gmapache: the input list must be an injective map and without repeated elements."))
#     if(not len(list(set([y for (x, y) in input_anchor]))) == len(input_anchor)):
#         raise(ValueError("gmapache: the input list must be an injective map and without repeated elements."))

#     # output holders
#     cdef list extensions = []
#     good_anchor = False

#     # local variables (cython)
#     cdef int node = 0
#     cdef int node1 = 0
#     cdef int node2 = 0
#     cdef int counter = 0
#     cdef int current_limit = 0
#     cdef int required_limit = 0
#     cdef float scalation_value = 0
#     cdef size_t expected_order = 0
#     cdef cpp_string comma
#     comma.push_back(44)
#     cdef cpp_string temp_str
#     cdef cpp_pair[int, cpp_unordered_set[int]] each_pair
#     cdef cpp_vector[int] next_level
#     cdef cpp_vector[int] current_level
#     cdef cpp_vector[cpp_pair[int, int]] encoded_anchor
#     cdef cpp_vector[cpp_pair[int, int]] each_extension
#     cdef cpp_vector[cpp_vector[cpp_pair[int, int]]] encoded_extensions
#     cdef cpp_unordered_map[int, int] total_order
#     cdef cpp_unordered_map[int, cpp_bool] visited
#     cdef cpp_unordered_map[int, cpp_unordered_set[int]] connectivity_neighbors
#     cdef partial_maps_directed_graph directed_G
#     cdef partial_maps_directed_graph directed_H
#     cdef partial_maps_undirected_graph undirected_G
#     cdef partial_maps_undirected_graph undirected_H

#     # local variables (python)
#     cdef list encoded_graphs = []
#     cdef dict info = dict()
#     cdef dict encoded_node_names = dict()
#     cdef dict encoded_node_label = dict()
#     cdef dict encoded_edge_label = dict()
#     undirected_copy_H = None

#     # encode graphs
#     encoded_graphs, encoded_node_names, encoded_node_label, encoded_edge_label = encode_graphs([nx_G, nx_H])

#     # encode match
#     encoded_anchor = encode_match(input_anchor, encoded_node_names)

#     # prepare nodes
#     if(nx.is_directed(nx_G)):
#         directed_G.nodes = {node:info["GMNL"] for (node, info) in encoded_graphs[0].nodes(data = True)}
#         directed_H.nodes = {node:info["GMNL"] for (node, info) in encoded_graphs[1].nodes(data = True)}
#     else:
#         undirected_G.nodes = {node:info["GMNL"] for (node, info) in encoded_graphs[0].nodes(data = True)}
#         undirected_H.nodes = {node:info["GMNL"] for (node, info) in encoded_graphs[1].nodes(data = True)}

#     # prepare edges
#     if(nx.is_directed(nx_G)):
#         directed_G.edges = {str(node1)+","+str(node2):info["GMEL"] for (node1, node2, info) in encoded_graphs[0].edges(data = True)}
#         directed_H.edges = {str(node1)+","+str(node2):info["GMEL"] for (node1, node2, info) in encoded_graphs[1].edges(data = True)}
#     else:
#         # prepare edges of undirected G
#         for (node1, node2, info) in encoded_graphs[0].edges(data = True):
#             if(node1 == node2):
#                 temp_str = to_string(node1) + comma + to_string(node1)
#                 undirected_G.edges[temp_str] = info["GMEL"]
#             else:
#                 # save the two label edges to simplify future access
#                 temp_str = to_string(node1) + comma + to_string(node2)
#                 undirected_G.edges[temp_str] = info["GMEL"]
#                 temp_str = to_string(node2) + comma + to_string(node1)
#                 undirected_G.edges[temp_str] = info["GMEL"]
#         # prepare edges of undirected H
#         for (node1, node2, info) in encoded_graphs[1].edges(data = True):
#             if(node1 == node2):
#                 temp_str = to_string(node1) + comma + to_string(node1)
#                 undirected_H.edges[temp_str] = info["GMEL"]
#             else:
#                 # save the two label edges to simplify future access
#                 temp_str = to_string(node1) + comma + to_string(node2)
#                 undirected_H.edges[temp_str] = info["GMEL"]
#                 temp_str = to_string(node2) + comma + to_string(node1)
#                 undirected_H.edges[temp_str] = info["GMEL"]

#     # prepare neighbors
#     if(nx.is_directed(nx_G)):
#         directed_G.in_neighbors = {node:set(encoded_graphs[0].predecessors(node)) for node in list(encoded_graphs[0].nodes())}
#         directed_H.in_neighbors = {node:set(encoded_graphs[1].predecessors(node)) for node in list(encoded_graphs[1].nodes())}
#         directed_G.out_neighbors = {node:set(encoded_graphs[0].neighbors(node)) for node in list(encoded_graphs[0].nodes())}
#         directed_H.out_neighbors = {node:set(encoded_graphs[1].neighbors(node)) for node in list(encoded_graphs[1].nodes())}
#     else:
#         undirected_G.neighbors = {node:set(encoded_graphs[0].neighbors(node)) for node in list(encoded_graphs[0].nodes())}
#         undirected_H.neighbors = {node:set(encoded_graphs[1].neighbors(node)) for node in list(encoded_graphs[1].nodes())}

#     # get total order for VF2-like analysis
#     # NOTE: for connected extensions this should be given by concentric neighborhoods around the anchor with a multisource BFS
#     undirected_copy_H = deepcopy(encoded_graphs[1])
#     if(nx.is_directed(nx_H)):
#         undirected_copy_H = undirected_copy_H.to_undirected()
#     connectivity_neighbors = {node:set(undirected_copy_H.neighbors(node)) for node in list(undirected_copy_H.nodes())}
#     current_level = [node2 for (node1, node2) in encoded_anchor]
#     for each_pair in connectivity_neighbors:
#         node = each_pair.first
#         visited[node] = False

#     for node in current_level:
#         counter = counter + 1
#         total_order[node] = counter
#         visited[node] = True

#     while(not current_level.empty()):
#         # reinitialize next_level
#         next_level.clear()
#         # iterate getting immediate neighbors not already ordered
#         for node1 in current_level:
#             for node2 in connectivity_neighbors[node1]:
#                 # only assign oreder of not visited yet
#                 if(not visited[node2]):
#                     # increase and assign counter
#                     counter = counter + 1
#                     total_order[node2] = counter
#                     visited[node2] = True
#                     # level management
#                     next_level.push_back(node2)
#         # update nodes to be ordered
#         current_level.clear()
#         current_level = next_level

#     for each_pair in connectivity_neighbors:
#         node = each_pair.first
#         if(not visited[node]):
#             counter = counter + 1
#             total_order[node] = counter
#             visited[node] = True

#     # get expected order
#     expected_order = min([nx_G.order(), nx_H.order()])

#     # set recursion limit if recursive version was requested
#     if(not iterative_search):
#         scalation_value = 1.5
#         required_limit = max([nx_G.order(), nx_H.order()])
#         current_limit = getrecursionlimit()
#         if(current_limit < (scalation_value * required_limit)):
#             setrecursionlimit(int(scalation_value * required_limit))

#     # get maximum extensions
#     if(nx.is_directed(nx_G)):
#         if(iterative_search):
#             directed_maximum_connected_extensions_iterative(node_labels,
#                                                             edge_labels,
#                                                             all_extensions,
#                                                             expected_order,
#                                                             encoded_anchor,
#                                                             total_order,
#                                                             directed_G,
#                                                             directed_H,
#                                                             encoded_extensions)
#         else:
#             directed_maximum_connected_extensions_recursive(node_labels,
#                                                             edge_labels,
#                                                             all_extensions,
#                                                             expected_order,
#                                                             encoded_anchor,
#                                                             total_order,
#                                                             directed_G,
#                                                             directed_H,
#                                                             encoded_extensions)
#     else:
#         if(iterative_search):
#             undirected_maximum_connected_extensions_iterative(node_labels,
#                                                               edge_labels,
#                                                               all_extensions,
#                                                               expected_order,
#                                                               encoded_anchor,
#                                                               total_order,
#                                                               undirected_G,
#                                                               undirected_H,
#                                                               encoded_extensions)
#         else:
#             undirected_maximum_connected_extensions_recursive(node_labels,
#                                                               edge_labels,
#                                                               all_extensions,
#                                                               expected_order,
#                                                               encoded_anchor,
#                                                               total_order,
#                                                               undirected_G,
#                                                               undirected_H,
#                                                               encoded_extensions)
#     # decode maximum extensions
#     for each_extension in encoded_extensions:
#         extensions.append(decode_match(list(each_extension), encoded_node_names))

#     # check if the anchor was a good partial map
#     if((len(extensions[0]) == nx_G.order()) and (len(extensions[0]) == nx_H.order())):
#         good_anchor = True

#     # end of function
#     return(extensions, good_anchor)




# # functions - maximum connected extensions - undirected ########################




# # function: core routine of VF2-like undirected approach - iterative -----------
# # NOTE: an iterative DFS version of this algorithm can be implemented without a "visited"
# # list, since the total order given to the vertices of the second graph guarantees that
# # the search space is actually a search tree, and thus cannot have repeated states.
# cdef void undirected_maximum_connected_extensions_iterative(cpp_bool node_labels,
#                                                             cpp_bool edge_labels,
#                                                             cpp_bool all_extensions,
#                                                             size_t expected_order,
#                                                             cpp_vector[cpp_pair[int, int]] input_anchor,
#                                                             cpp_unordered_map[int, int] & total_order,
#                                                             partial_maps_undirected_graph & G,
#                                                             partial_maps_undirected_graph & H,
#                                                             cpp_vector[cpp_vector[cpp_pair[int, int]]] & all_matches) noexcept:

#     # local variables (cython)
#     cdef size_t new_score = 0
#     cdef size_t old_score = 0
#     cdef cpp_bool semantic_feasibility_res = True
#     cdef cpp_bool syntactic_feasibility_res = True
#     cdef cpp_pair[int, int] each_pair
#     cdef cpp_vector[cpp_pair[int, int]] new_match
#     cdef cpp_vector[cpp_pair[int, int]] candidates
#     cdef cpp_vector[cpp_pair[int, int]] current_match
#     cdef cpp_unordered_set[int] current_match_G
#     cdef cpp_unordered_set[int] current_match_H
#     cdef cpp_unordered_map[int, int] forward_match
#     cdef cpp_stack[cpp_vector[cpp_pair[int, int]]] dfs_stack

#     # initialize stack with (non-empty) input anchor match
#     dfs_stack.push(input_anchor)

#     # iterative DFS search
#     while(not dfs_stack.empty()):
#         # get current state and pop it
#         current_match.clear()
#         current_match = dfs_stack.top()
#         dfs_stack.pop()

#         # test initial match and improve if possible
#         if(all_matches.empty()):
#             if(not current_match.empty()):
#                 all_matches.push_back(current_match)
#                 # if complete match and only one then return
#                 if(G.nodes.size() == H.nodes.size()):
#                     if(current_match.size() == G.nodes.size()):
#                         if(not all_extensions):
#                             return
#         else:
#             # test improvement in matches
#             new_score = current_match.size()
#             old_score = all_matches[0].size()
#             # save match if it has the same score
#             if(new_score == old_score):
#                 all_matches.push_back(current_match)
#             # overwrite everything with new match if it improves score
#             if(new_score > old_score):
#                 all_matches.clear()
#                 all_matches.push_back(current_match)
#                 # if complete match and only one then return
#                 if(G.nodes.size() == H.nodes.size()):
#                     if(current_match.size() == G.nodes.size()):
#                         if(not all_extensions):
#                             return

#         # if not optimal yet then obtain available pairs
#         if(current_match.size() < expected_order):
#             # reinitialize variables used for filtering candidates
#             current_match_G.clear()
#             current_match_H.clear()
#             forward_match.clear()
#             candidates.clear()

#             # generate auxiliary structures
#             for each_pair in current_match:
#                 current_match_G.insert(each_pair.first)
#                 current_match_H.insert(each_pair.second)
#                 forward_match[each_pair.first] = each_pair.second

#             # get candidate pairs
#             candidates = undirected_candidates(current_match,
#                                                current_match_G,
#                                                current_match_H,
#                                                G.neighbors,
#                                                H.neighbors,
#                                                total_order)

#             # evaluate candidates
#             for each_pair in candidates:
#                 # evaluate sintactic feasibility
#                 syntactic_feasibility_res = syntactic_feasibility(each_pair.first,
#                                                                   each_pair.second,
#                                                                   current_match_G,
#                                                                   current_match_H,
#                                                                   forward_match,
#                                                                   G.neighbors,
#                                                                   H.neighbors)

#                 if(syntactic_feasibility_res):
#                     # evaluate semantic feasibility
#                     if(node_labels or edge_labels):
#                         semantic_feasibility_res = semantic_feasibility(node_labels,
#                                                                         edge_labels,
#                                                                         each_pair.first,
#                                                                         each_pair.second,
#                                                                         current_match_G,
#                                                                         forward_match,
#                                                                         G.nodes,
#                                                                         H.nodes,
#                                                                         G.neighbors,
#                                                                         G.edges,
#                                                                         H.edges)

#                     # push to stack if valid
#                     if(semantic_feasibility_res):
#                         # build new match
#                         new_match.clear()
#                         new_match = current_match
#                         new_match.push_back(each_pair)
#                         # add new valid candidate states
#                         dfs_stack.push(new_match)
#     # end of function




# # function: core routine of VF2-like undirected approach - recursive -----------
# cdef void undirected_maximum_connected_extensions_recursive(cpp_bool node_labels,
#                                                             cpp_bool edge_labels,
#                                                             cpp_bool all_extensions,
#                                                             size_t expected_order,
#                                                             cpp_vector[cpp_pair[int, int]] current_match,
#                                                             cpp_unordered_map[int, int] & total_order,
#                                                             partial_maps_undirected_graph & G,
#                                                             partial_maps_undirected_graph & H,
#                                                             cpp_vector[cpp_vector[cpp_pair[int, int]]] & all_matches) noexcept:

#     # local variables (cython)
#     cdef size_t new_score = 0
#     cdef size_t old_score = 0
#     cdef cpp_bool semantic_feasibility_res = True
#     cdef cpp_bool syntactic_feasibility_res = True
#     cdef cpp_pair[int, int] each_pair
#     cdef cpp_vector[cpp_pair[int, int]] new_match
#     cdef cpp_vector[cpp_pair[int, int]] candidates
#     cdef cpp_unordered_set[int] current_match_G
#     cdef cpp_unordered_set[int] current_match_H
#     cdef cpp_unordered_map[int, int] forward_match

#     # test initial match and improve if possible
#     if(all_matches.empty()):
#         if(not current_match.empty()):
#             all_matches.push_back(current_match)
#             # if complete match and only one then return
#             if(G.nodes.size() == H.nodes.size()):
#                 if(current_match.size() == G.nodes.size()):
#                     if(not all_extensions):
#                         return
#     else:
#         # test improvement in matches
#         new_score = current_match.size()
#         old_score = all_matches[0].size()
#         # save match if it has the same score
#         if(new_score == old_score):
#             all_matches.push_back(current_match)
#         # overwrite everything with new match if it improves score
#         if(new_score > old_score):
#             all_matches.clear()
#             all_matches.push_back(current_match)
#             # if complete match and only one then return
#             if(G.nodes.size() == H.nodes.size()):
#                 if(current_match.size() == G.nodes.size()):
#                     if(not all_extensions):
#                         return

#     # if not optimal yet then obtain available pairs
#     if(current_match.size() < expected_order):
#         # generate auxiliary structures
#         for each_pair in current_match:
#             current_match_G.insert(each_pair.first)
#             current_match_H.insert(each_pair.second)
#             forward_match[each_pair.first] = each_pair.second

#         # get candidate pairs
#         candidates = undirected_candidates(current_match,
#                                            current_match_G,
#                                            current_match_H,
#                                            G.neighbors,
#                                            H.neighbors,
#                                            total_order)

#         # evaluate candidates
#         for each_pair in candidates:
#             # evaluate sintactic feasibility
#             syntactic_feasibility_res = syntactic_feasibility(each_pair.first,
#                                                               each_pair.second,
#                                                               current_match_G,
#                                                               current_match_H,
#                                                               forward_match,
#                                                               G.neighbors,
#                                                               H.neighbors)

#             if(syntactic_feasibility_res):
#                 # evaluate semantic feasibility
#                 if(node_labels or edge_labels):
#                     semantic_feasibility_res = semantic_feasibility(node_labels,
#                                                                     edge_labels,
#                                                                     each_pair.first,
#                                                                     each_pair.second,
#                                                                     current_match_G,
#                                                                     forward_match,
#                                                                     G.nodes,
#                                                                     H.nodes,
#                                                                     G.neighbors,
#                                                                     G.edges,
#                                                                     H.edges)

#                 if(semantic_feasibility_res):
#                     # build new match
#                     new_match.clear()
#                     new_match = current_match
#                     new_match.push_back(each_pair)
#                     # extend match
#                     undirected_maximum_connected_extensions_recursive(node_labels,
#                                                                       edge_labels,
#                                                                       all_extensions,
#                                                                       expected_order,
#                                                                       new_match,
#                                                                       total_order,
#                                                                       G,
#                                                                       H,
#                                                                       all_matches)

#                     # finish if only one complete extension was requested and it was already found
#                     if(G.nodes.size() == H.nodes.size()):
#                         # a superset of the anchor is always present in vector of all matches at this point
#                         if(all_matches[0].size() == G.nodes.size()):
#                             if(not all_extensions):
#                                 return
#     # end of function




# # function: get candidate pairs for undirected extension search ----------------
# cdef cpp_vector[cpp_pair[int, int]] undirected_candidates(cpp_vector[cpp_pair[int, int]] & current_match,
#                                                           cpp_unordered_set[int] & current_match_G,
#                                                           cpp_unordered_set[int] & current_match_H,
#                                                           cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_G,
#                                                           cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_H,
#                                                           cpp_unordered_map[int, int] & total_order) noexcept:

#     # output holders
#     cdef cpp_vector[cpp_pair[int, int]] candidate_pairs

#     # local variables (cython)
#     cdef int node = 0
#     cdef int node1 = 0
#     cdef int node2 = 0
#     cdef int reference_maximum = 0
#     cdef cpp_string comma
#     comma.push_back(44)
#     cdef cpp_string temp_string
#     cdef cpp_pair[int, int] temp_pair
#     cdef cpp_pair[int, int] each_pair
#     cdef cpp_vector[int] valid_G
#     cdef cpp_vector[int] valid_H
#     cdef cpp_unordered_set[cpp_string] candidate_pairs_member

#     # get maximum value of total order in match
#     for node in current_match_H:
#         if(total_order[node] > reference_maximum):
#             reference_maximum = total_order[node]

#     # get candidate pairs based on valid sets
#     for each_pair in current_match:
#         # reinitialize valid neighbors
#         valid_G.clear()
#         valid_H.clear()

#         # get valid neighbors in G
#         for node in neigh_G[each_pair.first]:
#             # if not yet in match
#             if(current_match_G.find(node) == current_match_G.end()):
#                 valid_G.push_back(node)

#         # get valid neighbors in H
#         for node in neigh_H[each_pair.second]:
#             # if total order greater than in match
#             if(total_order[node] > reference_maximum):
#                 # if not yet in match
#                 if(current_match_H.find(node) == current_match_H.end()):
#                     valid_H.push_back(node)

#         # make product of valid neighbors
#         if((not valid_G.empty()) and (not valid_H.empty())):
#             for node1 in valid_G:
#                 for node2 in valid_H:
#                     temp_string = to_string(node1) + comma + to_string(node2)
#                     if(candidate_pairs_member.find(temp_string) == candidate_pairs_member.end()):
#                         # add proper pair
#                         temp_pair.first = node1
#                         temp_pair.second = node2
#                         candidate_pairs.push_back(temp_pair)
#                         # add string version for constant look ups
#                         candidate_pairs_member.insert(temp_string)

#     # end of function
#     return(candidate_pairs)




# # functions - maximum connected extensions - directed ##########################




# # function: core routine of VF2-like directed approach - iterative -------------
# # NOTE: an iterative DFS version of this algorithm can be implemented without a "visited"
# # list, since the total order given to the vertices of the second graph guarantees that
# # the search space is actually a search tree, and thus cannot have repeated states.
# cdef void directed_maximum_connected_extensions_iterative(cpp_bool node_labels,
#                                                           cpp_bool edge_labels,
#                                                           cpp_bool all_extensions,
#                                                           size_t expected_order,
#                                                           cpp_vector[cpp_pair[int, int]] input_anchor,
#                                                           cpp_unordered_map[int, int] & total_order,
#                                                           partial_maps_directed_graph & G,
#                                                           partial_maps_directed_graph & H,
#                                                           cpp_vector[cpp_vector[cpp_pair[int, int]]] & all_matches) noexcept:

#     # local variables (cython)
#     cdef size_t new_score = 0
#     cdef size_t old_score = 0
#     cdef cpp_bool in_semantic_feasibility = True
#     cdef cpp_bool in_syntactic_feasibility = True
#     cdef cpp_bool out_semantic_feasibility = True
#     cdef cpp_bool out_syntactic_feasibility = True
#     cdef cpp_pair[int, int] each_pair
#     cdef cpp_vector[cpp_pair[int, int]] new_match
#     cdef cpp_vector[cpp_pair[int, int]] candidates
#     cdef cpp_vector[cpp_pair[int, int]] current_match
#     cdef cpp_unordered_set[int] current_match_G
#     cdef cpp_unordered_set[int] current_match_H
#     cdef cpp_unordered_map[int, int] forward_match
#     cdef cpp_stack[cpp_vector[cpp_pair[int, int]]] dfs_stack

#     # initialize stack with (non-empty) input anchor match
#     dfs_stack.push(input_anchor)

#     # iterative DFS search
#     while(not dfs_stack.empty()):
#         # get current state and pop it
#         current_match.clear()
#         current_match = dfs_stack.top()
#         dfs_stack.pop()

#         # test initial match and improve if possible
#         if(all_matches.empty()):
#             if(not current_match.empty()):
#                 all_matches.push_back(current_match)
#                 # if complete match and only one then return
#                 if(G.nodes.size() == H.nodes.size()):
#                     if(current_match.size() == G.nodes.size()):
#                         if(not all_extensions):
#                             return
#         else:
#             # test improvement in matches
#             new_score = current_match.size()
#             old_score = all_matches[0].size()
#             # save match if it has the same score
#             if(new_score == old_score):
#                 all_matches.push_back(current_match)
#             # overwrite everything with new match if it improves score
#             if(new_score > old_score):
#                 all_matches.clear()
#                 all_matches.push_back(current_match)
#                 # if complete match and only one then return
#                 if(G.nodes.size() == H.nodes.size()):
#                     if(current_match.size() == G.nodes.size()):
#                         if(not all_extensions):
#                             return

#         # if not optimal yet then obtain available pairs
#         if(current_match.size() < expected_order):
#             # reinitialize variables used for filtering candidates
#             current_match_G.clear()
#             current_match_H.clear()
#             forward_match.clear()
#             candidates.clear()

#             # generate auxiliary structures
#             for each_pair in current_match:
#                 current_match_G.insert(each_pair.first)
#                 current_match_H.insert(each_pair.second)
#                 forward_match[each_pair.first] = each_pair.second

#             # get candidate pairs
#             candidates = directed_candidates(current_match,
#                                              current_match_G,
#                                              current_match_H,
#                                              G.in_neighbors,
#                                              H.in_neighbors,
#                                              G.out_neighbors,
#                                              H.out_neighbors,
#                                              total_order)

#             # evaluate candidates
#             for each_pair in candidates:
#                 # evaluate syntactic feasibility with in-neighbors
#                 in_syntactic_feasibility = syntactic_feasibility(each_pair.first,
#                                                                  each_pair.second,
#                                                                  current_match_G,
#                                                                  current_match_H,
#                                                                  forward_match,
#                                                                  G.in_neighbors,
#                                                                  H.in_neighbors)

#                 if(in_syntactic_feasibility):
#                     # evaluate syntactic feasibility with out-neighbors
#                     out_syntactic_feasibility = syntactic_feasibility(each_pair.first,
#                                                                       each_pair.second,
#                                                                       current_match_G,
#                                                                       current_match_H,
#                                                                       forward_match,
#                                                                       G.out_neighbors,
#                                                                       H.out_neighbors)

#                     if(out_syntactic_feasibility):
#                         # evaluate semantic feasibility with in-neighbors
#                         if(node_labels or edge_labels):
#                             in_semantic_feasibility = semantic_feasibility(node_labels,
#                                                                            edge_labels,
#                                                                            each_pair.first,
#                                                                            each_pair.second,
#                                                                            current_match_G,
#                                                                            forward_match,
#                                                                            G.nodes,
#                                                                            H.nodes,
#                                                                            G.in_neighbors,
#                                                                            G.edges,
#                                                                            H.edges)

#                         if(in_semantic_feasibility):
#                             # evaluate semantic feasibility with out-neighbors
#                             if(node_labels or edge_labels):
#                                 out_semantic_feasibility = semantic_feasibility(node_labels,
#                                                                                 edge_labels,
#                                                                                 each_pair.first,
#                                                                                 each_pair.second,
#                                                                                 current_match_G,
#                                                                                 forward_match,
#                                                                                 G.nodes,
#                                                                                 H.nodes,
#                                                                                 G.out_neighbors,
#                                                                                 G.edges,
#                                                                                 H.edges)

#                             # push to stack if valid
#                             if(out_semantic_feasibility):
#                                 # build new match
#                                 new_match.clear()
#                                 new_match = current_match
#                                 new_match.push_back(each_pair)
#                                 # add new valid candidate states
#                                 dfs_stack.push(new_match)
#     # end of function




# # function: core routine of VF2-like directed approach - recursive -------------
# cdef void directed_maximum_connected_extensions_recursive(cpp_bool node_labels,
#                                                           cpp_bool edge_labels,
#                                                           cpp_bool all_extensions,
#                                                           size_t expected_order,
#                                                           cpp_vector[cpp_pair[int, int]] current_match,
#                                                           cpp_unordered_map[int, int] & total_order,
#                                                           partial_maps_directed_graph & G,
#                                                           partial_maps_directed_graph & H,
#                                                           cpp_vector[cpp_vector[cpp_pair[int, int]]] & all_matches) noexcept:

#     # local variables (cython)
#     cdef size_t new_score = 0
#     cdef size_t old_score = 0
#     cdef cpp_bool in_semantic_feasibility = True
#     cdef cpp_bool in_syntactic_feasibility = True
#     cdef cpp_bool out_semantic_feasibility = True
#     cdef cpp_bool out_syntactic_feasibility = True
#     cdef cpp_pair[int, int] each_pair
#     cdef cpp_vector[cpp_pair[int, int]] new_match
#     cdef cpp_vector[cpp_pair[int, int]] candidates
#     cdef cpp_unordered_set[int] current_match_G
#     cdef cpp_unordered_set[int] current_match_H
#     cdef cpp_unordered_map[int, int] forward_match

#     # test initial match and consecutive matches
#     if(all_matches.empty()):
#         if(not current_match.empty()):
#             all_matches.push_back(current_match)
#             # if complete match and only one then return
#             if(G.nodes.size() == H.nodes.size()):
#                 if(current_match.size() == G.nodes.size()):
#                     if(not all_extensions):
#                         return
#     else:
#         # test improvement in matches
#         new_score = current_match.size()
#         old_score = all_matches[0].size()
#         # save match if it has the same score
#         if(new_score == old_score):
#             all_matches.push_back(current_match)
#         # overwrite with new match if it improves score
#         if(new_score > old_score):
#             all_matches.clear()
#             all_matches.push_back(current_match)
#             # if complete match and only one then return
#             if(G.nodes.size() == H.nodes.size()):
#                 if(current_match.size() == G.nodes.size()):
#                     if(not all_extensions):
#                         return

#     # if not optimal yet then obtain available pairs
#     if(current_match.size() < expected_order):
#         # generate auxiliary structures
#         for each_pair in current_match:
#             current_match_G.insert(each_pair.first)
#             current_match_H.insert(each_pair.second)
#             forward_match[each_pair.first] = each_pair.second

#         # get candidate pairs
#         candidates = directed_candidates(current_match,
#                                          current_match_G,
#                                          current_match_H,
#                                          G.in_neighbors,
#                                          H.in_neighbors,
#                                          G.out_neighbors,
#                                          H.out_neighbors,
#                                          total_order)

#         # evaluate candidates
#         for each_pair in candidates:
#             # evaluate syntactic feasibility with in-neighbors
#             in_syntactic_feasibility = syntactic_feasibility(each_pair.first,
#                                                              each_pair.second,
#                                                              current_match_G,
#                                                              current_match_H,
#                                                              forward_match,
#                                                              G.in_neighbors,
#                                                              H.in_neighbors)

#             if(in_syntactic_feasibility):
#                 # evaluate syntactic feasibility with out-neighbors
#                 out_syntactic_feasibility = syntactic_feasibility(each_pair.first,
#                                                                   each_pair.second,
#                                                                   current_match_G,
#                                                                   current_match_H,
#                                                                   forward_match,
#                                                                   G.out_neighbors,
#                                                                   H.out_neighbors)

#                 if(out_syntactic_feasibility):
#                     # evaluate semantic feasibility with in-neighbors
#                     if(node_labels or edge_labels):
#                         in_semantic_feasibility = semantic_feasibility(node_labels,
#                                                                        edge_labels,
#                                                                        each_pair.first,
#                                                                        each_pair.second,
#                                                                        current_match_G,
#                                                                        forward_match,
#                                                                        G.nodes,
#                                                                        H.nodes,
#                                                                        G.in_neighbors,
#                                                                        G.edges,
#                                                                        H.edges)

#                     if(in_semantic_feasibility):
#                         # evaluate semantic feasibility with out-neighbors
#                         if(node_labels or edge_labels):
#                             out_semantic_feasibility = semantic_feasibility(node_labels,
#                                                                             edge_labels,
#                                                                             each_pair.first,
#                                                                             each_pair.second,
#                                                                             current_match_G,
#                                                                             forward_match,
#                                                                             G.nodes,
#                                                                             H.nodes,
#                                                                             G.out_neighbors,
#                                                                             G.edges,
#                                                                             H.edges)

#                         # push to stack if valid
#                         if(out_semantic_feasibility):
#                             # build new match
#                             new_match.clear()
#                             new_match = current_match
#                             new_match.push_back(each_pair)
#                             # extend match
#                             directed_maximum_connected_extensions_recursive(node_labels,
#                                                                             edge_labels,
#                                                                             all_extensions,
#                                                                             expected_order,
#                                                                             new_match,
#                                                                             total_order,
#                                                                             G,
#                                                                             H,
#                                                                             all_matches)

#                             # finish if only one complete extension was requested and it was already found
#                             if(G.nodes.size() == H.nodes.size()):
#                                 # anchor is always present in vector of all matches at this point
#                                 if(all_matches[0].size() == G.nodes.size()):
#                                     if(not all_extensions):
#                                         return
#     # end of function




# # function: get candidate pairs for directed extension search ------------------
# cdef cpp_vector[cpp_pair[int, int]] directed_candidates(cpp_vector[cpp_pair[int, int]] & current_match,
#                                                         cpp_unordered_set[int] & current_match_G,
#                                                         cpp_unordered_set[int] & current_match_H,
#                                                         cpp_unordered_map[int, cpp_unordered_set[int]] & in_neigh_G,
#                                                         cpp_unordered_map[int, cpp_unordered_set[int]] & in_neigh_H,
#                                                         cpp_unordered_map[int, cpp_unordered_set[int]] & out_neigh_G,
#                                                         cpp_unordered_map[int, cpp_unordered_set[int]] & out_neigh_H,
#                                                         cpp_unordered_map[int, int] & total_order) noexcept:

#     # output holders
#     cdef cpp_vector[cpp_pair[int, int]] candidate_pairs

#     # local variables (cython)
#     cdef int node = 0
#     cdef int node1 = 0
#     cdef int node2 = 0
#     cdef int reference_maximum = 0
#     cdef cpp_string comma
#     comma.push_back(44)
#     cdef cpp_string temp_string
#     cdef cpp_pair[int, int] temp_pair
#     cdef cpp_pair[int, int] each_pair
#     cdef cpp_vector[int] valid_G
#     cdef cpp_vector[int] valid_H
#     cdef cpp_unordered_set[cpp_string] candidate_pairs_member

#     # get maximum value of total order in match
#     for node in current_match_H:
#         if(total_order[node] > reference_maximum):
#             reference_maximum = total_order[node]

#     # get candidates based on in-neighbors
#     for each_pair in current_match:
#         # reinitialize valid neighbors
#         valid_G.clear()
#         valid_H.clear()

#         # get valid in-neighbors in G
#         for node in in_neigh_G[each_pair.first]:
#             # if not yet in match
#             if(current_match_G.find(node) == current_match_G.end()):
#                 valid_G.push_back(node)

#         # get valid in-neighbors in H
#         for node in in_neigh_H[each_pair.second]:
#             # if total order greater than in match
#             if(total_order[node] > reference_maximum):
#                 # if not yet in match
#                 if(current_match_H.find(node) == current_match_H.end()):
#                     valid_H.push_back(node)

#         # make product of valid in-neighbors
#         if((not valid_G.empty()) and (not valid_H.empty())):
#             for node1 in valid_G:
#                 for node2 in valid_H:
#                     temp_string = to_string(node1) + comma + to_string(node2)
#                     if(candidate_pairs_member.find(temp_string) == candidate_pairs_member.end()):
#                         # add proper pair
#                         temp_pair.first = node1
#                         temp_pair.second = node2
#                         candidate_pairs.push_back(temp_pair)
#                         # add string version for constant look ups
#                         candidate_pairs_member.insert(temp_string)

#     # alternatively get candidate pairs from out-neighbors
#     if(candidate_pairs.empty()):
#         # get candidates based on valid sets
#         for each_pair in current_match:
#             # reinitialize valid neighbors
#             valid_G.clear()
#             valid_H.clear()

#             # get valid out-neighbors in G
#             for node in out_neigh_G[each_pair.first]:
#                 # if not yet in match
#                 if(current_match_G.find(node) == current_match_G.end()):
#                     valid_G.push_back(node)

#             # get valid out-neighbors in H
#             for node in out_neigh_H[each_pair.second]:
#                 # if total order greater than in match
#                 if(total_order[node] > reference_maximum):
#                     # if not yet in match
#                     if(current_match_H.find(node) == current_match_H.end()):
#                         valid_H.push_back(node)

#             # make product of valid out-neighbors
#             if((not valid_G.empty()) and (not valid_H.empty())):
#                 for node1 in valid_G:
#                     for node2 in valid_H:
#                         temp_string = to_string(node1) + comma + to_string(node2)
#                         if(candidate_pairs_member.find(temp_string) == candidate_pairs_member.end()):
#                             # add proper pair
#                             temp_pair.first = node1
#                             temp_pair.second = node2
#                             candidate_pairs.push_back(temp_pair)
#                             # add string version for constant look ups
#                             candidate_pairs_member.insert(temp_string)

#     # end of function
#     return(candidate_pairs)




# # functions - feasability of matches - undirected and directed #################




# # function: evaluate syntactic feasability for connected extension -------------
# cdef cpp_bool syntactic_feasibility(int node1,
#                                     int node2,
#                                     cpp_unordered_set[int] & current_match_G,
#                                     cpp_unordered_set[int] & current_match_H,
#                                     cpp_unordered_map[int, int] & forward_match,
#                                     cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_G,
#                                     cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_H) noexcept:

#     # local variables (cython)
#     cdef int node = 0
#     cdef int mapped = 0
#     cdef cpp_unordered_set[int] neighbors_match_G
#     cdef cpp_unordered_set[int] neighbors_match_H

#     # loop-consistency-test
#     if(neigh_G[node1].find(node1) != neigh_G[node1].end()):
#         if(neigh_H[node2].find(node2) == neigh_H[node2].end()):
#             # node1 has a loop in G but node2 has no loop in H
#             return(False)
#     if(neigh_G[node1].find(node1) == neigh_G[node1].end()):
#         if(neigh_H[node2].find(node2) != neigh_H[node2].end()):
#             # node1 has no loop in G but node2 has a loop in H
#             return(False)

#     # look ahead 0: consistency of neighbors in match
#     for node in neigh_G[node1]:
#         if(current_match_G.find(node) != current_match_G.end()):
#             mapped = forward_match[node]
#             neighbors_match_G.insert(mapped)

#     for node in neigh_H[node2]:
#         if(current_match_H.find(node) != current_match_H.end()):
#             neighbors_match_H.insert(node)

#     if(neighbors_match_G.size() != neighbors_match_H.size()):
#         # one node has more neighbors in the match than the other
#         return(False)
#     else:
#         for mapped in neighbors_match_G:
#             if(neighbors_match_H.find(mapped) == neighbors_match_H.end()):
#                 # the neighbors dont respect the match
#                 return(False)

#     # end of function
#     return(True)




# # function: evaluate semantic feasability for connected extension --------------
# cdef cpp_bool semantic_feasibility(cpp_bool node_labels,
#                                    cpp_bool edge_labels,
#                                    int node1,
#                                    int node2,
#                                    cpp_unordered_set[int] & current_match_G,
#                                    cpp_unordered_map[int, int] & forward_match,
#                                    cpp_unordered_map[int, int] & nodes_G,
#                                    cpp_unordered_map[int, int] & nodes_H,
#                                    cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_G,
#                                    cpp_unordered_map[cpp_string, int] & edges_G,
#                                    cpp_unordered_map[cpp_string, int] & edges_H) noexcept:

#     # local variables (cython)
#     cdef int node = 0
#     cdef cpp_string comma
#     comma.push_back(44)
#     cdef cpp_string labeled_edge_G
#     cdef cpp_string labeled_edge_H

#     if(node_labels):
#         # compare vertex-labels
#         if(nodes_G[node1] != nodes_H[node2]):
#             return(False)

#     if(edge_labels):
#         # compare loop-labels
#         if(neigh_G[node1].find(node1) != neigh_G[node1].end()):
#             # loop in G
#             labeled_edge_G = to_string(node1) + comma + to_string(node1)
#             # loop in H
#             labeled_edge_H = to_string(node2) + comma + to_string(node2)
#             # compare edge labels
#             if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
#                 return(False)

#         # compare non-loop edge-labels
#         for node in neigh_G[node1]:
#             if(current_match_G.find(node) != current_match_G.end()):
#                 # edge in G with only one end in match
#                 labeled_edge_G = to_string(node1) + comma + to_string(node)
#                 # edge in H with only one end in match
#                 labeled_edge_H = to_string(node2) + comma + to_string(forward_match[node])
#                 # compare edge labels
#                 if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
#                     return(False)

#     # end of function
#     return(True)




# ################################################################################
# ################################################################################

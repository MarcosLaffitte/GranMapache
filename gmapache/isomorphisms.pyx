################################################################################
#                                                                              #
# - GranMapache: GRAphs-and-Networks MAPping Applications                      #
#   with Cython and HEuristics                                                 #
#                                                                              #
# - Module: isomorphisms                                                       #
#                                                                              #
# - Description: analysis of isomorphisms and automorphisms.                   #
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
from .integerization import encode_graphs, decode_graphs, decode_match




# C/C++ structs ################################################################




# struct: undirected graph ---------------------------------------------------
# using unordered maps and unordered sets for scalability
cdef struct isomorphisms_undirected_graph:
    cpp_unordered_map[int, int] nodes
    cpp_unordered_map[cpp_string, int] edges
    cpp_unordered_map[int, cpp_unordered_set[int]] neighbors




# struct: directed graph -------------------------------------------------------
# using unordered maps and unordered sets for scalability
cdef struct isomorphisms_directed_graph:
    cpp_unordered_map[int, int] nodes
    cpp_unordered_map[cpp_string, int] edges
    cpp_unordered_map[int, cpp_unordered_set[int]] in_neighbors
    cpp_unordered_map[int, cpp_unordered_set[int]] out_neighbors




# structure: container for information of candidate matches undirected ---------
# using unordered maps and unordered sets for scalability
cdef struct isomorphisms_candidates_struct_undirected:
    cpp_vector[cpp_pair[int, int]] candidates
    cpp_unordered_set[int] ring_G
    cpp_unordered_set[int] ring_H




# structure: container for information of candidate matches directed -----------
# using unordered maps and unordered sets for scalability
cdef struct isomorphisms_candidates_struct_directed:
    cpp_vector[cpp_pair[int, int]] candidates
    cpp_unordered_set[int] in_ring_G
    cpp_unordered_set[int] in_ring_H
    cpp_unordered_set[int] out_ring_G
    cpp_unordered_set[int] out_ring_H




# algorithms ###################################################################




# functions - search isomorphisms - wrapper ####################################




# function: callable wrapper for searching for isomorphisms --------------------
def search_isomorphisms(nx_G = nx.Graph(),           # can also be a networkx DiGraph
                        nx_H = nx.Graph(),           # can also be a networkx DiGraph
                        all_isomorphisms = False,    # by default stops when finding one isomorphism (if any)
                        iterative_search = True):    # by default an iterative search is used, otherwise a recursive version is called
    # description
    """
    > description: receives two networkx (di-)graphs G and H both directed or both
    undirected, with the same number of vertices and edges, and a boolean variable
    indicating if the function should finish when finding only one isomorphism from
    G to H (if any), or if it should search for all possible isomorphisms from G to
    H and return them. An additional boolean variable determines if the search should
    be done iteratively (by default), or recursively.

    > input:
    * nx_G - first networkx (di)graph being matched.
    * nx_H - second networkx (di)graph being matched.
    * all_isomorphisms - boolean variable indicating if the function should stop as soon
    as one isomorphism is found (if any) - this is the default behavior- or if it should
    search for all possible isomorphisms between the graphs.
    * iterative_search - boolean indicating if the iterative version of this algorithm
    should be used (the default), or if a recursive version of it should be used instead.

    > output:
    * isomorphisms - (possibly empty) list of isomorphisms, each as a list of 2-tuples
    (x, y) of nodes x from G and y from H representing the one-to-one correspondences
    preserving adjacency and labels.
    * found_isomorphism - boolean value indicating if any isomorphism was found and
    equivalently if G and H are isomorphic (labeled) (di-)graphs.

    > calls:
    * .integerization.encode_graphs
    * .integerization.decode_graphs
    * .integerization.decode_match
    * undirected_search_isomorphisms
    * directed_search_isomorphisms
    """

    # exception handling and input correctness
    test_undir = nx.Graph()
    test_dir = nx.DiGraph()
    test_bool = False
    # check that first argument is networkx graph or digraph
    if(type(nx_G) not in [type(test_undir), type(test_dir)]):
        raise(ValueError("gmapache: first argument must be a networkx graph or digraph."))
    # check that second argument is networkx graph or digraph
    if(type(nx_H) not in [type(test_undir), type(test_dir)]):
        raise(ValueError("gmapache: second argument must be a networkx graph or digraph."))
    # check that third argument is networkx graph or digraph
    if(type(all_isomorphisms) not in [type(test_bool)]):
        raise(ValueError("gmapache: third argument must be a boolean variable."))
    # check that fourth argument is networkx graph or digraph
    if(type(iterative_search) not in [type(test_bool)]):
        raise(ValueError("gmapache: fourth argument must be a boolean variable."))
    # check that input graphs are of the same type
    if((nx.is_directed(nx_G)) and (not nx.is_directed(nx_H))):
        raise(ValueError("gmapache: input graphs must be both directed or both undirected."))
    if((not nx.is_directed(nx_G)) and (nx.is_directed(nx_H))):
        raise(ValueError("gmapache: input graphs must be both directed or both undirected."))
    # check that input graphs are not null graphs
    if((nx_G.order() == 0) or (nx_H.order() == 0)):
        raise(ValueError("gmapache: input graphs must have at least one node each."))
    # check that input graphs have the same number of nodes
    if(not nx_G.order() == nx_H.order()):
        raise(ValueError("gmapache: input graphs must have the same number of nodes."))
    # check that input graphs have the same number of edges
    if(not nx_G.size() == nx_H.size()):
        raise(ValueError("gmapache: input graphs must have the same number of edges."))

    # output holders
    cdef list isomorphisms = []
    found_isomorphism = False

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
    cdef cpp_vector[cpp_pair[int, int]] each_isomorphism
    cdef cpp_vector[cpp_pair[int, int]] current_match
    cdef cpp_vector[cpp_vector[cpp_pair[int, int]]] encoded_isomorphisms
    cdef cpp_unordered_map[int, int] total_order
    cdef isomorphisms_directed_graph directed_G
    cdef isomorphisms_directed_graph directed_H
    cdef isomorphisms_undirected_graph undirected_G
    cdef isomorphisms_undirected_graph undirected_H

    # local variables (python)
    cdef list encoded_graphs = []
    cdef dict info = dict()
    cdef dict encoded_node_names = dict()
    cdef dict encoded_node_label = dict()
    cdef dict encoded_edge_label = dict()

    # encode graphs
    encoded_graphs, encoded_node_names, encoded_node_label, encoded_edge_label = encode_graphs([nx_G, nx_H])

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
    expected_order = nx_G.order()

    # set recursion limit
    if(not iterative_search):
        scalation_value = 1.5
        required_limit = nx_G.order()
        current_limit = getrecursionlimit()
        if(current_limit < (scalation_value * required_limit)):
            setrecursionlimit(int(scalation_value * required_limit))

    # prepare current match
    current_match.clear()

    # evaluate isomorphisms
    if(nx.is_directed(nx_G)):
        if(iterative_search):
            directed_search_isomorphisms_iterative(all_isomorphisms,
                                                   expected_order,
                                                   current_match,
                                                   total_order,
                                                   directed_G,
                                                   directed_H,
                                                   encoded_isomorphisms)
        else:
            directed_search_isomorphisms_recursive(all_isomorphisms,
                                                   expected_order,
                                                   current_match,
                                                   total_order,
                                                   directed_G,
                                                   directed_H,
                                                   encoded_isomorphisms)
    else:
        if(iterative_search):
            undirected_search_isomorphisms_iterative(all_isomorphisms,
                                                     expected_order,
                                                     current_match,
                                                     total_order,
                                                     undirected_G,
                                                     undirected_H,
                                                     encoded_isomorphisms)
        else:
            undirected_search_isomorphisms_recursive(all_isomorphisms,
                                                     expected_order,
                                                     current_match,
                                                     total_order,
                                                     undirected_G,
                                                     undirected_H,
                                                     encoded_isomorphisms)

    # decode isomorphisms
    for each_isomorphism in encoded_isomorphisms:
        isomorphisms.append(decode_match(list(each_isomorphism), encoded_node_names))

    # check if the graphs were isomorphic
    if(len(isomorphisms) > 0):
        found_isomorphism = True

    # end of function
    return(isomorphisms, found_isomorphism)




# functions - search isomorphisms - undirected #################################




# function: core routine of VF2-like undirected approach - iterative -----------
# NOTE: an iterative DFS version of this algorithm can be implemented without a "visited"
# list, since the total order given to the vertices of the second graph guarantees that
# the search space is actually a search tree, and thus cannot have repeated states.
cdef void undirected_search_isomorphisms_iterative(cpp_bool all_isomorphisms,
                                                   size_t expected_order,
                                                   cpp_vector[cpp_pair[int, int]] input_anchor,
                                                   cpp_unordered_map[int, int] & total_order,
                                                   isomorphisms_undirected_graph & G,
                                                   isomorphisms_undirected_graph & H,
                                                   cpp_vector[cpp_vector[cpp_pair[int, int]]] & all_matches) noexcept:

    # local variables (cython)
    cdef size_t new_score = 0
    cdef size_t old_score = 0
    cdef cpp_bool semantic_feasibility_res = False
    cdef cpp_bool syntactic_feasibility_res = False
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[cpp_pair[int, int]] new_match
    cdef cpp_vector[cpp_pair[int, int]] current_match
    cdef cpp_unordered_set[int] current_match_G
    cdef cpp_unordered_set[int] current_match_H
    cdef cpp_unordered_map[int, int] forward_match
    cdef cpp_stack[cpp_vector[cpp_pair[int, int]]] dfs_stack
    cdef isomorphisms_candidates_struct_undirected candidates_struct

    # initialize stack with (non-empty) input anchor match
    dfs_stack.push(input_anchor)

    # iterative DFS search
    while(not dfs_stack.empty()):
        # get current state and pop it
        current_match.clear()
        current_match = dfs_stack.top()
        dfs_stack.pop()

        # save if isomorphism was reached
        # NOTE: in the future we might want score-optimal isomorphisms
        if(current_match.size() == expected_order):
            all_matches.push_back(current_match)
            if(not all_isomorphisms):
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
                                                      G.nodes,
                                                      H.nodes,
                                                      G.neighbors,
                                                      H.neighbors,
                                                      total_order)

            # evaluate candidates
            for each_pair in candidates_struct.candidates:
                # evaluate sintactic feasibility
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
                    semantic_feasibility_res = semantic_feasibility(each_pair.first,
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
cdef void undirected_search_isomorphisms_recursive(cpp_bool all_isomorphisms,
                                                   size_t expected_order,
                                                   cpp_vector[cpp_pair[int, int]] current_match,
                                                   cpp_unordered_map[int, int] & total_order,
                                                   isomorphisms_undirected_graph & G,
                                                   isomorphisms_undirected_graph & H,
                                                   cpp_vector[cpp_vector[cpp_pair[int, int]]] & all_matches) noexcept:

    # local variables (cython)
    cdef size_t new_score = 0
    cdef size_t old_score = 0
    cdef cpp_bool semantic_feasibility_res = False
    cdef cpp_bool syntactic_feasibility_res = False
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[cpp_pair[int, int]] new_match
    cdef cpp_unordered_set[int] current_match_G
    cdef cpp_unordered_set[int] current_match_H
    cdef cpp_unordered_map[int, int] forward_match
    cdef isomorphisms_candidates_struct_undirected candidates_struct

    # save if isomorphism was reached
    # NOTE: in the future we might want score-optimal isomorphisms
    if(current_match.size() == expected_order):
        all_matches.push_back(current_match)
        if(not all_isomorphisms):
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
                                                  G.nodes,
                                                  H.nodes,
                                                  G.neighbors,
                                                  H.neighbors,
                                                  total_order)

        # evaluate candidates
        for each_pair in candidates_struct.candidates:
            # evaluate sintactic feasibility
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
                semantic_feasibility_res = semantic_feasibility(each_pair.first,
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
                    undirected_search_isomorphisms_recursive(all_isomorphisms,
                                                             expected_order,
                                                             new_match,
                                                             total_order,
                                                             G,
                                                             H,
                                                             all_matches)

                    # finish if only one isosmorphism was requested and it was already found
                    if(not all_matches.empty()):
                        if(not all_isomorphisms):
                            return
    # end of function




# function: get candidate pairs for undirected isomorphism search --------------
cdef isomorphisms_candidates_struct_undirected undirected_candidates(size_t expected_order,
                                                                     cpp_vector[cpp_pair[int, int]] & current_match,
                                                                     cpp_unordered_set[int] & current_match_G,
                                                                     cpp_unordered_set[int] & current_match_H,
                                                                     cpp_unordered_map[int, int] & nodes_G,
                                                                     cpp_unordered_map[int, int] & nodes_H,
                                                                     cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_G,
                                                                     cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_H,
                                                                     cpp_unordered_map[int, int] & total_order) noexcept:

    # output holders
    cdef isomorphisms_candidates_struct_undirected candidates_struct

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

    # alternatively get all unmatched nodes also not adjacent with match
    if(candidate_pairs.empty()):
        # reinitialize valid neighbors
        valid_G.clear()
        valid_H.clear()

        # get alternatives in G
        for each_pair in nodes_G:
            node = each_pair.first
            if(current_match_G.find(node) == current_match_G.end()):
                # should only match vertices outside of ring at this point
                if(ring_G.find(node) == ring_G.end()):
                    valid_G.push_back(node)

        # get alternatives in H
        for each_pair in nodes_H:
            node = each_pair.first
            if(current_match_H.find(node) == current_match_H.end()):
                # should only match vertices outside of ring at this point
                if(ring_H.find(node) == ring_H.end()):
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

    # pack return structure
    candidates_struct.candidates = candidate_pairs
    candidates_struct.ring_G = ring_G
    candidates_struct.ring_H = ring_H
    # end of function
    return(candidates_struct)




# functions - search isomorphisms - directed ###################################




# function: core routine of VF2-like directed approach - iterative -------------
# NOTE: an iterative DFS version of this algorithm can be implemented without a "visited"
# list, since the total order given to the vertices of the second graph guarantees that
# the search space is actually a search tree, and thus cannot have repeated states.
cdef void directed_search_isomorphisms_iterative(cpp_bool all_isomorphisms,
                                                 size_t expected_order,
                                                 cpp_vector[cpp_pair[int, int]] input_anchor,
                                                 cpp_unordered_map[int, int] & total_order,
                                                 isomorphisms_directed_graph & G,
                                                 isomorphisms_directed_graph & H,
                                                 cpp_vector[cpp_vector[cpp_pair[int, int]]] & all_matches) noexcept:

    # local variables (cython)
    cdef size_t new_score = 0
    cdef size_t old_score = 0
    cdef cpp_bool in_semantic_feasibility = False
    cdef cpp_bool in_syntactic_feasibility = False
    cdef cpp_bool out_semantic_feasibility = False
    cdef cpp_bool out_syntactic_feasibility = False
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[cpp_pair[int, int]] new_match
    cdef cpp_vector[cpp_pair[int, int]] current_match
    cdef cpp_unordered_set[int] current_match_G
    cdef cpp_unordered_set[int] current_match_H
    cdef cpp_unordered_map[int, int] forward_match
    cdef cpp_stack[cpp_vector[cpp_pair[int, int]]] dfs_stack
    cdef isomorphisms_candidates_struct_directed candidates_struct

    # initialize stack with (non-empty) input anchor match
    dfs_stack.push(input_anchor)

    # iterative DFS search
    while(not dfs_stack.empty()):
        # get current state and pop it
        current_match.clear()
        current_match = dfs_stack.top()
        dfs_stack.pop()

        # save if isomorphism was reached
        # NOTE: in the future we might want score-optimal isomorphisms
        if(current_match.size() == expected_order):
            all_matches.push_back(current_match)
            if(not all_isomorphisms):
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
                                                    G.nodes,
                                                    H.nodes,
                                                    G.in_neighbors,
                                                    H.in_neighbors,
                                                    G.out_neighbors,
                                                    H.out_neighbors,
                                                    total_order)

            # evaluate candidates
            for each_pair in candidates_struct.candidates:
                # evaluate sintactic feasibility of in-neighbors
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
                    # evaluate sintactic feasibility of out-neighbors
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
                        in_semantic_feasibility = semantic_feasibility(each_pair.first,
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
                            out_semantic_feasibility = semantic_feasibility(each_pair.first,
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
cdef void directed_search_isomorphisms_recursive(cpp_bool all_isomorphisms,
                                                 size_t expected_order,
                                                 cpp_vector[cpp_pair[int, int]] current_match,
                                                 cpp_unordered_map[int, int] & total_order,
                                                 isomorphisms_directed_graph & G,
                                                 isomorphisms_directed_graph & H,
                                                 cpp_vector[cpp_vector[cpp_pair[int, int]]] & all_matches) noexcept:

    # local variables (cython)
    cdef size_t new_score = 0
    cdef size_t old_score = 0
    cdef cpp_bool in_semantic_feasibility = False
    cdef cpp_bool in_syntactic_feasibility = False
    cdef cpp_bool out_semantic_feasibility = False
    cdef cpp_bool out_syntactic_feasibility = False
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[cpp_pair[int, int]] new_match
    cdef cpp_unordered_set[int] current_match_G
    cdef cpp_unordered_set[int] current_match_H
    cdef cpp_unordered_map[int, int] forward_match
    cdef isomorphisms_candidates_struct_directed candidates_struct

    # save if isomorphism was reached
    # NOTE: in the future we might want score-optimal isomorphisms
    if(current_match.size() == expected_order):
        all_matches.push_back(current_match)
        if(not all_isomorphisms):
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
                                                G.nodes,
                                                H.nodes,
                                                G.in_neighbors,
                                                H.in_neighbors,
                                                G.out_neighbors,
                                                H.out_neighbors,
                                                total_order)

        # evaluate candidates
        for each_pair in candidates_struct.candidates:
            # evaluate sintactic feasibility of in-neighbors
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
                # evaluate sintactic feasibility of out-neighbors
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
                    in_semantic_feasibility = semantic_feasibility(each_pair.first,
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
                        out_semantic_feasibility = semantic_feasibility(each_pair.first,
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
                            directed_search_isomorphisms_recursive(all_isomorphisms,
                                                                   expected_order,
                                                                   new_match,
                                                                   total_order,
                                                                   G,
                                                                   H,
                                                                   all_matches)

                            # finish if only one isosmorphism was requested and it was already found
                            if(not all_matches.empty()):
                                if(not all_isomorphisms):
                                    return
    # end of function




# function: get candidate pairs for directed isomorphism search ----------------
cdef isomorphisms_candidates_struct_directed directed_candidates(size_t expected_order,
                                                                 cpp_vector[cpp_pair[int, int]] & current_match,
                                                                 cpp_unordered_set[int] & current_match_G,
                                                                 cpp_unordered_set[int] & current_match_H,
                                                                 cpp_unordered_map[int, int] & nodes_G,
                                                                 cpp_unordered_map[int, int] & nodes_H,
                                                                 cpp_unordered_map[int, cpp_unordered_set[int]] & in_neigh_G,
                                                                 cpp_unordered_map[int, cpp_unordered_set[int]] & in_neigh_H,
                                                                 cpp_unordered_map[int, cpp_unordered_set[int]] & out_neigh_G,
                                                                 cpp_unordered_map[int, cpp_unordered_set[int]] & out_neigh_H,
                                                                 cpp_unordered_map[int, int] & total_order) noexcept:

    # output holders
    cdef isomorphisms_candidates_struct_directed candidates_struct

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

        # additionally save neighbors in out-ring around match in G
        for node in valid_G:
            # save if unrepeated
            if(out_ring_G.find(node) == out_ring_G.end()):
                out_ring_G.insert(node)

        # additionally save neighbors in out-ring around match in H
        for node in valid_H:
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

    # last alternative is to get all unmatched nodes also not adjacent with match
    if(candidate_pairs.empty()):
        # reinitialize valid neighbors
        valid_G.clear()
        valid_H.clear()

        # get alternatives in G
        for each_pair in nodes_G:
            node = each_pair.first
            if(current_match_G.find(node) == current_match_G.end()):
                # should only match vertices outside of ring at this point
                if(in_ring_G.find(node) == in_ring_G.end()):
                    if(out_ring_G.find(node) == out_ring_G.end()):
                        valid_G.push_back(node)

        # get alternatives in H
        for each_pair in nodes_H:
            node = each_pair.first
            if(current_match_H.find(node) == current_match_H.end()):
                # should only match vertices outside of ring at this point
                if(in_ring_H.find(node) == in_ring_H.end()):
                    if(out_ring_H.find(node) == out_ring_H.end()):
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

    # pack return structure
    candidates_struct.candidates = candidate_pairs
    candidates_struct.in_ring_G = in_ring_G
    candidates_struct.in_ring_H = in_ring_H
    candidates_struct.out_ring_G = out_ring_G
    candidates_struct.out_ring_H = out_ring_H
    # end of function
    return(candidates_struct)




# functions - feasability of matches - undirected and directed #################




# function: evaluate syntactic feasability for isomorphism search --------------
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




# function: evaluate semantic feasability for isomorphism search ---------------
cdef cpp_bool semantic_feasibility(int node1,
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

    # compare vertex-labels
    if(nodes_G[node1] != nodes_H[node2]):
        return(False)

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
            # get node label
            node_label = nodes_G[node]
            # get edge label (possibly loop label)
            labeled_edge_G = to_string(node1) + comma + to_string(node)
            edge_label = edges_G[labeled_edge_G]
            # count node label
            if(count_node_ring_G.find(node_label) != count_node_ring_G.end()):
                count_node_ring_G[node_label] = count_node_ring_G[node_label] + 1
            else:
                count_node_ring_G[node_label] = 1
            # count edge label
            if(count_edge_ring_G.find(edge_label) != count_edge_ring_G.end()):
                count_edge_ring_G[edge_label] = count_edge_ring_G[edge_label] + 1
            else:
                count_edge_ring_G[edge_label] = 1

        # count in H
        for node in neighbors_ring_H:
            # get node label
            node_label = nodes_H[node]
            # get edge label (possibly loop label)
            labeled_edge_H = to_string(node2) + comma + to_string(node)
            edge_label = edges_H[labeled_edge_H]
            # count node label
            if(count_node_ring_H.find(node_label) != count_node_ring_H.end()):
                count_node_ring_H[node_label] = count_node_ring_H[node_label] + 1
            else:
                count_node_ring_H[node_label] = 1
            # count edge label
            if(count_edge_ring_H.find(edge_label) != count_edge_ring_H.end()):
                count_edge_ring_H[edge_label] = count_edge_ring_H[edge_label] + 1
            else:
                count_edge_ring_H[edge_label] = 1

        # compare number of types of adjacent nodes
        if(count_node_ring_G.size() != count_node_ring_H.size()):
            return(False)

        # compare number of types of incident edges
        if(count_edge_ring_G.size() != count_edge_ring_H.size()):
            return(False)

        # compare types of adjacent nodes
        for each_pair in count_node_ring_G:
            if(count_node_ring_H.find(each_pair.first) == count_node_ring_H.end()):
                return(False)
            else:
                if(each_pair.second != count_node_ring_H[each_pair.first]):
                    return(False)

        # compare types of incident edges
        for each_pair in count_edge_ring_G:
            if(count_edge_ring_H.find(each_pair.first) == count_edge_ring_H.end()):
                return(False)
            else:
                if(each_pair.second != count_edge_ring_H[each_pair.first]):
                    return(False)

    # label look ahead 2: compare labels of extern neighbors (neither in match nor adjacent to match)
    if(not neighbors_extern_G.empty()):
        for node in neighbors_extern_G:
            # get node label
            node_label = nodes_G[node]
            # get edge label (possibly loop label)
            labeled_edge_G = to_string(node1) + comma + to_string(node)
            edge_label = edges_G[labeled_edge_G]
            # count node label
            if(count_node_extern_G.find(node_label) != count_node_extern_G.end()):
                count_node_extern_G[node_label] = count_node_extern_G[node_label] + 1
            else:
                count_node_extern_G[node_label] = 1
            # count edge label
            if(count_edge_extern_G.find(edge_label) != count_edge_extern_G.end()):
                count_edge_extern_G[edge_label] = count_edge_extern_G[edge_label] + 1
            else:
                count_edge_extern_G[edge_label] = 1

        for node in neighbors_extern_H:
            # get node label
            node_label = nodes_H[node]
            # get edge label (possibly loop label)
            labeled_edge_H = to_string(node2) + comma + to_string(node)
            edge_label = edges_H[labeled_edge_H]
            # count node label
            if(count_node_extern_H.find(node_label) != count_node_extern_H.end()):
                count_node_extern_H[node_label] = count_node_extern_H[node_label] + 1
            else:
                count_node_extern_H[node_label] = 1
            # count edge label
            if(count_edge_extern_H.find(edge_label) != count_edge_extern_H.end()):
                count_edge_extern_H[edge_label] = count_edge_extern_H[edge_label] + 1
            else:
                count_edge_extern_H[edge_label] = 1

        # compare number of types of adjacent nodes
        if(count_node_extern_G.size() != count_node_extern_H.size()):
            return(False)

        # compare number of types of incident edges
        if(count_edge_extern_G.size() != count_edge_extern_H.size()):
            return(False)

        # compare types of adjacent nodes
        for each_pair in count_node_extern_G:
            if(count_node_extern_H.find(each_pair.first) == count_node_extern_H.end()):
                return(False)
            else:
                if(each_pair.second != count_node_extern_H[each_pair.first]):
                    return(False)

        # compare types of incident edges
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

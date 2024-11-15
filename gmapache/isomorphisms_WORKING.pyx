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
from .integerization import encode_graphs, decode_graphs, decode_match



# C/C++ structs ################################################################



# struct: undirected graph ---------------------------------------------------
cdef struct isomorphisms_undirected_graph:
    cpp_map[int, int] nodes
    cpp_map[cpp_pair[int, int], int] edges
    cpp_map[int, cpp_vector[int]] neighbors



# struct: directed graph -------------------------------------------------------
cdef struct isomorphisms_directed_graph:
    cpp_map[int, int] nodes
    cpp_map[cpp_pair[int, int], int] edges
    cpp_map[int, cpp_vector[int]] in_neighbors
    cpp_map[int, cpp_vector[int]] out_neighbors



# structure: container for information of candidate matches undirected ---------
cdef struct isomorphisms_candidates_struct_undirected:
    cpp_vector[cpp_pair[int, int]] candidates
    cpp_vector[int] ring_G
    cpp_vector[int] ring_H



# structure: container for information of candidate matches directed -----------
cdef struct isomorphisms_candidates_struct_directed:
    cpp_vector[cpp_pair[int, int]] candidates
    cpp_vector[int] in_ring_G
    cpp_vector[int] in_ring_H
    cpp_vector[int] out_ring_G
    cpp_vector[int] out_ring_H



# algorithms ###################################################################



# functions - search isomorphisms - wrapper ####################################



# function: callable wrapper for searching for isomorphisms --------------------
def search_isomorphisms(nx_G = nx.Graph(),           # can also be a networkx DiGraph
                        nx_H = nx.Graph(),           # can also be a networkx DiGraph
                        all_isomorphisms = False):   # by default stops when finding one isomorphism (if any)
    # description
    """
    > description: receives two networkx (di-)graphs G and H both directed or both
    undirected, with the same number of vertices and edges, and a boolean variable
    indicating if the function should finish when finding only one isomorphism from
    G to H (if any), or if it should search for all possible isomorphisms from G to
    H and return them.

    > input:
    * nx_G - first networkx (di)graph being matched.
    * nx_H - second networkx (di)graph being matched.
    * all_isomorphisms - boolean variable indicating if the function should stop as soon
    as one isomorphism is found (if any) - this is the default behavior- or if it should
    search for all possible isomorphisms between the graphs.

    > output:
    * isomorphisms - (possibly empty) list of isomorphisms, each as a list of 2-tuples
    (x, y) of nodes x from G and y from H representing the one-to-one correspondences
    preserving adyacency and labels.
    * found_isomorphism - boolean value indicating if any isomorphism was found and
    equivalently if G and H are isomorphic labeled (di-)graphs.

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
    if(type(nx_G) not in [type(test_undir), type(test_dir)]):
        raise(ValueError("gmapache: first argument must be a networkx graph or digraph."))
    if(type(nx_H) not in [type(test_undir), type(test_dir)]):
        raise(ValueError("gmapache: second argument must be a networkx graph or digraph."))
    if(type(all_isomorphisms) not in [type(test_bool)]):
        raise(ValueError("gmapache: third argument must be a boolean variable."))
    if((nx.is_directed(nx_G)) and (not nx.is_directed(nx_H))):
        raise(ValueError("gmapache: input graphs must be both directed or both undirected."))
    if((not nx.is_directed(nx_G)) and (nx.is_directed(nx_H))):
        raise(ValueError("gmapache: input graphs must be both directed or both undirected."))
    if((nx_G.order() == 0) or (nx_H.order() == 0)):
        raise(ValueError("gmapache: input graphs must have at least one node each."))
    if(not nx_G.order() == nx_H.order()):
        raise(ValueError("gmapache: input graphs must have the same number of nodes."))
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
    cdef cpp_vector[cpp_pair[int, int]] each_isomorphism
    cdef cpp_vector[cpp_pair[int, int]] current_match
    cdef cpp_vector[cpp_vector[cpp_pair[int, int]]] encoded_isomorphisms
    cdef cpp_map[int, int] total_order
    cdef isomorphisms_undirected_graph undirected_G
    cdef isomorphisms_undirected_graph undirected_H
    cdef isomorphisms_directed_graph directed_G
    cdef isomorphisms_directed_graph directed_H
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
        directed_G.edges = {(node1, node2):info["GMEL"] for (node1, node2, info) in encoded_graphs[0].edges(data = True)}
        directed_H.edges = {(node1, node2):info["GMEL"] for (node1, node2, info) in encoded_graphs[1].edges(data = True)}
    else:
        undirected_G.edges = {tuple(sorted([node1, node2])):info["GMEL"] for (node1, node2, info) in encoded_graphs[0].edges(data = True)}
        undirected_H.edges = {tuple(sorted([node1, node2])):info["GMEL"] for (node1, node2, info) in encoded_graphs[1].edges(data = True)}
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
    for node in list(encoded_graphs[1].nodes()):
        counter = counter + 1
        total_order[node] = counter
    # get expected order
    expected_order = nx_G.order()
    # set recursion limit
    scalation_value = 1.5
    required_limit = nx_G.order()
    current_limit = getrecursionlimit()
    if(current_limit < (scalation_value * required_limit)):
        setrecursionlimit(int(scalation_value * required_limit))
    # prepare current match
    current_match.clear()
    # evaluate isomorphisms
    if(nx.is_directed(nx_G)):
        directed_search_isomorphisms(expected_order,
                                     current_match,
                                     all_isomorphisms,
                                     directed_G.nodes,
                                     directed_H.nodes,
                                     directed_G.in_neighbors,
                                     directed_H.in_neighbors,
                                     directed_G.out_neighbors,
                                     directed_H.out_neighbors,
                                     directed_G.edges,
                                     directed_H.edges,
                                     total_order,
                                     encoded_isomorphisms)
    else:
        undirected_search_isomorphisms(expected_order,
                                       current_match,
                                       all_isomorphisms,
                                       undirected_G.nodes,
                                       undirected_H.nodes,
                                       undirected_G.neighbors,
                                       undirected_H.neighbors,
                                       undirected_G.edges,
                                       undirected_H.edges,
                                       total_order,
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



# function: core routine of VF2-like undirected approach -----------------------
cdef void undirected_search_isomorphisms(size_t expected_order,
                                         cpp_vector[cpp_pair[int, int]] current_match,
                                         cpp_bool & all_isomorphisms,
                                         cpp_map[int, int] & nodes_G,
                                         cpp_map[int, int] & nodes_H,
                                         cpp_map[int, cpp_vector[int]] & neigh_G,
                                         cpp_map[int, cpp_vector[int]] & neigh_H,
                                         cpp_map[cpp_pair[int, int], int] & edges_G,
                                         cpp_map[cpp_pair[int, int], int] & edges_H,
                                         cpp_map[int, int] & total_order,
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
    cdef cpp_map[int, int] forward_match
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
            current_match_G.push_back(each_pair.first)
            current_match_H.push_back(each_pair.second)
            forward_match[each_pair.first] = each_pair.second
        # get candidates information
        candidates_struct = undirected_candidates(current_match,
                                                  current_match_G,
                                                  current_match_H,
                                                  nodes_G,
                                                  nodes_H,
                                                  neigh_G,
                                                  neigh_H,
                                                  total_order)
        # evaluate candidates
        for each_pair in candidates_struct.candidates:
            # evaluate sintactic feasibility
            syntactic_feasibility = undirected_syntactic_feasibility(each_pair.first,
                                                                     each_pair.second,
                                                                     candidates_struct.ring_G,
                                                                     candidates_struct.ring_H,
                                                                     current_match_G,
                                                                     current_match_H,
                                                                     forward_match,
                                                                     neigh_G,
                                                                     neigh_H)
            if(syntactic_feasibility):
                # evaluate semantic feasibility
                semantic_feasibility = undirected_semantic_feasibility(each_pair.first,
                                                                       each_pair.second,
                                                                       candidates_struct.ring_G,
                                                                       candidates_struct.ring_H,
                                                                       current_match_G,
                                                                       current_match_H,
                                                                       forward_match,
                                                                       nodes_G,
                                                                       nodes_H,
                                                                       neigh_G,
                                                                       neigh_H,
                                                                       edges_G,
                                                                       edges_H)
                if(semantic_feasibility):
                    # build new match
                    new_match = current_match
                    new_match.push_back(each_pair)
                    # extend match
                    undirected_search_isomorphisms(expected_order,
                                                   new_match,
                                                   all_isomorphisms,
                                                   nodes_G,
                                                   nodes_H,
                                                   neigh_G,
                                                   neigh_H,
                                                   edges_G,
                                                   edges_H,
                                                   total_order,
                                                   all_matches)
                    # finish if only one isosmorphism was requested and it was already found
                    if(not all_isomorphisms):
                        if(not all_matches.empty()):
                            return
    # end of function



# function: get candidate pairs for undirected isomorphism search --------------
cdef isomorphisms_candidates_struct_undirected undirected_candidates(cpp_vector[cpp_pair[int, int]] & current_match,
                                                                     cpp_vector[int] & current_match_G,
                                                                     cpp_vector[int] & current_match_H,
                                                                     cpp_map[int, int] & nodes_G,
                                                                     cpp_map[int, int] & nodes_H,
                                                                     cpp_map[int, cpp_vector[int]] & neigh_G,
                                                                     cpp_map[int, cpp_vector[int]] & neigh_H,
                                                                     cpp_map[int, int] & total_order) noexcept:
    # local variables
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int reference_minimum = -1
    cdef cpp_vector[int] ring_G
    cdef cpp_vector[int] ring_H
    cdef cpp_vector[int] valid_G
    cdef cpp_vector[int] valid_H
    cdef cpp_vector[cpp_pair[int, int]] candidate_pairs
    cdef cpp_vector[cpp_pair[int, int]] filtered_candidate_pairs
    cdef cpp_pair[int, int] each_pair
    cdef cpp_pair[int, int] temp_pair
    cdef isomorphisms_candidates_struct_undirected candidates_struct
    # get valid sets
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
        # additionally save neighbors in ring around match in G
        for node in valid_G:
            # save if unrepeated
            if(find(ring_G.begin(), ring_G.end(), node) == ring_G.end()):
                ring_G.push_back(node)
        # additionally save neighbors in ring around match in H
        for node in valid_H:
            # save if unrepeated
            if(find(ring_H.begin(), ring_H.end(), node) == ring_H.end()):
                ring_H.push_back(node)
    # alternatively get all unmatched nodes also not adjacent with match
    if(candidate_pairs.empty()):
        # reinitialize valid neighbors
        valid_G.clear()
        valid_H.clear()
        # get alternatives in G
        for each_pair in nodes_G:
            node = each_pair.first
            if(find(current_match_G.begin(), current_match_G.end(), node) == current_match_G.end()):
                # can only match vertices outside of ring at this point
                if(find(ring_G.begin(), ring_G.end(), node) == ring_G.end()):
                    valid_G.push_back(node)
        # get alternatives in H
        for each_pair in nodes_H:
            node = each_pair.first
            if(find(current_match_H.begin(), current_match_H.end(), node) == current_match_H.end()):
                # can only match vertices outside of ring at this point
                if(find(ring_H.begin(), ring_H.end(), node) == ring_H.end()):
                    valid_H.push_back(node)
        # make product of valid neighbors
        if((not valid_G.empty()) and (not valid_H.empty())):
            for node1 in valid_G:
                for node2 in valid_H:
                    temp_pair.first = node1
                    temp_pair.second = node2
                    # during isomorphism search the pair cannot be repeated at this point
                    candidate_pairs.push_back(temp_pair)
    # filter candidate pairs
    if(not candidate_pairs.empty()):
        # get minimum total order and filter candidates
        for each_pair in candidate_pairs:
            node = each_pair.second
            if(reference_minimum == -1):
                reference_minimum = total_order[node]
                filtered_candidate_pairs.push_back(each_pair)
            else:
                if(total_order[node] == reference_minimum):
                    filtered_candidate_pairs.push_back(each_pair)
                if(total_order[node] < reference_minimum):
                    reference_minimum = total_order[node]
                    filtered_candidate_pairs.clear()
                    filtered_candidate_pairs.push_back(each_pair)
        # copy just for clean output coding
        candidate_pairs = filtered_candidate_pairs
    # pack return structure
    candidates_struct.candidates = candidate_pairs
    candidates_struct.ring_G = ring_G
    candidates_struct.ring_H = ring_H
    # end of function
    return(candidates_struct)



# function: evaluate syntactic feasability for undirected isomorphisms ---------
cdef cpp_bool undirected_syntactic_feasibility(int node1,
                                               int node2,
                                               cpp_vector[int] & ring_G,
                                               cpp_vector[int] & ring_H,
                                               cpp_vector[int] & current_match_G,
                                               cpp_vector[int] & current_match_H,
                                               cpp_map[int, int] & forward_match,
                                               cpp_map[int, cpp_vector[int]] & neigh_G,
                                               cpp_map[int, cpp_vector[int]] & neigh_H) noexcept:
    # local variables
    cdef int node = 0
    cdef int mapped = 0
    cdef cpp_vector[int] neighbors_ring_G
    cdef cpp_vector[int] neighbors_ring_H
    cdef cpp_vector[int] neighbors_match_G
    cdef cpp_vector[int] neighbors_match_H
    cdef cpp_vector[int] neighbors_extern_G
    cdef cpp_vector[int] neighbors_extern_H
    # loop-consistency-test
    if(find(neigh_G[node1].begin(), neigh_G[node1].end(), node1) != neigh_G[node1].end()):
        if(find(neigh_H[node2].begin(), neigh_H[node2].end(), node2) == neigh_H[node2].end()):
            # node1 has a loop in G but node2 has no loop in H
            return(False)
    if(find(neigh_G[node1].begin(), neigh_G[node1].end(), node1) == neigh_G[node1].end()):
        if(find(neigh_H[node2].begin(), neigh_H[node2].end(), node2) != neigh_H[node2].end()):
            # node1 has no loop in G but node2 has a loop in H
            return(False)
    # obtain tripartition of neighbors
    for node in neigh_G[node1]:
        if(find(current_match_G.begin(), current_match_G.end(), node) != current_match_G.end()):
            # save map of neighbor to be compared later
            mapped = forward_match[node]
            neighbors_match_G.push_back(mapped)
        else:
            if(find(ring_G.begin(), ring_G.end(), node) != ring_G.end()):
                # save neighbor since we are just comparing numbers later
                neighbors_ring_G.push_back(node)
            else:
                # save neighbor since we are just comparing numbers later
                neighbors_extern_G.push_back(node)
    for node in neigh_H[node2]:
        if(find(current_match_H.begin(), current_match_H.end(), node) != current_match_H.end()):
            neighbors_match_H.push_back(node)
        else:
            if(find(ring_H.begin(), ring_H.end(), node) != ring_H.end()):
                neighbors_ring_H.push_back(node)
            else:
                neighbors_extern_H.push_back(node)
    # look ahead 0: consistency of neighbors in match
    if(neighbors_match_G.size() != neighbors_match_H.size()):
        # one node has more neighbors in the match than the other
        return(False)
    else:
        for mapped in neighbors_match_G:
            if(find(neighbors_match_H.begin(), neighbors_match_H.end(), mapped) == neighbors_match_H.end()):
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



# function: evaluate semantic feasability for undirected isomorphisms ----------
cdef cpp_bool undirected_semantic_feasibility(int node1,
                                              int node2,
                                              cpp_vector[int] & ring_G,
                                              cpp_vector[int] & ring_H,
                                              cpp_vector[int] & current_match_G,
                                              cpp_vector[int] & current_match_H,
                                              cpp_map[int, int] & forward_match,
                                              cpp_map[int, int] & nodes_G,
                                              cpp_map[int, int] & nodes_H,
                                              cpp_map[int, cpp_vector[int]] & neigh_G,
                                              cpp_map[int, cpp_vector[int]] & neigh_H,
                                              cpp_map[cpp_pair[int, int], int] & edges_G,
                                              cpp_map[cpp_pair[int, int], int] & edges_H) noexcept:
    # local variables
    cdef int node = 0
    cdef int mapped = 0
    cdef int node_label = 0
    cdef int edge_label = 0
    cdef cpp_pair[int, int] each_pair
    cdef cpp_pair[int, int] labeled_edge_G
    cdef cpp_pair[int, int] labeled_edge_H
    cdef cpp_vector[int] neighbors_ring_G
    cdef cpp_vector[int] neighbors_ring_H
    cdef cpp_vector[int] neighbors_match_G
    cdef cpp_vector[int] neighbors_match_H
    cdef cpp_vector[int] neighbors_extern_G
    cdef cpp_vector[int] neighbors_extern_H
    cdef cpp_map[int, int] count_node_ring_G
    cdef cpp_map[int, int] count_node_ring_H
    cdef cpp_map[int, int] count_edge_ring_G
    cdef cpp_map[int, int] count_edge_ring_H
    cdef cpp_map[int, int] count_node_extern_G
    cdef cpp_map[int, int] count_node_extern_H
    cdef cpp_map[int, int] count_edge_extern_G
    cdef cpp_map[int, int] count_edge_extern_H
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
    # obtain tripartition of neighbors in G
    for node in neigh_G[node1]:
        if(find(current_match_G.begin(), current_match_G.end(), node) != current_match_G.end()):
            neighbors_match_G.push_back(node)
        else:
            if(find(ring_G.begin(), ring_G.end(), node) != ring_G.end()):
                # save neighbor since we are just comparing numbers later
                neighbors_ring_G.push_back(node)
            else:
                # save neighbor since we are just comparing numbers later
                neighbors_extern_G.push_back(node)
    # label look ahead 0: compare non-loop edge-labels in match
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
    # obtain tripartition of neighbors in H
    for node in neigh_H[node2]:
        if(find(current_match_H.begin(), current_match_H.end(), node) != current_match_H.end()):
            neighbors_match_H.push_back(node)
        else:
            if(find(ring_H.begin(), ring_H.end(), node) != ring_H.end()):
                neighbors_ring_H.push_back(node)
            else:
                neighbors_extern_H.push_back(node)
    # label look ahead 1: compare labels of neighbors in ring (not in match but adjacent to match)
    if(not neighbors_ring_G.empty()):
        for node in neighbors_ring_G:
            # get node label
            node_label = nodes_G[node]
            # get edge label (possibly loop label)
            if(node1 <= node):
                labeled_edge_G.first = node1
                labeled_edge_G.second = node
            else:
                labeled_edge_G.first = node
                labeled_edge_G.second = node1
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
        for node in neighbors_ring_H:
            # get node label
            node_label = nodes_H[node]
            # get edge label (possibly loop label)
            if(node2 <= node):
                labeled_edge_H.first = node2
                labeled_edge_H.second = node
            else:
                labeled_edge_H.first = node
                labeled_edge_H.second = node2
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
            if(node1 <= node):
                labeled_edge_G.first = node1
                labeled_edge_G.second = node
            else:
                labeled_edge_G.first = node
                labeled_edge_G.second = node1
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
            if(node2 <= node):
                labeled_edge_H.first = node2
                labeled_edge_H.second = node
            else:
                labeled_edge_H.first = node
                labeled_edge_H.second = node2
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



# functions - search isomorphisms - directed ###################################



# function: core routine of VF2-like directed approach -------------------------
cdef void directed_search_isomorphisms(size_t expected_order,
                                       cpp_vector[cpp_pair[int, int]] current_match,
                                       cpp_bool & all_isomorphisms,
                                       cpp_map[int, int] & nodes_G,
                                       cpp_map[int, int] & nodes_H,
                                       cpp_map[int, cpp_vector[int]] & in_neigh_G,
                                       cpp_map[int, cpp_vector[int]] & in_neigh_H,
                                       cpp_map[int, cpp_vector[int]] & out_neigh_G,
                                       cpp_map[int, cpp_vector[int]] & out_neigh_H,
                                       cpp_map[cpp_pair[int, int], int] & edges_G,
                                       cpp_map[cpp_pair[int, int], int] & edges_H,
                                       cpp_map[int, int] & total_order,
                                       cpp_vector[cpp_vector[cpp_pair[int, int]]] & all_matches) noexcept:
    # local variables
    cdef size_t new_score = 0
    cdef size_t old_score = 0
    cdef cpp_bool semantic_feasibility = False
    cdef cpp_bool syntactic_feasibility = False
    cdef cpp_vector[int] in_ring_G
    cdef cpp_vector[int] in_ring_H
    cdef cpp_vector[int] out_ring_G
    cdef cpp_vector[int] out_ring_H
    cdef cpp_vector[int] current_match_G
    cdef cpp_vector[int] current_match_H
    cdef cpp_vector[cpp_pair[int, int]] new_match
    cdef cpp_vector[cpp_pair[int, int]] candidates
    cdef cpp_map[int, int] forward_match
    cdef cpp_pair[int, int] each_pair
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
            current_match_G.push_back(each_pair.first)
            current_match_H.push_back(each_pair.second)
            forward_match[each_pair.first] = each_pair.second
        # get candidate pairs
        candidates_tuple = directed_candidates(current_match,
                                               current_match_G,
                                               current_match_H,
                                               nodes_G,
                                               nodes_H,
                                               in_neigh_G,
                                               in_neigh_H,
                                               out_neigh_G,
                                               out_neigh_H,
                                               total_order)
        # unpack candidates information
        candidates = candidates_tuple.first
        ring_G = candidates_tuple.second.first
        ring_H = candidates_tuple.second.second
        # evaluate candidates
        for each_pair in candidates:
            # evaluate sintactic feasibility
            syntactic_feasibility = directed_syntactic_feasibility(each_pair.first,
                                                                   each_pair.second,
                                                                   ring_G,
                                                                   ring_H,
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
                                                                     ring_G,
                                                                     ring_H,
                                                                     current_match_G,
                                                                     current_match_H,
                                                                     forward_match,
                                                                     nodes_G,
                                                                     nodes_H,
                                                                     in_neigh_G,
                                                                     in_neigh_H,
                                                                     out_neigh_G,
                                                                     out_neigh_H,
                                                                     edges_G,
                                                                     edges_H)
                if(semantic_feasibility):
                    # build new match
                    new_match = current_match
                    new_match.push_back(each_pair)
                    # extend match
                    directed_search_isomorphisms(expected_order,
                                                 new_match,
                                                 all_isomorphisms,
                                                 nodes_G,
                                                 nodes_H,
                                                 in_neigh_G,
                                                 in_neigh_H,
                                                 out_neigh_G,
                                                 out_neigh_H,
                                                 edges_G,
                                                 edges_H,
                                                 total_order,
                                                 all_matches)
                    # finish if only one isosmorphism was requested and it was already found
                    if(not all_isomorphisms):
                        if(not all_matches.empty()):
                            return
    # end of function



# function: get candidate pairs for directed isomorphism search ----------------
cdef isomorphisms_candidates_struct_directed directed_candidates(cpp_vector[cpp_pair[int, int]] current_match,
                                                                 cpp_vector[int] & current_match_G,
                                                                 cpp_vector[int] & current_match_H,
                                                                 cpp_map[int, int] & nodes_G,
                                                                 cpp_map[int, int] & nodes_H,
                                                                 cpp_map[int, cpp_vector[int]] & in_neigh_G,
                                                                 cpp_map[int, cpp_vector[int]] & in_neigh_H,
                                                                 cpp_map[int, cpp_vector[int]] & out_neigh_G,
                                                                 cpp_map[int, cpp_vector[int]] & out_neigh_H,
                                                                 cpp_map[int, int] & total_order) noexcept:
    # local variables
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int reference_minimum = -1
    cdef cpp_vector[int] in_ring_G
    cdef cpp_vector[int] in_ring_H
    cdef cpp_vector[int] out_ring_G
    cdef cpp_vector[int] out_ring_H
    cdef cpp_vector[int] valid_G
    cdef cpp_vector[int] valid_H
    cdef cpp_vector[cpp_pair[int, int]] candidate_pairs
    cdef cpp_vector[cpp_pair[int, int]] filtered_candidate_pairs
    cdef cpp_pair[int, int] temp_pair
    cdef cpp_pair[int, int] each_pair
    cdef isomorphisms_candidates_struct_directed candidates_struct
    # get valid sets
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
            # if not yet in match
            if(find(current_match_H.begin(), current_match_H.end(), node) == current_match_H.end()):
                # if total order greater than in match
                if(total_order[node] > reference_maximum):
                    valid_H.push_back(node)
        # make product of valid in-neighbors
        if((not valid_G.empty()) and (not valid_H.empty())):
            for node1 in valid_G:
                for node2 in valid_H:
                    temp_pair.first = node1
                    temp_pair.second = node2
                    if(find(candidate_pairs.begin(), candidate_pairs.end(), temp_pair) == candidate_pairs.end()):
                        candidate_pairs.push_back(temp_pair)
        # additionally save in-neighbors in ring around match in G
        for node in valid_G:
            # save if unrepeated
            if(find(ring_G.begin(), ring_G.end(), node) == ring_G.end()):
                ring_G.push_back(node)
        # additionally save in-neighbors in ring around match in H
        for node in valid_H:
            # save if unrepeated
            if(find(ring_H.begin(), ring_H.end(), node) == ring_H.end()):
                ring_H.push_back(node)
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
        # additionally save neighbors in ring around match in G
        for node in valid_G:
            # save if unrepeated
            if(find(ring_G.begin(), ring_G.end(), node) == ring_G.end()):
                ring_G.push_back(node)
        # additionally save neighbors in ring around match in H
        for node in valid_H:
            # save if unrepeated
            if(find(ring_H.begin(), ring_H.end(), node) == ring_H.end()):
                ring_H.push_back(node)
    # alternatively get all unmatched nodes also not adjacent with match
    if(candidate_pairs.empty()):
        # reinitialize valid neighbors
        valid_G.clear()
        valid_H.clear()
        # get alternatives in G
        for each_pair in nodes_G:
            node = each_pair.first
            if(find(current_match_G.begin(), current_match_G.end(), node) == current_match_G.end()):
                if(find(ring_G.begin(), ring_G.end(), node) == ring_G.end()):
                    valid_G.push_back(node)
        # get alternatives in H
        for each_pair in nodes_H:
            node = each_pair.first
            if(find(current_match_H.begin(), current_match_H.end(), node) == current_match_H.end()):
                if(find(ring_H.begin(), ring_H.end(), node) == ring_H.end()):
                    valid_H.push_back(node)
        # make product of valid neighbors
        if((not valid_G.empty()) and (not valid_H.empty())):
            for node1 in valid_G:
                for node2 in valid_H:
                    temp_pair.first = node1
                    temp_pair.second = node2
                    # during isomorphism search the pair cannot be repeated at this point
                    candidate_pairs.push_back(temp_pair)
    # filter candidate pairs
    if(not candidate_pairs.empty()):
        # get minimum total order
        for each_pair in candidate_pairs:
            node = each_pair.second
            if(reference_minimum == -1):
                reference_minimum = total_order[node]
                filtered_candidate_pairs.push_back(each_pair)
            else:
                if(total_order[node] == reference_minimum):
                    filtered_candidate_pairs.push_back(each_pair)
                if(total_order[node] < reference_minimum):
                    reference_minimum = total_order[node]
                    filtered_candidate_pairs.clear()
                    filtered_candidate_pairs.push_back(each_pair)
        # copy just for clean output coding
        candidate_pairs = filtered_candidate_pairs
    # pack return structure
    candidates_struct.candidates = candidate_pairs
    candidates_struct.in_ring_G = in_ring_G
    candidates_struct.in_ring_H = in_ring_H
    candidates_struct.out_ring_G = out_ring_G
    candidates_struct.out_ring_H = out_ring_H
    # end of function
    return(candidates_struct)



# function: evaluate syntactic feasability for directed isomorphisms -----------
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



# function: evaluate semantic feasability for directed isomorphisms ------------
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

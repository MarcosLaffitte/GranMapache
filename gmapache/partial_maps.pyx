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
# - Note: currently broadly using (typed) python lists and dictionaries since  #
#   we need the dynamic allocation, but we want to slowly but surely migrate   #
#   into pure C and C++ structures and objects.                                #
#                                                                              #
################################################################################



# dependencies #################################################################



# already in python ------------------------------------------------------------
import time
from copy import deepcopy
from itertools import product
from sys import getrecursionlimit, setrecursionlimit



# not in python ----------------------------------------------------------------
import networkx as nx



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
    # local variables
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int counter = 0
    cdef int current_limit = 0
    cdef int expected_order = 0
    cdef int required_limit = 0
    cdef float scalation_value = 0
    cdef list encoded_graphs = []
    cdef list encoded_anchor = []
    cdef list inside_anchor_H = []
    cdef list outside_anchor_H = []
    cdef list encoded_extensions = []
    cdef dict info = dict()
    cdef dict nodes_G = dict()
    cdef dict edges_G = dict()
    cdef dict nodes_H = dict()
    cdef dict edges_H = dict()
    cdef dict total_order = dict()
    cdef dict encoded_node_names = dict()
    cdef dict encoded_node_label = dict()
    cdef dict encoded_edge_label = dict()
    cdef dict in_neigh_G = dict()
    cdef dict in_neigh_H = dict()
    cdef dict out_neigh_G = dict()
    cdef dict out_neigh_H = dict()
    cdef dict neigh_G = dict()
    cdef dict neigh_H = dict()
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
        encoded_extensions = directed_maximum_connected_extensions(nodes_G, edges_G, in_neigh_G, out_neigh_G,
                                                                   nodes_H, edges_H, in_neigh_H, out_neigh_H,
                                                                   expected_order,
                                                                   encoded_anchor,
                                                                   [],
                                                                   total_order)
    else:
        encoded_extensions = undirected_maximum_connected_extensions(nodes_G, edges_G, neigh_G,
                                                                     nodes_H, edges_H, neigh_H,
                                                                     expected_order,
                                                                     encoded_anchor,
                                                                     [],
                                                                     total_order)
    # decode maximum extensions
    for each_extension in encoded_extensions:
        extensions.append(decode_match(each_extension, encoded_node_names))
    # check if the anchor was a good partial map
    if((len(extensions[0]) == len(nodes_G)) and (len(extensions[0]) == len(nodes_H))):
        good_anchor = True
    # end of function
    return(extensions, good_anchor)



# functions - maximum connected extensions - undirected ########################



# function: core routine of VF2-like approach ----------------------------------
cdef undirected_maximum_connected_extensions(nodes_G = dict(), edges_G = dict(), neigh_G = dict(),
                                             nodes_H = dict(), edges_H = dict(), neigh_H = dict(),
                                             expected_order = 0,
                                             current_match = [],
                                             all_matches = [],
                                             total_order = dict()):
    # local variables
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int new_score = 0
    cdef int old_score = 0
    cdef list new_match = []
    cdef list candidates = []
    cdef list current_match_G = []
    cdef list current_match_H = []
    cdef dict forward_match = dict()
    cdef dict inverse_match = dict()
    syntactic_feasibility = False
    semantic_feasibility = False
    # test initial match and consecutive matches
    if(len(all_matches) == 0):
        if(len(current_match) > 0):
            all_matches = [current_match]
    else:
        # test improvement in matches
        new_score = len(current_match)
        old_score = len(all_matches[0])
        # save match if it has the same score
        if(new_score == old_score):
            all_matches.append(current_match)
        # overwrite with new match if it improves score
        if(new_score > old_score):
            all_matches = [current_match]
    # if not optimal yet then obtain available pairs
    if(len(current_match) < expected_order):
        # generate auxiliary structures
        current_match_G = [node1 for (node1, node2) in current_match]
        current_match_H = [node2 for (node1, node2) in current_match]
        forward_match = {node1:node2 for (node1, node2) in current_match}
        inverse_match = {node2:node1 for (node1, node2) in current_match}
        # get candidate pairs
        candidates = undirected_candidates(current_match_G, current_match_H,
                                           neigh_G, neigh_H,
                                           total_order)
        # evaluate candidates
        for (node1, node2) in candidates:
            # evaluate sintactic feasibility
            syntactic_feasibility = undirected_syntactic_feasibility(node1, node2,
                                                                     current_match_G, current_match_H,
                                                                     forward_match, inverse_match,
                                                                     neigh_G, neigh_H)
            if(syntactic_feasibility):
                # evaluate semantic feasibility
                semantic_feasibility = undirected_semantic_feasibility(node1, node2,
                                                                       current_match_G,
                                                                       forward_match,
                                                                       nodes_G, edges_G, neigh_G,
                                                                       nodes_H, edges_H)
                if(semantic_feasibility):
                    # extend match
                    new_match = current_match + [(node1, node2)]
                    all_matches = undirected_maximum_connected_extensions(nodes_G, edges_G, neigh_G,
                                                                          nodes_H, edges_H, neigh_H,
                                                                          expected_order,
                                                                          new_match,
                                                                          all_matches,
                                                                          total_order)
    # end of function
    return(all_matches)



# function: get candidate pairs for undir extension search ---------------------
cdef undirected_candidates(current_match_G, current_match_H,
                           neigh_G, neigh_H,
                           total_order):
    # local variables
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int reference_maximum = 0
    cdef list valid_G = []
    cdef list valid_H = []
    cdef list candidate_pairs = []
    # get maximum value of total order in match
    reference_maximum = max([total_order[node] for node in current_match_H])
    # get valid sets
    for node1 in current_match_G:
        for node2 in neigh_G[node1]:
            if(node2 not in valid_G):
                if(node2 not in current_match_G):
                    valid_G.append(node2)
    for node1 in current_match_H:
        for node2 in neigh_H[node1]:
            if(node2 not in valid_H):
                if(node2 not in current_match_H):
                    if(total_order[node2] > reference_maximum):
                        valid_H.append(node2)
    # get candidates
    candidate_pairs = list(product(valid_G, valid_H))
    # end of function
    return(candidate_pairs)



# function: evaluate syntactic feasability for undir extension -----------------
cdef undirected_syntactic_feasibility(node1, node2,
                                      current_match_G, current_match_H,
                                      forward_match, inverse_match,
                                      neigh_G, neigh_H):
    # local variables
    cdef int node = 0
    cdef list neighbors_match_G = []
    cdef list neighbors_match_H = []
    # loop-consistency-test
    if((node1 in neigh_G[node1]) and (node2 not in neigh_H[node2])):
        return(False)
    if((node1 not in neigh_G[node1]) and (node2 in neigh_H[node2])):
        return(False)
    # look ahead 0: consistency of neighbors in match
    neighbors_match_G = [node for node in neigh_G[node1] if(node in current_match_G)]
    neighbors_match_H = [node for node in neigh_H[node2] if(node in current_match_H)]
    for node in neighbors_match_G:
        if(forward_match[node] not in neighbors_match_H):
            return(False)
    for node in neighbors_match_H:
        if(inverse_match[node] not in neighbors_match_G):
            return(False)
    # end of function
    return(True)



# function: evaluate semantic feasability for undir extension ------------------
cdef undirected_semantic_feasibility(node1, node2,
                                     current_match_G,
                                     forward_match,
                                     nodes_G, edges_G, neigh_G,
                                     nodes_H, edges_H):
    # local variables
    cdef int a1 = 0
    cdef int a2 = 0
    cdef int b1 = 0
    cdef int b2 = 0
    cdef int node = 0
    cdef list neighbors_match_G = []
    # compare vertex-labels
    if(not nodes_G[node1] == nodes_H[node2]):
        return(False)
    # compare loop-labels
    if(node1 in neigh_G[node1]):
        if(not edges_G[(node1, node1)] == edges_H[(node2, node2)]):
            return(False)
    # compare non-loop edge-labels
    neighbors_match_G = [node for node in neigh_G[node1] if(node in current_match_G)]
    for node in neighbors_match_G:
        a1 = min([node1, node])
        b1 = max([node1, node])
        a2 = min([node2, forward_match[node]])
        b2 = max([node2, forward_match[node]])
        if(not edges_G[(a1, b1)] == edges_H[(a2, b2)]):
            return(False)
    # end of function
    return(True)



# functions - maximum connected extensions - directed ##########################



# function: core routine of VF2-like approach ----------------------------------
cdef directed_maximum_connected_extensions(nodes_G = dict(), edges_G = dict(), in_neigh_G = dict(), out_neigh_G = dict(),
                                           nodes_H = dict(), edges_H = dict(), in_neigh_H = dict(), out_neigh_H = dict(),
                                           expected_order = 0,
                                           current_match = [],
                                           all_matches = [],
                                           total_order = dict()):
    # local variables
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int new_score = 0
    cdef int old_score = 0
    cdef list new_match = []
    cdef list candidates = []
    cdef list current_match_G = []
    cdef list current_match_H = []
    cdef dict forward_match = dict()
    cdef dict inverse_match = dict()
    syntactic_feasibility = False
    semantic_feasibility = False
    # test initial match and consecutive matches
    if(len(all_matches) == 0):
        if(len(current_match) > 0):
            all_matches = [current_match]
    else:
        # test improvement in matches
        new_score = len(current_match)
        old_score = len(all_matches[0])
        # save match if it has the same score
        if(new_score == old_score):
            all_matches.append(current_match)
        # overwrite with new match if it improves score
        if(new_score > old_score):
            all_matches = [current_match]
    # if not optimal yet then obtain available pairs
    if(len(current_match) < expected_order):
        # generate auxiliary structures
        current_match_G = [node1 for (node1, node2) in current_match]
        current_match_H = [node2 for (node1, node2) in current_match]
        forward_match = {node1:node2 for (node1, node2) in current_match}
        inverse_match = {node2:node1 for (node1, node2) in current_match}
        # get candidate pairs
        candidates = directed_candidates(current_match_G, current_match_H,
                                         in_neigh_G, out_neigh_G,
                                         in_neigh_H, out_neigh_H,
                                         total_order)
        # evaluate candidates
        for (node1, node2) in candidates:
            # evaluate sintactic feasibility
            syntactic_feasibility = directed_syntactic_feasibility(node1, node2,
                                                                   current_match_G, current_match_H,
                                                                   forward_match, inverse_match,
                                                                   in_neigh_G, out_neigh_G,
                                                                   in_neigh_H, out_neigh_H)
            if(syntactic_feasibility):
                # evaluate semantic feasibility
                semantic_feasibility = directed_semantic_feasibility(node1, node2,
                                                                     current_match_G,
                                                                     forward_match,
                                                                     nodes_G, edges_G, in_neigh_G, out_neigh_G,
                                                                     nodes_H, edges_H)
                if(semantic_feasibility):
                    # extend match
                    new_match = current_match + [(node1, node2)]
                    all_matches = directed_maximum_connected_extensions(nodes_G, edges_G, in_neigh_G, out_neigh_G,
                                                                        nodes_H, edges_H, in_neigh_H, out_neigh_H,
                                                                        expected_order,
                                                                        new_match,
                                                                        all_matches,
                                                                        total_order)
    # end of function
    return(all_matches)



# function: get candidate pairs for undir extension search ---------------------
cdef directed_candidates(current_match_G, current_match_H,
                         in_neigh_G, out_neigh_G,
                         in_neigh_H, out_neigh_H,
                         total_order):
    # local variables
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int reference_maximum = 0
    cdef list valid_G = []
    cdef list valid_H = []
    cdef list candidate_pairs = []
    # get maximum value of total order in match
    reference_maximum = max([total_order[node] for node in current_match_H])
    # get valid sets
    for node1 in current_match_G:
        for node2 in list(set(in_neigh_G[node1] + out_neigh_G[node1])):
            if(node2 not in valid_G):
                if(node2 not in current_match_G):
                    valid_G.append(node2)
    for node1 in current_match_H:
        for node2 in list(set(in_neigh_H[node1] + out_neigh_H[node1])):
            if(node2 not in valid_H):
                if(node2 not in current_match_H):
                    if(total_order[node2] > reference_maximum):
                        valid_H.append(node2)
    # get candidates
    candidate_pairs = list(product(valid_G, valid_H))
    # end of function
    return(candidate_pairs)



# function: evaluate syntactic feasability for undir extension -----------------
cdef directed_syntactic_feasibility(node1, node2,
                                    current_match_G, current_match_H,
                                    forward_match, inverse_match,
                                    in_neigh_G, out_neigh_G,
                                    in_neigh_H, out_neigh_H):
    # local variables
    cdef int node = 0
    cdef list in_neighbors_match_G = []
    cdef list in_neighbors_match_H = []
    cdef list out_neighbors_match_G = []
    cdef list out_neighbors_match_H = []
    # loop-consistency-test
    if((node1 in in_neigh_G[node1]) and (node2 not in in_neigh_H[node2])):
        return(False)
    if((node1 not in in_neigh_G[node1]) and (node2 in in_neigh_H[node2])):
        return(False)
    # look ahead 0: consistency of neighbors in match
    in_neighbors_match_G = [node for node in in_neigh_G[node1] if(node in current_match_G)]
    in_neighbors_match_H = [node for node in in_neigh_H[node2] if(node in current_match_H)]
    for node in in_neighbors_match_G:
        if(forward_match[node] not in in_neighbors_match_H):
            return(False)
    for node in in_neighbors_match_H:
        if(inverse_match[node] not in in_neighbors_match_G):
            return(False)
    out_neighbors_match_G = [node for node in out_neigh_G[node1] if(node in current_match_G)]
    out_neighbors_match_H = [node for node in out_neigh_H[node2] if(node in current_match_H)]
    for node in out_neighbors_match_G:
        if(forward_match[node] not in out_neighbors_match_H):
            return(False)
    for node in out_neighbors_match_H:
        if(inverse_match[node] not in out_neighbors_match_G):
            return(False)
    # end of function
    return(True)



# function: evaluate semantic feasability for undir extension ------------------
cdef directed_semantic_feasibility(node1, node2,
                                   current_match_G,
                                   forward_match,
                                   nodes_G, edges_G, in_neigh_G, out_neigh_G,
                                   nodes_H, edges_H):
    # local variables
    cdef int a1 = 0
    cdef int a2 = 0
    cdef int b1 = 0
    cdef int b2 = 0
    cdef int node = 0
    cdef list in_neighbors_match_G = []
    cdef list out_neighbors_match_G = []
    # compare vertex-labels
    if(not nodes_G[node1] == nodes_H[node2]):
        return(False)
    # compare loop-labels
    if(node1 in in_neigh_G[node1]):
        if(not edges_G[(node1, node1)] == edges_H[(node2, node2)]):
            return(False)
    # compare non-loop edge-labels
    in_neighbors_match_G = [node for node in in_neigh_G[node1] if(node in current_match_G)]
    for node in in_neighbors_match_G:
        if(not edges_G[(node, node1)] == edges_H[(forward_match[node], node2)]):
            return(False)
    out_neighbors_match_G = [node for node in out_neigh_G[node1] if(node in current_match_G)]
    for node in out_neighbors_match_G:
        if(not edges_G[(node1, node)] == edges_H[(node2, forward_match[node])]):
            return(False)
    # end of function
    return(True)



################################################################################
################################################################################

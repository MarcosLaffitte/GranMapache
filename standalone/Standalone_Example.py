################################################################################
#                                                                              #
# - GranMapache: GRAphs-and-Networks MAPping Applications                      #
#   with Cython and HEuristics                                                 #
#                                                                              #
# - Description: suite for the analysis of bijections, morphisms, alignments   #
#   and other maps between graphs.                                             #
#                                                                              #
# - Repository (privite): https://github.com/MarcosLaffitte/GranMapache        #
#                                                                              #
# - Example of usage                                                           #
#                                                                              #
# - Instructions:                                                              #
#                                                                              #
#   * installing required python version into anaconda environment             #
#     conda create -n [env_name] python=3.11.0                                 #
#                                                                              #
#   * for all the following, activate anaconda environment with                #
#     conda activate [env_name]                                                #
#                                                                              #
#   * in order to install this package you need C and C++ compilers. If you    #
#     dont have them already in your system you can install them INSIDE the    #
#     anaconda environment with                                                #
#     conda install -c conda-forge cxx-compiler                                #
#                                                                              #
#   * installing the package and dependencies INSIDE anaconda environment,     #
#     go to the folder where setup.py is located and run                       #
#     python -m pip install .                                                  #
#                                                                              #
#   * running Standalone_Example.py (make sure to "activate conda [env_name]") #
#     python Standalone_Example.py                                             #
#                                                                              #
#   * uninstalling gmapache from anaconda environment                          #
#     python -m pip uninstall gmapache                                         #
#                                                                              #
#   * removing anaconda environmnet completely                                 #
#     conda remove -n [env_name] --all                                         #
#                                                                              #
# - NOTE:                                                                      #
#                                                                              #
#   * the module is not yet available online in the pip system, but that       #
#     is our intention for this package in the future, which will make the     #
#     the instalation / removal simpler.                                       #
#                                                                              #
#   * not installing the package inside an anaconda enviornment or python      #
#     environment may produce unexpected dependency erros, which my produce    #
#     incompatibility warnings and also runtime errors, therefore our          #
#     recommendation for installing this package inside an environment.        #
#                                                                              #
#   * when installing additional packages not required by gmapache the         #
#     installers may suggest different python versions, thus always check the  #
#     compatibility of those version with the required by gmapache.            #
#                                                                              #
#   * after installing gmapache you can call it from a script in any location  #
#     in your computer, except for a script inside the parent folder of the    #
#     gmapache package. This will change once we made the package public.      #
#                                                                              #
################################################################################


# dependencies #################################################################


# already in python ------------------------------------------------------------
import time


# not in python ----------------------------------------------------------------
import numpy
import cython
import networkx as nx
import matplotlib.pyplot as plt


# custom dependencies ----------------------------------------------------------
import gmapache as gm


# examples #####################################################################


# example: stable extension ----------------------------------------------------


# We build two networkx graphs G and H with arbitrary vertex labels and edge labels,
# and define a "partial map" as an injective match between the graphs, also refered
# to as "anchor" given as a list of unrepeated 2-tuples (x, y) of vertices x in G and
# y in H. This data is passed to the method gm.search_stable_extension(G, H, partial_map),
# which extends the partial map and returns two results: (1) a list of all the possible
# stable extensions of the reaction center also expressed as matches, and
# (2) a boolean value indicating if the extensions form bijections and, equivalentely,
# if the "anchor" was a good partial map.


# buid first test graph for stable extension
G = nx.Graph()
G.add_node("g1", color = "blue", charge = "0")
G.add_node("g2", color = "blue", charge = "0")
G.add_node("g3", color = "blue", charge = "0")
G.add_node("g4", color = "blue", charge = "0")
G.add_node("g5", color = "blue", charge = "0")
G.add_node("g6", color = "blue", charge = "0")
G.add_node("g7", color = "green", charge = "+1")
G.add_node("g8", color = "green", charge = "+1")
G.add_node("g9", color = "green", charge = "+1")
G.add_node("g10", color = "red", charge = "-1")
G.add_node("g11", color = "red", charge = "-1")
G.add_node("g12", color = "red", charge = "-1")
G.add_edge("g1", "g6", bond = "single", strength = "neutral")
G.add_edge("g2", "g3", bond = "single", strength = "neutral")
G.add_edge("g4", "g5", bond = "single", strength = "neutral")
G.add_edge("g1", "g7", bond = "double", strength = "strong")
G.add_edge("g7", "g8", bond = "double", strength = "strong")
G.add_edge("g7", "g9", bond = "double", strength = "strong")
G.add_edge("g4", "g10", bond = "triple", strength = "weak")
G.add_edge("g10", "g11", bond = "triple", strength = "weak")
G.add_edge("g11", "g12", bond = "triple", strength = "weak")


# buid second test graph for stable extension
H = nx.Graph()
H.add_node("h1", color = "blue", charge = "0")
H.add_node("h2", color = "blue", charge = "0")
H.add_node("h3", color = "blue", charge = "0")
H.add_node("h4", color = "blue", charge = "0")
H.add_node("h5", color = "blue", charge = "0")
H.add_node("h6", color = "blue", charge = "0")
H.add_node("h10", color = "green", charge = "+1")
H.add_node("h11", color = "green", charge = "+1")
H.add_node("h12", color = "green", charge = "+1")
H.add_node("h7", color = "red", charge = "-1")
H.add_node("h8", color = "red", charge = "-1")
H.add_node("h9", color = "red", charge = "-1")
H.add_edge("h1", "h2", bond = "single", strength = "neutral")
H.add_edge("h3", "h4", bond = "single", strength = "neutral")
H.add_edge("h5", "h6", bond = "single", strength = "neutral")
H.add_edge("h1", "h10", bond = "double", strength = "strong")
H.add_edge("h10", "h11", bond = "double", strength = "strong")
H.add_edge("h10", "h12", bond = "double", strength = "strong")
H.add_edge("h4", "h7", bond = "triple", strength = "weak")
H.add_edge("h7", "h8", bond = "triple", strength = "weak")
H.add_edge("h8", "h9", bond = "triple", strength = "weak")


# buid test anchor for stable extension
partial_map = [("g1", "h1"),
               ("g2", "h2"),
               ("g3", "h3"),
               ("g4", "h4"),
               ("g5", "h5"),
               ("g6", "h6")]


# testing stable extension
initial_time = time.time()
all_extensions, good_map = gm.search_stable_extension(G, H, partial_map)
final_time = time.time()


# print results
print("--------------------------------------------------")
print("> Stable Extension")
print("\n")
print("***** Given anchor:")
print(partial_map)
print("***** Got extensions:")
print(all_extensions)
print("***** Was anchor a good partial map?")
print(good_map)
print("***** Running time [s]")
print(final_time - initial_time)
print("\n")


# example: maximum common induced anchored subgraphs ---------------------------


# The function gm.search_maximum_common_anchored_subgraphs(G, H, input_anchor) receives
# two non-null networkx (di-)graphs G and H (possibly with different number of nodes), and a
# non-empty injective map between them (here called anchor), and searches for the maximum
# common induced subgraphs that extend the matches in the anchor. The function specifically
# searches for proper extensions of the anchor, i.e., containing at least one more match
# than the anchor itself. The parameter "reachability", controls if the candidate common
# subgraphs should be "connected" to the anchor, or better put "reachable" from the anchor
# in the sense that such subgraphs should contain a path between every node and at least one
# node from the anchor. Thus, if the anchor is connected and reachability is set to True,
# the function will equivalently get a maximum-common-induced-connected-anchored-subgraph.
# If both graphs have the same order, the function first will search for a complete-induced-
# extension and only if no such extension is found it will continue with the search for the
# maximum common induced anchored subgraphs, thus this function can be more time consuming
# that simply running the search for the stable extension. Moreover, if both graphs
# have the same order and a stable extension exists between them, this function
# will return such extension independently of the parameter "reachability" and regardless
# if the complete extension induces a connected ITS.


# first graph
G = nx.Graph()
G.add_node(1, color = "blue")
G.add_node(2, color = "red")
G.add_node(3, color = "red")
G.add_node(4, color = "orange")
G.add_node(5, color = "purple")
G.add_node(6, color = "pink")
G.add_edge(1, 2, weight = 1)
G.add_edge(1, 3, weight = 1)
G.add_edge(2, 4, weight = 1)
G.add_edge(3, 4, weight = 1)
G.add_edge(4, 5, weight = 1)
G.add_edge(5, 6, weight = 1)


# second graph
H = nx.Graph()
H.add_node("a", color = "blue")
H.add_node("b", color = "red")
H.add_node("c", color = "red")
H.add_node("d", color = "green")
H.add_node("e", color = "yellow")
H.add_node("f", color = "purple")
H.add_node("g", color = "pink")
H.add_edge("a", "b", weight = 1)
H.add_edge("a", "c", weight = 1)
H.add_edge("b", "d", weight = 1)
H.add_edge("c", "d", weight = 1)
H.add_edge("d", "e", weight = 1)
H.add_edge("d", "f", weight = 1)
H.add_edge("f", "g", weight = 1)


# partial map, with nodes x from first graph and y from second graph
initial_map = [(1, "a"), (2, "b")]


# get maximum common subgraph with reachability constraint
initial_time = time.time()
extensions, found_extensions = gm.search_maximum_common_anchored_subgraphs(G, H, input_anchor = initial_map,
                                                                           node_labels = True, edge_labels = True,
                                                                           all_extensions = False, reachability = True)
final_time = time.time()


# print resutls
print("--------------------------------------------------")
print("> Maximum Common Induced Subgraph (I - with reachability constraint)")
print("\n")
print("***** Received anchor")
print(initial_map)
print("***** Were any proper extensions found?")
print(found_extensions)
print("***** Anchor or Extensions:")
print(extensions)
print("***** Running time [s]")
print(final_time - initial_time)
print("\n")


# get maximum common subgraph without reachability constraint
initial_time = time.time()
extensions, found_extensions = gm.search_maximum_common_anchored_subgraphs(G, H, input_anchor = initial_map,
                                                                           node_labels = True, edge_labels = True,
                                                                           all_extensions = False, reachability = False)
final_time = time.time()


# print resutls
print("--------------------------------------------------")
print("> Maximum Common Induced Subgraph (II - without reachability constraint)")
print("\n")
print("***** Received anchor")
print(initial_map)
print("***** Were any proper extensions found?")
print(found_extensions)
print("***** Anchor or Extensions:")
print(extensions)
print("***** Running time [s]")
print(final_time - initial_time)
print("\n")


# example: isomorphisms --------------------------------------------------------


# We build four networkx graphs G, H, F and L. They are passed to the method
# gm.search_isomorphisms(G, H), which determines if there exist at least one
# isomorphism between these graphs and returns it. A boolean variable is also
# returned indicating if the graphs were indeed ismorphic. This method has an
# optional boolean argument "all_isomorphisms". When False the function returns
# only one isomorphism (if any), this is the default behavior since it is faster
# than enumerating all isomorphisms. If set to True the function will search and
# return all possible isomorphisms, which is an exhaustive procedure and thus can
# be more time consuming depending on the input graphs. Moreover, the function can
# take into acount all the node-labels and/or edge-labels of the graphs. This is
# controlled by the boolean variables node_labels and edge_labels which by default
# are set to False, thus ignoring labels. The function can optionally also take into
# account a total order for the nodes of the graphs. This is required by the all
# the VF2-like algorithms to carry on the search. If this is not provided then
# the functions computes a totoal order internally, based on the (out-) degrees
# of the nodes in the underlying (unlabled) versions of the inout graphs.


# build first test graph for isomorphism search
G = nx.Graph()
G.add_edge(1, 2)
G.add_edge(2, 3)
G.add_edge(3, 4)
G.add_edge(4, 5)
G.add_edge(5, 6)
G.add_edge(5, 7)


# build second test graph for isomorphism search
H = nx.Graph()
H.add_edge("a", "b")
H.add_edge("b", "c")
H.add_edge("c", "d")
H.add_edge("d", "e")
H.add_edge("e", "f")
H.add_edge("e", "g")


# build third test graph for isomorphism search
F = nx.Graph()
F.add_edge("x1", "x2")
F.add_edge("x1", "x3")
F.add_edge("x1", "x4")
F.add_edge("x1", "x5")
F.add_edge("x1", "x6")
F.add_edge("x1", "x7")


# build fourth test graph for isomorphism search
L = nx.Graph()
L.add_node(1, color = "blue")
L.add_node(2, color = "blue")
L.add_node(3, color = "red")
L.add_node(4, color = "red")
L.add_node(5, color = "green")
L.add_node(6, color = "purple")
L.add_node(7, color = "purple")
L.add_edge(1, 2, bond = "single")
L.add_edge(2, 3, bond = "double")
L.add_edge(3, 4, bond = "single")
L.add_edge(4, 5, bond = "double")
L.add_edge(5, 6, bond = "single")
L.add_edge(5, 7, bond = "triple")


# testing isomorphisms
initial_time = time.time()
isomorphisms, are_isomorphic = gm.search_isomorphisms(G, H, all_isomorphisms = True)
final_time = time.time()


# print resutls
print("--------------------------------------------------")
print("> Isomorphism Search (underlying isomorphic)")
print("\n")
print("***** Got isomorphisms:")
print(isomorphisms)
print("***** Were graphs isomorphic?")
print(are_isomorphic)
print("***** Running time [s]")
print(final_time - initial_time)
print("\n")


# testing isomorphisms
initial_time = time.time()
isomorphisms, are_isomorphic = gm.search_isomorphisms(G, F)
final_time = time.time()


# print resutls
print("--------------------------------------------------")
print("> Isomorphism Search (underlying non-isomorphic)")
print("\n")
print("***** Got isomorphisms:")
print(isomorphisms)
print("***** Were graphs isomorphic?")
print(are_isomorphic)
print("***** Running time [s]")
print(final_time - initial_time)
print("\n")


# testing isomorphisms
initial_time = time.time()
isomorphisms, are_isomorphic = gm.search_isomorphisms(G, L, all_isomorphisms = True)
final_time = time.time()


# print resutls
print("--------------------------------------------------")
print("> Isomorphism Search (underlying isomorphic and ignoring labels)")
print("\n")
print("***** Got isomorphisms:")
print(isomorphisms)
print("***** Were graphs isomorphic?")
print(are_isomorphic)
print("***** Running time [s]")
print(final_time - initial_time)
print("\n")


# testing isomorphisms
initial_time = time.time()
isomorphisms, are_isomorphic = gm.search_isomorphisms(G, L, node_labels = True, edge_labels = True, all_isomorphisms = True)
final_time = time.time()


# print resutls
print("--------------------------------------------------")
print("> Isomorphism Search (underlying isomorphic but not isomorphic when considering labels)")
print("\n")
print("***** Got isomorphisms:")
print(isomorphisms)
print("***** Were graphs isomorphic?")
print(are_isomorphic)
print("***** Running time [s]")
print(final_time - initial_time)
print("\n")


# testing isomorphisms
initial_time = time.time()
isomorphisms, are_isomorphic = gm.search_isomorphisms(L, L, node_labels = True, edge_labels = True, all_isomorphisms = True)
final_time = time.time()


# print resutls
print("--------------------------------------------------")
print("> Isomorphism Search (A - underlying have two isomorphisms but with labels only one because of edge labels)")
print("\n")
print("***** Got isomorphisms:")
print(isomorphisms)
print("***** Were graphs isomorphic?")
print(are_isomorphic)
print("***** Running time [s]")
print(final_time - initial_time)
print("\n")


# testing isomorphisms
initial_time = time.time()
isomorphisms, are_isomorphic = gm.search_isomorphisms(L, L, node_labels = True, edge_labels = False, all_isomorphisms = True)
final_time = time.time()


# print resutls
print("--------------------------------------------------")
print("> Isomorphism Search (B - underlying have two isomorphisms but with labels only one because of edge labels; ignoring here only edge labels)")
print("\n")
print("***** Got isomorphisms:")
print(isomorphisms)
print("***** Were graphs isomorphic?")
print(are_isomorphic)
print("***** Running time [s]")
print(final_time - initial_time)
print("\n")


# example: subgraph isomorphisms -----------------------------------------------


# We build two labeled graphs G and H, where G is isomorphic to a subgraph of H.
# Here we consider both induced and general subgraphs. In other words, we can
# search for subgraph isomorphisms from G to arbitrary subgraphs of H (also
# called monomorphisms in CS), or for subgraph isomorphisms from G to induced
# subgraphs of H (subgraph isomorphisms in CS). We call the function
# gm.search_subgraph_isomorphisms(G, H) to search for embeddings of G in H.
# The function by default searches for unlabeld subgraph isomorphisms, but similar
# to the other functions in this package it can receive arguments to explicitly
# use node labels, edge labels, or both. Similarly it searches by default for only
# one embedding, but an argument can be passed to obtain all possible embeddings
# (if any). Moreover, this function searches for subgraph isomorphisms of G with
# induced subgraphs of H by default, but an argument can be passed to search for
# monomorphisms instead, i.e., subgraph isomorphisms of G with arbitrary subgraphs
# of H. If both graphs have the same number of vertices and edges then the function
# runs a graph isomorphism search instead of the subgraph search.


# build (smaller) first graph for subgraph isomorphism test
G = nx.Graph()
# add blue square with edges of weight one
G.add_node(1, color = "blue")
G.add_node(2, color = "blue")
G.add_node(3, color = "blue")
G.add_node(4, color = "blue")
G.add_edge(1, 2, weight = 1)
G.add_edge(2, 3, weight = 1)
G.add_edge(3, 4, weight = 1)
G.add_edge(4, 1, weight = 1)
# add pendant node
G.add_node(5, color = "root")
G.add_edge(1, 5, weight = 0)


# build (bigger) second graph for subgraph isomorphism test
H = nx.Graph()
# add blue square with edges of weight one
H.add_node("a", color = "blue")
H.add_node("b", color = "blue")
H.add_node("c", color = "blue")
H.add_node("d", color = "blue")
H.add_edge("a", "b", weight = 1)
H.add_edge("b", "c", weight = 1)
H.add_edge("c", "d", weight = 1)
H.add_edge("d", "a", weight = 1)
# add blue square with edges of weight 2.5
H.add_node("e", color = "blue")
H.add_node("f", color = "blue")
H.add_node("g", color = "blue")
H.add_node("h", color = "blue")
H.add_edge("e", "f", weight = 2.5)
H.add_edge("f", "g", weight = 2.5)
H.add_edge("g", "h", weight = 2.5)
H.add_edge("h", "e", weight = 2.5)
# add green square with edges of weight 2.5
H.add_node("w", color = "green")
H.add_node("x", color = "green")
H.add_node("y", color = "green")
H.add_node("z", color = "green")
H.add_edge("w", "x", weight = 2.5)
H.add_edge("x", "y", weight = 2.5)
H.add_edge("y", "z", weight = 2.5)
H.add_edge("z", "w", weight = 2.5)
# connect the three squares to a single node
H.add_node("center", color = "root")
H.add_edge("a", "center", weight = 0)
H.add_edge("e", "center", weight = 0)
H.add_edge("w", "center", weight = 0)


# build square graph for non-induced subgraph isomorphism test
S = nx.Graph()
# add blue square with edges of weight one
S.add_node(1, color = "blue")
S.add_node(2, color = "blue")
S.add_node(3, color = "blue")
S.add_node(4, color = "blue")
S.add_edge(1, 2, weight = 1)
S.add_edge(2, 3, weight = 1)
S.add_edge(3, 4, weight = 1)
S.add_edge(4, 1, weight = 1)


# build kite graph for non-induced subgraph isomorphism test
K = nx.Graph()
# add blue square with edges of weight one
K.add_node("a", color = "blue")
K.add_node("b", color = "blue")
K.add_node("c", color = "blue")
K.add_node("d", color = "blue")
K.add_node("e", color = "red")
K.add_edge("a", "b", weight = 1)
K.add_edge("b", "c", weight = 1)
K.add_edge("c", "d", weight = 1)
K.add_edge("d", "a", weight = 1)
K.add_edge("b", "d", weight = 1)
K.add_edge("c", "e", weight = 1)


# testing subgraph isomorphisms
initial_time = time.time()
embeddings, subgraph_isomorphic = gm.search_subgraph_isomorphisms(G, H, node_labels = True, edge_labels = True, all_isomorphisms = False)
final_time = time.time()


# print resutls
print("--------------------------------------------------")
print("> Subgraph Isomorphism Search (A - considering node and edge labels; only one embedding)")
print("\n")
print("***** Got subgraph isomorphisms:")
print(embeddings)
print("***** Were the graphs subgraph-isomorphic?")
print(subgraph_isomorphic)
print("***** Running time [s]")
print(final_time - initial_time)
print("\n")


# testing subgraph isomorphisms
initial_time = time.time()
embeddings, subgraph_isomorphic = gm.search_subgraph_isomorphisms(G, H, node_labels = True, edge_labels = True, all_isomorphisms = True)
final_time = time.time()


# print resutls
print("--------------------------------------------------")
print("> Subgraph Isomorphism Search (B - considering node and edge labels; all embeddings)")
print("\n")
print("***** Got subgraph isomorphisms:")
print(embeddings)
print("***** Were the graphs subgraph-isomorphic?")
print(subgraph_isomorphic)
print("***** Running time [s]")
print(final_time - initial_time)
print("\n")


# testing subgraph isomorphisms
initial_time = time.time()
embeddings, subgraph_isomorphic = gm.search_subgraph_isomorphisms(G, H, node_labels = True, edge_labels = False, all_isomorphisms = True)
final_time = time.time()


# print resutls
print("--------------------------------------------------")
print("> Subgraph Isomorphism Search (C - considering only node labels; all embeddings)")
print("\n")
print("***** Got subgraph isomorphisms:")
print(embeddings)
print("***** Were the graphs subgraph-isomorphic?")
print(subgraph_isomorphic)
print("***** Running time [s]")
print(final_time - initial_time)
print("\n")


# testing subgraph isomorphisms
initial_time = time.time()
embeddings, subgraph_isomorphic = gm.search_subgraph_isomorphisms(G, H, node_labels = False, edge_labels = False, all_isomorphisms = True)
final_time = time.time()


# print resutls
print("--------------------------------------------------")
print("> Subgraph Isomorphism Search (D - no labels; all embeddings)")
print("\n")
print("***** Got subgraph isomorphisms:")
print(embeddings)
print("***** Were the graphs subgraph-isomorphic?")
print(subgraph_isomorphic)
print("***** Running time [s]")
print(final_time - initial_time)
print("\n")


# testing subgraph isomorphisms
initial_time = time.time()
embeddings, subgraph_isomorphic = gm.search_subgraph_isomorphisms(G, H, node_labels = False, edge_labels = False, all_isomorphisms = False)
final_time = time.time()


# print resutls
print("--------------------------------------------------")
print("> Subgraph Isomorphism Search (E - no labels; only one embedding)")
print("\n")
print("***** Got subgraph isomorphisms:")
print(embeddings)
print("***** Were the graphs subgraph-isomorphic?")
print(subgraph_isomorphic)
print("***** Running time [s]")
print(final_time - initial_time)
print("\n")


# testing non-induced subgraph isomorphisms
initial_time = time.time()
embeddings, subgraph_isomorphic = gm.search_subgraph_isomorphisms(S, K, node_labels = True, edge_labels = True, all_isomorphisms = True)
final_time = time.time()


# print resutls
print("--------------------------------------------------")
print("> Subgraph Isomorphism Search (test for induced - default option; all isomorphisms)")
print("\n")
print("***** Got subgraph isomorphisms:")
print(embeddings)
print("***** Were the graphs induced subgraph-isomorphic?")
print(subgraph_isomorphic)
print("***** Running time [s]")
print(final_time - initial_time)
print("\n")


# testing non-induced subgraph isomorphisms
initial_time = time.time()
embeddings, subgraph_isomorphic = gm.search_subgraph_isomorphisms(S, K, node_labels = True, edge_labels = True, all_isomorphisms = True, induced = False)
final_time = time.time()


# print resutls
print("--------------------------------------------------")
print("> Subgraph Isomorphism Search (test for non-induced; all isomorphisms)")
print("\n")
print("***** Got subgraph isomorphisms:")
print(embeddings)
print("***** Were the graphs monomorphic?")
print(subgraph_isomorphic)
print("***** Running time [s]")
print(final_time - initial_time)
print("\n")


################################################################################
################################################################################

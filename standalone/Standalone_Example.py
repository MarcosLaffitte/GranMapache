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
#     conda create -n [env_name] python=3.11.10                                #
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
#     pip install .                                                            #
#                                                                              #
#   * running Standalone_Example.py (make sure to "activate conda [env_name]") #
#     python Standalone_Example.py                                             #
#                                                                              #
#   * uninstalling gmapache from anaconda environment                          #
#     pip uninstall gmapache                                                   #
#                                                                              #
#   * removing anaconda environmnet completely                                 #
#     conda remove -n [env_name] --all                                         #
#                                                                              #
# - NOTE:                                                                      #
#                                                                              #
#   * the module is not yet available online in the pip system, but that       #
#     is our intention for this package in the future.                         #
#                                                                              #
#   * not installing the package inside an anaconda enviornment or python      #
#     environment may produce unexpected dependency erros, which my produce    #
#     incompatibility warnings and also runtime errors, therefore our          #
#     recommendation for installing this package inside an environment.        #
#                                                                              #
#   * when installing additional packages not required by gmapache, the        #
#     installers may suggest different python versions, thus always check the  #
#     compatibility of those version with the required by gmapache             #
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


# example: induced connected extensions ----------------------------------------


# We build two networkx graphs G and H with arbitrary vertex labels and edge labels, and define
# a "partial map" as an injective match between the graphs, also refered to as "anchor" given
# as a list of unrepeated 2-tuples (x, y) of vertices x in G and y in H. This data is passed to the
# method gm.induced_connected_extensions(G, H, partial_map), which extends the partial map
# and returns two results: (1) a list of all the possible induced extensions of the reaction center
# also expressed as matches, and (2) a boolean value indicating if the extensions form bijections
# and, equivalentely, if the "anchor" was a good partial map.


# buid first test graph for induced connected extension
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


# buid second test graph for induced connected extension
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


# buid test anchor for induced connected extension
partial_map = [("g1", "h1"),
               ("g2", "h2"),
               ("g3", "h3"),
               ("g4", "h4"),
               ("g5", "h5"),
               ("g6", "h6")]


# testing induced connected extension
initial_time = time.time()
all_extensions, good_map = gm.induced_connected_extensions(G, H, partial_map)
final_time = time.time()


# print results
print("--------------------------------------------------")
print("> Induced Connected Extension")
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


# example: isomorphisms --------------------------------------------------------


# We build four networkx graphs G, H, F and L. They are passed to the method
# gm.search_isomorphisms(G, H), which determines if there exist at least one
# isomorphism between these graphs and returns it. A boolean variable is also
# returned indicating if the graphs were indeed ismorphic. This method has an
# optional boolean argument "all_isomorphisms". When False the function returns
# only one isomorphism (if any), this is the default behavior since it is faster
# than enumerating all isomorphism. If set to True the function will search and
# return all possible isomorphisms, which is an exhaustive and can be more time
# consuming depending on the input graphs. Moreover, the function can take into
# acount all the node-labels and/or edge-labels of the graphs. This is controlled
# by the booleand variables node_labels and edge_labels which by default are set
# to False thus ignoring labels. The function can optionally also take into account
# a total order for the nodes of the graphs. This is required by the all the VF2-like
# algorithms to carry on the search. If this is not provided then the functions
# computes a totoal order internally, based on the degrees of the nodes in the
# underlying (unlabled) versions of the inout graphs.


# buid first test graph for isomorphism search
G = nx.Graph()
G.add_edge(1, 2)
G.add_edge(2, 3)
G.add_edge(3, 4)
G.add_edge(4, 5)
G.add_edge(5, 6)
G.add_edge(5, 7)


# buid second test graph for isomorphism search
H = nx.Graph()
H.add_edge("a", "b")
H.add_edge("b", "c")
H.add_edge("c", "d")
H.add_edge("d", "e")
H.add_edge("e", "f")
H.add_edge("e", "g")


# buid third test graph for isomorphism search
F = nx.Graph()
F.add_edge("x1", "x2")
F.add_edge("x1", "x3")
F.add_edge("x1", "x4")
F.add_edge("x1", "x5")
F.add_edge("x1", "x6")
F.add_edge("x1", "x7")


# buid fourth test graph for isomorphism search
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


################################################################################
################################################################################

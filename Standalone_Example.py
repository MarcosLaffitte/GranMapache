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
#   * builing cython code in place (if dependencies are already satisfied)     #
#     python setup.py build_ext --inplace                                      #
#                                                                              #
#   * running main of the module if it was built in place                      #
#     python -m gmapache                                                       #
#                                                                              #
#   * installing the package and dependencies inside anaconda environment      #
#     got to the folder where setup.py is located and run                      #
#     pip install .                                                            #
#                                                                              #
#   * running Standalone_Example.py if the package was installed in anaconda   #
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
#     environment may poduce dependency erros, both incompatibilities but      #
#     also runtime errors.                                                     #
#                                                                              #
################################################################################


# dependencies #################################################################


# already in python ------------------------------------------------------------
import time
from copy import deepcopy


# not in python ----------------------------------------------------------------
import numpy
import cython
import networkx as nx
import matplotlib.pyplot as plt


# custom dependencies ----------------------------------------------------------
import gmapache as gm


# example #########################################################################


# We build two networkx graphs G and H with arbitrary vertex labels and edge labels, and define
# a "reaction center" as an injective match between the graphs, also refered to as "anchor", given
# as a list of unrepeated 2-tuples (x, y) of vertices x in G and y in H. This data is passed to the
# method gm.partial_maps.maximum_connected_extensions(G, H, reaction_center), which extends the
# reaction center and returns two results: (1) a list of all the possible maximum extensions of the
# reaction center also expressed as matches, and (2) a boolean value indicating if the extension for
# a bijection and equivalentely if the "anchor" was a good partial map.


# buid first test graph for maximum connected extension
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


# buid second test graph for maximum connected extension
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


# buid test anchor for maximum connected extension
reaction_center = [("g1", "h1"),
                   ("g2", "h2"),
                   ("g3", "h3"),
                   ("g4", "h4"),
                   ("g5", "h5"),
                   ("g6", "h6")]


# testing maximum connected extension
all_extensions, good_center = gm.partial_maps.maximum_connected_extensions(G, H, reaction_center)
print("*** Given anchor:")
print(reaction_center)
print("*** Got extensions:")
print(all_extensions)
print("*** Was anchor a good reaction center?")
print(good_center)


################################################################################
################################################################################

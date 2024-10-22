################################################################################
#                                                                              #
# - GranMapache: GRAphs-and-Networks MAPping Applications                      #
#   with Cython and HEuristics                                                 #
#                                                                              #
# - Description: suite for the analysis of bijections, morphisms, alignments   #
#   and other maps between graphs.                                             #
#                                                                              #
# - Example of usage                                                           #
#                                                                              #
################################################################################


# dependencies #################################################################


# already in python ------------------------------------------------------------
import time
from copy import deepcopy


# not in python ----------------------------------------------------------------
import cython
import networkx as nx
import matplotlib.pyplot as plt
import numpy


# custom dependencies ----------------------------------------------------------
from .integerization import encode_graphs, decode_graphs, encode_match, decode_match
from .partial_maps import maximum_connected_extensions
# from .common_subgraphs import *
# from .verbosity import *


# functions ####################################################################


# main #########################################################################


if __name__ == "__main__":


    # buid first test graph
    G = nx.Graph()
    G.add_node("hola1", color = "blue", charge = -1, radius = 0.5)
    G.add_node("hola2", color = "red", charge = 1, radius = 0.75)
    G.add_node("hola3", color = "red", charge = 1, radius = 0.9)
    G.add_node("hola4", color = "blue", charge = -1, radius = 0.5)
    G.add_node("hola5", color = "green", charge = -1, radius = 0.25)
    G.add_edge("hola1", "hola2", weight = 10, bond = "single", strength = "weak")
    G.add_edge("hola2", "hola3", weight = 8, bond = "double", strength = "strong")
    G.add_edge("hola3", "hola4", weight = 1, bond = "triple", strength = "strong")
    G.add_edge("hola4", "hola5", weight = 10, bond = "single", strength = "weak")


    # buid second test graph
    H = nx.Graph()
    H.add_node("hola1", color = "red", charge = 1, radius = 0.5)
    H.add_node("hola6", color = "red", charge = 1, radius = 0.75)
    H.add_node("hola7", color = "red", charge = 1, radius = 0.9)
    H.add_node("hola8", color = "blue", charge = -1, radius = 0.5)
    H.add_node("hola9", color = "green", charge = -1, radius = 0.25)
    H.add_edge("hola1", "hola6", weight = 10, bond = "single", strength = "weak")
    H.add_edge("hola1", "hola7", weight = 7, bond = "triple", strength = "strong")
    H.add_edge("hola1", "hola8", weight = 10, bond = "triple", strength = "strong")
    H.add_edge("hola1", "hola9", weight = 3, bond = "single", strength = "strong")


    # build input graph
    input_list = [G, H]


    # testing encoding
    graph_encoding = encode_graphs(input_list)
    encoded_graph_list = graph_encoding[0]
    encoded_node_names = graph_encoding[1]
    encoded_node_label = graph_encoding[2]
    encoded_edge_label = graph_encoding[3]


    # testing decoding
    recovered_graphs = decode_graphs(encoded_graph_list,
                                     encoded_node_names,
                                     encoded_node_label,
                                     encoded_edge_label)

    # encode a match
    match = [("hola1", "hola1"),
             ("hola2", "hola6"),
             ("hola3", "hola7"),
             ("hola4", "hola8")]
    match_encoding = encode_match(match, encoded_node_names)


    # decode a match
    recovered_match = decode_match(match_encoding, encoded_node_names)


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
    G.add_edge("g9", "g10", bond = "triple", strength = "weak")
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
    print(maximum_connected_extensions(G, H, reaction_center))





################################################################################
################################################################################

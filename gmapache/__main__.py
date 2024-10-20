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
# - Observations and ToDos:                                                    #
#   * for a large number of graphs, integerization with cython is kinda faster #
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
from .partial_maps import maximum_connected_extension
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
    encoded_node_name = graph_encoding[1]
    encoded_node_label = graph_encoding[2]
    encoded_edge_label = graph_encoding[3]


    # testing decoding
    recovered_graphs = decode_graphs(encoded_graph_list,
                                     encoded_node_name,
                                     encoded_node_label,
                                     encoded_edge_label)

    # encode a match
    match = [("hola1", "hola1"),
             ("hola2", "hola6"),
             ("hola3", "hola7"),
             ("hola4", "hola8")]
    match_encoding = encode_match(match, encoded_node_name)


    # decode a match
    recovered_match = decode_match(match_encoding, encoded_node_name)


    # testing maximum connected extension
    maximum_connected_extension()


################################################################################
################################################################################

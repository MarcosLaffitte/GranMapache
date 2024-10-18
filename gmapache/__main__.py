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


# not in python ----------------------------------------------------------------
import cython
import networkx as nx


# custom dependencies ----------------------------------------------------------
from .integerization import encode_graphs
# from .common_subgraphs import bla
# from .verbosity import bli


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

    # testing encoding
    encoded_graphs = encode_graphs([G, H])



################################################################################
################################################################################

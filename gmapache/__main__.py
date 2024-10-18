################################################################################
#                                                                              #
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
    
    # testing encoding
    G = nx.path_graph(5)
    the_list = encode_graphs([G])
    
    

################################################################################
################################################################################

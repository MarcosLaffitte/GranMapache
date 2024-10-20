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
################################################################################


# dependencies #################################################################


# already in python ------------------------------------------------------------
from copy import deepcopy


# not in python ----------------------------------------------------------------
import cython
import networkx as nx


# numpy in cython --------------------------------------------------------------
cimport numpy as cnp


# custom dependencies ----------------------------------------------------------
from .integerization import encode_graphs, decode_graphs, encode_match, decode_match


# algorithms ###################################################################


# functions - maximum connected extension  #####################################


# function: callable wrapper for the maximum connected extension ---------------
def maximum_connected_extension():
    # description
    """
    > description: receives two graphs G, and H, and a match between them, here called
    anchor, and uses a VF2-like approach to obtain a maximum extension of the anchor
    producing a connected common subgraph (not necessarily maximum itslef). The anchor
    alone also produces a subgraph, which may not be an induced common subgraph, but
    the subgraph produced by the extension after removing the achor is always induced.

    > input:
    * G - first networkx (di)graph benig matched
    * H - second networkx (di)graph benig matched
    * anchor - list of 2-tuples (x, y) of nodes x from G and y from H.

    > output:
    * extension - list of 2-tuples (x, y) of nodes x from G and y from H representing
    the maximum connected extension of the anchor.
    * good_extension - boolean indicating if the extension covers now all nodes, i.e.,
    it is a bijection, between G and H.

    > calls:
    * .integerization.encode_graphs
    * .integerization.decode_graphs
    * .integerization.encode_match
    * .integerization.decode_match
    """
    # output holders
    # cython variables
    # local variables
    # decode match
    # end of function
    return(0)


################################################################################
################################################################################

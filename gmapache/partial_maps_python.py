################################################################################
#                                                                              #
# - GranMapache: GRAphs-and-Networks MAPping Applications                      #
#   with Cython and HEuristics                                                 #
#                                                                              #
# - Module: integerization                                                     #
#                                                                              #
# - Description: convert node "names", node attributes and edge attributes     #
#   into integers to simplify their comparison in the other rutines            #
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


# functions ####################################################################


# function: --------------------------------------------------------------------


################################################################################
################################################################################

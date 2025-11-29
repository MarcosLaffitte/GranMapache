################################################################################
#                                                                              #
#  - Generation of Random Hydrogen Extensions                                  #
#                                                                              #
#  - for AMB - BMC - extension of WABI 2025                                    #
#                                                                              #
#  - Made by Marcos Laffitte - Github @MarcosLaffitte                          #
#                                                                              #
################################################################################




# dependencies #################################################################




# already in python ------------------------------------------------------------
import sys
import math
import time
import random
import pickle
from copy import deepcopy
from itertools import permutations
from operator import eq, itemgetter




# additional dependencies ------------------------------------------------------
import numpy
import cython
import networkx as nx
import pysmiles as ps
import cmasher as cmr
from rdkit import Chem
import matplotlib.pyplot as plt




# networkx isomorphism with node and edge match --------------------------------
from operator import eq
from networkx.algorithms import isomorphism
from networkx.algorithms.isomorphism import generic_node_match, generic_edge_match




# long's synkit ----------------------------------------------------------------
from synkit.IO import rsmi_to_its, rsmi_to_graph, load_from_pickle
from synkit.Graph.Hyrogen._misc import check_hcount_change
from synkit.Graph.ITS.its_decompose import get_rc, its_decompose




# custom dependencies ----------------------------------------------------------
import gmapache as gm




# warnings ---------------------------------------------------------------------
import logging
logging.getLogger("pysmiles").setLevel(logging.CRITICAL)




# parameters ###################################################################




# parameters -------------------------------------------------------------------
reactions = 100
carbons_G = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
carbons_H = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
hydrogens = [2, 3, 4, 5]




# analysis #####################################################################




# output data holders
random_hydrogen_extensions = []




# traversal of hydrogen amounts
for H_num in hydrogens:

    # generation of reactions
    for i in range(reactions):

        # choose G-side adjacent carbons
        carbon_neighbors_G = random.sample(carbons_G, H_num)

        # choose H-side adjacent carbons
        carbon_neighbors_H = random.sample(carbons_H, H_num)

        # save reaction
        random_hydrogen_extensions.append((deepcopy(carbon_neighbors_G), deepcopy(carbon_neighbors_H)))




# pickle data
output_file = open("random_reactions.pkl", "wb")
pickle.dump(random_hydrogen_extensions, output_file)
output_file.close()




################################################################################
################################################################################
